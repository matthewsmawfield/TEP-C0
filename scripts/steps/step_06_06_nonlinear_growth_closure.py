#!/usr/bin/env python3
"""Step 06-06: TEP Nonlinear Growth Closure (Halo-Model, Gradient-Screened Benchmark)

Evaluates the gradient-dependent screening envelope f(g) = [1 + (g/g_t)^n]^-1
on the EdS acoustic-sector benchmark using the characteristic Newtonian
acceleration at the virial radius of each halo (canonical threshold
g_t = c H_0 / 2 ~= 3.4e-10 m/s^2 under strictly constant M_Pl in S_grav;
the canonical two-branch scale-transition and window-weighted sigma_8
decomposition are evaluated in Paper 18 Step 21 Part D).

Background status: the EdS acoustic-sector benchmark used here is the
counterfactual matter-only branch (Step 05-11). The corpus's physical
single-history background is the native Jordan-frame branch
H = M(z) H_LCDM, which fits Pantheon+, DESI DR2 BAO and the Planck acoustic
angle jointly (Step 03-14); this step remains a screening-envelope benchmark
and does not define the physical background.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from c0_common import (
    read_json,
    TEPLogger,
    ensure_dirs,
    print_status,
    rounded,
    set_step_logger,
    step_csv_path,
    step_json_path,
    write_csv,
    write_json,
)

STEP_ID = "step_06_06_nonlinear_growth_closure"

# ============================================================================
# Cosmological parameters (EdS acoustic-sector benchmark)
# ============================================================================
H0 = 70.0                           # km/s/Mpc
OMEGA_M = 1.0
OMEGA_L = 0.0
SIGMA8_PLANCK = 0.812
SIGMA8_PLANCK_ERR = 0.007

# Halo-model parameters
RHO_CRIT = 2.775e11                 # h² M_sun / Mpc³ (critical density at z=0)
DELTA_VIR = 200.0                   # Overdensity for virialized halos


def halo_mass_function_dn_dM(M: float, z: float = 0.0) -> float:
    """Assigned power-law mass-function surrogate in physical-mass coordinates.

    Returns the approximate number density of halos per unit mass.
    """
    # For simplicity, use a power-law approximation with exponential cutoff
    # dn/dM ∝ M^(-1.9) exp(-M/M_*) where M_* is the characteristic mass
    M_star = 1e13  # physical M_sun, assigned characteristic mass at z=0
    norm = 1e-15
    return norm * (M / M_star) ** (-1.9) * np.exp(-M / M_star)


def halo_bias_ST(M: float, z: float = 0.0) -> float:
    """Sheth-Tormen halo bias (simplified).

    For EdS universe at z=0, the bias scales roughly as:
        b(M) ≈ 1 + (ν - 1) / δ_c
    where ν = δ_c² / σ²(M) and δ_c ≈ 1.686.
    """
    delta_c = 1.686
    # σ(M) ∝ M^(-1/6) for a power-law spectrum with n_s ≈ 1
    sigma_M = 5.0 * (M / 1e12) ** (-1.0 / 6.0)
    nu = (delta_c / sigma_M) ** 2
    bias = 1.0 + (nu - 1.0) / delta_c
    return max(bias, 0.5)


# ============================================================================
# Gradient-dependent screening (TEP v3)
# ============================================================================
from c0_screening_constants import SCREENING_ACCELERATION_M_S2
G_T = SCREENING_ACCELERATION_M_S2  # master c H0/(2 beta_A^2), m/s^2
N_SCREEN = 2.0


def gradient_screening_envelope(g: float, g_t: float = G_T, n: float = N_SCREEN) -> float:
    """TEP gradient-dependent screening envelope f(g) = [1 + (g/g_t)^n]^-1."""
    ratio = g / g_t
    return 1.0 / (1.0 + ratio ** n)


def halo_characteristic_acceleration(M_msun: float, z: float = 0.0,
                                     delta_vir: float = 200.0) -> float:
    """Characteristic Newtonian acceleration at the virial radius of a halo.

    g_vir = G * M / R_vir^2  [m/s^2]
    """
    # Critical density at z=0 in h^2 M_sun / Mpc^3
    rho_crit_0 = 2.775e11
    Ez2 = (1.0 + z) ** 3
    rho_crit_z = rho_crit_0 * (H0 / 100.0) ** 2 * Ez2
    rho_vir = delta_vir * rho_crit_z
    R_vir = (3.0 * M_msun / (4.0 * np.pi * rho_vir)) ** (1.0 / 3.0)  # physical Mpc

    G = 6.67430e-11
    M_sun_kg = 1.98847e30
    Mpc_m = 3.08567758e22

    M_kg = M_msun * M_sun_kg
    R_m = R_vir * Mpc_m
    g_vir = G * M_kg / (R_m ** 2)
    return g_vir


def screened_halo_bias(M: float, z: float = 0.0) -> float:
    """Halo bias with gradient-dependent TEP screening.

    For typical halos g_vir << g_t, so f(g) ≈ 1 (unscreened).
    Only the very central regions of the most massive clusters approach
    g_t, giving marginal suppression.
    """
    b_lin = halo_bias_ST(M, z)
    g_halo = halo_characteristic_acceleration(M, z)
    f_g = gradient_screening_envelope(g_halo)
    # In halos where g < g_t, the coupling is active, preserving
    # fifth-force enhancement of structure formation.
    return 1.0 + (b_lin - 1.0) * f_g


def compute_screened_sigma8(sigma8_linear: float | None = None) -> dict:
    """Compute σ₈ from the screened halo-model power spectrum.

    Uses a simplified two-halo model where the nonlinear power is
    suppressed on small scales due to halo-scale screening.
    """
    if sigma8_linear is None:
        upstream = read_json(step_json_path("step_06_04_growth_validation"))
        sigma8_linear = float(upstream["sigma_8"]["tep_eds"])
    if not np.isfinite(sigma8_linear) or sigma8_linear <= 0:
        raise ValueError("Current EdS diagnostic sigma8 must be finite and positive")

    # Integration limits (M_sun)
    M_min = 1e10
    M_max = 1e16

    # Split at the declared bias-floor kink for controlled quadrature.
    log_bias_floor = np.log(1e12 * ((1 - .5 * 1.686) * 25 / 1.686 ** 2) ** 3)

    # 1. Compute the unscreened 2-halo amplitude
    # I_b = ∫ dn/dM b(M) dM
    def integrand_bias(M):
        return halo_mass_function_dn_dM(M) * halo_bias_ST(M)

    I_b_unscreened, err_unscreened = quad(
        lambda logM: integrand_bias(np.exp(logM)) * np.exp(logM),
        np.log(M_min), np.log(M_max), epsabs=1e-12, epsrel=1e-10, limit=100, points=[log_bias_floor])

    # 2. Compute the screened 2-halo amplitude
    def integrand_bias_screened(M):
        return halo_mass_function_dn_dM(M) * screened_halo_bias(M)

    I_b_screened, err_screened = quad(
        lambda logM: integrand_bias_screened(np.exp(logM)) * np.exp(logM),
        np.log(M_min), np.log(M_max), epsabs=1e-12, epsrel=1e-10, limit=100, points=[log_bias_floor])
    if not (I_b_unscreened > 0 and I_b_screened > 0):
        raise ValueError("Positive halo integrands require positive integrals")
    if max(err_unscreened / I_b_unscreened, err_screened / I_b_screened) > 1e-8:
        raise ValueError("Halo-mass quadrature did not converge")

    # 3. The ratio of screened to unscreened bias integrals gives
    #    the suppression of the 2-halo term.
    bias_suppression = I_b_screened / I_b_unscreened

    # 4. In the halo model, σ₈² ∝ ∫ P(k) W²(kR) k² dk / (2π²)
    #    where W is the top-hat window.  The 2-halo term dominates
    #    on large scales (k < 1 h/Mpc), while the 1-halo term
    #    dominates on small scales.  For σ₈ (R = 8 Mpc/h), the
    #    2-halo term is dominant.
    #
    #    A simple approximation: σ₈_screened ≈ σ₈_linear × bias_suppression
    sigma8_screened = sigma8_linear * bias_suppression

    # 5. Compare to Planck
    deviation = abs(sigma8_screened - SIGMA8_PLANCK) / SIGMA8_PLANCK_ERR

    return {
        "I_b_unscreened": rounded(I_b_unscreened, 6),
        "I_b_screened": rounded(I_b_screened, 6),
        "bias_suppression_factor": rounded(bias_suppression, 6),
        "integration_coordinate": "ln(physical M/Msun), with Jacobian M",
        "quadrature_relative_error_max": max(err_unscreened / I_b_unscreened, err_screened / I_b_screened),
        "sigma8_linear": sigma8_linear,
        "sigma8_screened_halo_model": rounded(sigma8_screened, 6),
        "sigma8_planck": SIGMA8_PLANCK,
        "sigma8_planck_err": SIGMA8_PLANCK_ERR,
        "deviation_from_planck_sigma": rounded(deviation, 2),
        "passes_planck": deviation < 3,
    }


def run():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    print_status("Computing halo-model screened σ₈", "PROCESS")
    result = compute_screened_sigma8()

    print_status(f"  Linear σ₈ (EdS):            {result['sigma8_linear']:.4f}", "INFO")
    print_status(f"  Halo-model screened σ₈:     {result['sigma8_screened_halo_model']:.4f}", "INFO")
    print_status(f"  Planck σ₈:                  {result['sigma8_planck']:.4f} ± {result['sigma8_planck_err']:.4f}", "INFO")
    print_status(f"  Bias suppression factor:    {result['bias_suppression_factor']:.4f}", "INFO")
    print_status(f"  Deviation from Planck:      {result['deviation_from_planck_sigma']:.2f}σ", "INFO")

    if result["passes_planck"]:
        print_status("  Halo-model σ₈ is within 3σ of Planck", "SUCCESS")
    else:
        print_status("  Halo-model σ₈ still in tension with Planck", "WARNING")
        print_status("  Full nonlinear TEP-CLASS closure still required", "INFO")

    # ------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------
    # Characteristic acceleration for a reference halo
    g_ref = halo_characteristic_acceleration(1e13)
    g_cluster = halo_characteristic_acceleration(1e15)
    f_ref = gradient_screening_envelope(g_ref)
    f_cluster = gradient_screening_envelope(g_cluster)

    payload = {
        "step": STEP_ID,
        "status": "completed",
        "description": (
            "Gradient-dependent screened halo-model: f(g) replaces S(ρ). "
            "Cosmic halos have g << g_t, so f ≈ 1 and the halo model is "
            "essentially unscreened.  The ~0.55 growth suppression comes from "
            "α_M running in hi_class, not from environmental gradient screening."
        ),
        "method": "gradient_screened_halo_model_v3",
        "input_source": "step_06_04_growth_validation.json:sigma_8.tep_eds",
        "scope": "Simplified assigned mass-function/bias surrogate for the EdS diagnostic; not a native nonlinear perturbation prediction",
        "screening": {
            "operator": "f(g) = [1 + (g/g_t)^n]^-1",
            "g_t_m_s2": G_T,
            "n": N_SCREEN,
            "g_ref_1e13_msun_m_s2": f"{g_ref:.3e}",
            "g_cluster_1e15_msun_m_s2": f"{g_cluster:.3e}",
            "f_ref": rounded(f_ref, 6),
            "f_cluster": rounded(f_cluster, 6),
            "note": "g_vir << g_t for all cosmological halos; f ≈ 1 (unscreened)",
        },
        "parameters": {
            "H0": H0,
            "Omega_m": OMEGA_M,
            "delta_vir": DELTA_VIR,
            "rho_crit_h2_Msun_Mpc3": RHO_CRIT,
        },
        "result": result,
        "comparison_to_phenomenological": {
            "phenomenological_factor": 0.55,
            "gradient_screened_halo_model_suppression": result["bias_suppression_factor"],
            "phenomenological_sigma8": rounded(result["sigma8_linear"] * 0.55, 4),
            "gradient_screened_halo_model_sigma8": result["sigma8_screened_halo_model"],
            "note": "Gradient screening gives f≈1; suppression is from α_M running in full hi_class",
        },
        "conclusion": (
            "gradient_screened_unsuppressed_halo_model"
            if result["bias_suppression_factor"] > 0.95
            else "conditional_pass_halo_model"
        ),
    }

    write_json(step_json_path(STEP_ID), payload)
    write_csv(step_csv_path(STEP_ID), [{
        "quantity": "sigma8_screened",
        "value": result["sigma8_screened_halo_model"],
        "unit": "dimensionless",
        "note": "conditional_simplified_halo_model",
    }])

    print_status(f"Results written to {step_json_path(STEP_ID)}", "SUCCESS")
    print_status(f"{STEP_ID} complete", "TITLE")
    return payload


if __name__ == "__main__":
    run()
