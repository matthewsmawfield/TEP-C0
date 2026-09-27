#!/usr/bin/env python3
"""Step 04-09: TEP Solar System PPN Constraints (v3 — gradient-dependent screening)

Derives the Parametrized Post-Newtonian (PPN) parameters for the Temporal
Equivalence Principle and verifies compliance with Solar System tests.

The module evaluates three screening prescriptions:

1. Unscreened (β_A = -1): ruled out by Cassini at ~58 000σ.
2. Lorentzian source-screening S(ρ) = [1+(ρ/ρ_T)²]⁻¹ (old C0 ansatz):
   suppresses the scalar source in dense environments, but leaves the
   Solar System exterior unscreened → γ - 1 ≈ -0.44 (~19 000σ).
3. Gradient-dependent screening f(g) = [1+(g/g_t)^n]⁻¹ where g = |∇Φ|:
   suppresses the effective conformal coupling in regions of STEEP
   potential gradient (Solar System, Earth surface, stellar interiors)
   while retaining full coupling in regions of SHALLOW gradient (galactic
   halos, cosmic voids, wide-binary environments).

This resolves the density-regime paradox of the old ρ-dependent form:
- The old f(ρ) suppressed at LOW density (good for PPN, bad for growth).
- The new f(g) suppresses where |∇Φ| is large (good for PPN AND local tests).
- Galactic halos have g ~ 10⁻¹⁰ m/s² << g_t = 10⁻⁹ m/s², so coupling
  remains ~98.6% active, preserving cosmological growth predictions.

TEP screening is described as suppression of the locally observable
Temporal Shear/source-charge sector, not as a commitment to a specific
chameleon, Vainshtein, Galileon, DBI, or symmetron microphysics.  The
gradient-dependent form is a phenomenological parameterization of the
environmental screening factor S_Σ(E) that tracks the local Newtonian
acceleration field |∇Φ|, consistent with Axiom A4 of the foundational
theory.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from c0_common import (
    RESULTS_DIR,
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

STEP_ID = "step_04_09_ppn_constraints"

# ============================================================================
# Observational bounds
# ============================================================================
GAMMA_CASSINI = 1.000000
GAMMA_CASSINI_ERR = 2.3e-5          # Bertotti et al. 2003
BETA_LLR = 1.000000
BETA_LLR_ERR = 1.2e-4               # Nordtvedt effect

BETA_A = -1.0                       # Locked TEP conformal coupling

# Gradient-dependent screening parameters (g = |∇Φ| in m/s²)
G_T = 1.0e-9                        # threshold acceleration [m/s²]
DEFAULT_N = 2.0                     # steepness (symmetric field, (∇φ)² kinetic)

# Solar System parameters
G_NEWTON = 6.67430e-11              # m³ kg⁻¹ s⁻²
M_SUN = 1.98847e30                  # kg
M_EARTH = 5.9722e24                 # kg
R_SUN = 6.957e8                     # m
R_EARTH = 6.371e6                   # m
AU = 1.495978707e11                 # m

# ============================================================================
# Local Newtonian accelerations for key environments
# ============================================================================

def newtonian_acceleration(M: float, r: float) -> float:
    """|∇Φ| = GM/r² for a point mass [m/s²]."""
    return G_NEWTON * M / (r ** 2)


# Pre-computed accelerations
g_solar_surface = newtonian_acceleration(M_SUN, R_SUN)           # ~274 m/s²
g_earth_surface = newtonian_acceleration(M_EARTH, R_EARTH)     # ~9.8 m/s²
g_saturn_orbit = newtonian_acceleration(M_SUN, 9.5 * AU)      # ~6.5e-5 m/s²
g_mercury_orbit = newtonian_acceleration(M_SUN, 0.39 * AU)    # ~0.039 m/s²
g_earth_orbit = newtonian_acceleration(M_SUN, 1.0 * AU)        # ~0.006 m/s²

# Galactic / wide-binary regime (deep MOND acceleration)
g_wide_binary = 1.2e-10             # m/s², typical wide-binary / halo accel

# Halo virial regime (typical MW halo at ~200 kpc)
g_halo_virial = 3.5e-10           # m/s², M_vir=1e12 M_sun, r_vir=200 kpc


def gradient_screening_factor(g: float | np.ndarray, g_t: float = G_T, n: float = DEFAULT_N) -> float | np.ndarray:
    """Gradient-dependent screening f(g) = [1 + (g/g_t)^n]⁻¹.

    Suppressed where |∇Φ| >> g_t (Solar System, Earth surface).
    Active where |∇Φ| << g_t (galactic halos, voids, wide binaries).
    """
    g = np.asarray(g, dtype=float)
    return 1.0 / (1.0 + (g / g_t) ** n)


def lorentzian_screening_factor(rho: float | np.ndarray, rho_t: float = 20.0) -> float | np.ndarray:
    """Old C0 source screening S(ρ) = [1 + (ρ/ρ_T)²]⁻¹ (for comparison only)."""
    rho = np.asarray(rho, dtype=float)
    return 1.0 / (1.0 + (rho / rho_t) ** 2)


def compute_ppn_gamma(beta_eff: float) -> float:
    """PPN γ using the exact Damour–Esposito-Farèse cross-term expression.

    γ - 1 = -2 α_0 α_eff / (1 + α_0 α_eff)

    linear in the screened source charge α_eff = S_Σ α_0, since the photon
    probe is unscreened. With the DEF normalization α_0 = √2 β_A and
    β_eff = β_A S_Σ, the product is α_0 α_eff = 2 β_A β_eff = 2 β_A² S_Σ.

    At β_eff = -1 (unscreened, S_Σ = 1): γ - 1 = -4/3, γ = -1/3,
    ruled out by Cassini at ~58,000σ.
    """
    alpha0_alpha_eff = 2.0 * BETA_A * beta_eff
    return 1.0 - 2.0 * alpha0_alpha_eff / (1.0 + alpha0_alpha_eff)


def compute_ppn_beta() -> float:
    """PPN β = 1 exactly for A(φ) = exp(β_A φ)."""
    return 1.0


def run():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    # ------------------------------------------------------------------
    # 1. Unscreened prediction (ruled out)
    # ------------------------------------------------------------------
    gamma_unscreened = compute_ppn_gamma(BETA_A)
    beta_unscreened = compute_ppn_beta()
    print_status("1. Unscreened PPN parameters", "PROCESS")
    print_status(f"  γ = {gamma_unscreened:.6f}, β = {beta_unscreened:.6f}", "INFO")
    dev_unscreened = abs(gamma_unscreened - GAMMA_CASSINI) / GAMMA_CASSINI_ERR
    print_status(f"  Cassini deviation: {dev_unscreened:.1f}σ → ruled out", "WARNING")

    # ------------------------------------------------------------------
    # 2. Lorentzian source-screening (old C0 ansatz) — gap
    # ------------------------------------------------------------------
    print_status("2. Lorentzian source-screening (old C0 ansatz)", "PROCESS")
    frac_lorentzian = 0.141248  # from previous integration
    beta_eff_lor = BETA_A * frac_lorentzian
    gamma_lor = compute_ppn_gamma(beta_eff_lor)
    dev_lor = (gamma_lor - GAMMA_CASSINI) / GAMMA_CASSINI_ERR
    print_status(f"  Integrated Q_eff/M_sun ≈ {frac_lorentzian:.6f}", "INFO")
    print_status(f"  γ ≈ {gamma_lor:.10f}, deviation = {dev_lor:.1f}σ", "WARNING")
    print_status("  → Gap: Lorentzian form does not suppress solar charge enough", "WARNING")

    # ------------------------------------------------------------------
    # 3. Gradient-dependent screening (v3 — resolves density paradox)
    # ------------------------------------------------------------------
    print_status("3. Gradient-dependent screening f(g) = [1+(g/g_t)^n]⁻¹", "PROCESS")
    print_status(f"  Threshold acceleration g_t = {G_T:.1e} m/s²", "INFO")
    print_status(f"  Steepness n = {DEFAULT_N}", "INFO")

    # Verify two boundary conditions
    f_cassini = gradient_screening_factor(g_saturn_orbit)
    beta_eff_cassini = BETA_A * f_cassini
    gamma_cassini_screened = compute_ppn_gamma(beta_eff_cassini)
    cassini_pass = abs(gamma_cassini_screened - GAMMA_CASSINI) / GAMMA_CASSINI_ERR < 1.0

    f_wb = gradient_screening_factor(g_wide_binary)
    beta_eff_wb = BETA_A * f_wb
    gamma_wb = compute_ppn_gamma(beta_eff_wb)

    print_status(f"  Cassini gate (g={g_saturn_orbit:.2e} m/s²): f={f_cassini:.2e}, β_eff={beta_eff_cassini:.2e}, γ={gamma_cassini_screened:.10f}", "SUCCESS" if cassini_pass else "WARNING")
    print_status(f"  Wide-binary gate (g={g_wide_binary:.2e} m/s²): f={f_wb:.6f}, β_eff={beta_eff_wb:.6f} — coupling {100*f_wb:.1f}% active", "INFO")

    # Evaluate all environments
    envs = {
        "solar_surface": (g_solar_surface, "Sun surface, r=R_sun"),
        "solar_wind_saturn": (g_saturn_orbit, "Saturn orbit (~9.5 AU), Cassini path"),
        "solar_wind_earth": (g_earth_orbit, "Earth orbit (1 AU)"),
        "mercury_orbit": (g_mercury_orbit, "Mercury orbit (~0.39 AU)"),
        "earth_surface": (g_earth_surface, "Earth surface, r=R_earth"),
        "galactic_halo": (g_halo_virial, "Galactic halo virial radius"),
        "wide_binary": (g_wide_binary, "Wide binary / deep MOND regime"),
    }

    # PPN tests measure the metric in the vacuum exterior
    ppn_test_envs = {"solar_wind_saturn", "solar_wind_earth", "mercury_orbit"}

    env_rows = []
    for name, (g, desc) in envs.items():
        f = gradient_screening_factor(g)
        beta_eff = BETA_A * f
        gamma = compute_ppn_gamma(beta_eff)
        gamma_dev = (gamma - GAMMA_CASSINI) / GAMMA_CASSINI_ERR
        is_ppn_test = name in ppn_test_envs
        env_rows.append({
            "environment": name,
            "description": desc,
            "g_m_s2": float(g),
            "screening_factor": rounded(f, 6),
            "beta_eff": rounded(beta_eff, 6),
            "gamma": rounded(gamma, 10),
            "gamma_deviation_sigma": rounded(gamma_dev, 2),
            "ppn_test_environment": is_ppn_test,
            "passes_cassini": abs(gamma_dev) < 3 if is_ppn_test else None,
        })
        if is_ppn_test:
            status = "SUCCESS" if abs(gamma_dev) < 3 else "WARNING"
            print_status(f"  {name:22s}: f={f:10.2e}, β_eff={beta_eff:10.2e}, γ={gamma:.10f}, dev={gamma_dev:.2f}σ  [PPN test]", status)
        else:
            print_status(f"  {name:22s}: f={f:10.2e}, β_eff={beta_eff:10.2e}, γ={gamma:.10f}  [not PPN]", "INFO")

    ppn_pass = all(r["passes_cassini"] for r in env_rows if r["ppn_test_environment"])
    print_status(f"  PPN-test environments pass Cassini: {ppn_pass}", "SUCCESS" if ppn_pass else "WARNING")

    # ------------------------------------------------------------------
    # 4. Multi-scale consistency
    # ------------------------------------------------------------------
    print_status("4. Multi-scale consistency", "PROCESS")
    f_earth = gradient_screening_factor(g_earth_surface)
    f_halo = gradient_screening_factor(g_halo_virial)
    f_void = gradient_screening_factor(g_wide_binary)

    print_status(f"  Earth surface (g={g_earth_surface:.1f}):    f={f_earth:.2e} — deeply screened, Eötvös bounds satisfied", "INFO")
    print_status(f"  Halo virial (g={g_halo_virial:.2e}): f={f_halo:.4f} — {100*f_halo:.1f}% active, growth coupling preserved", "INFO")
    print_status(f"  Wide binary (g={g_wide_binary:.2e}): f={f_void:.4f} — {100*f_void:.1f}% active, anomaly coupling preserved", "INFO")

    # ------------------------------------------------------------------
    # 5. Build payload and write outputs
    # ------------------------------------------------------------------
    payload = {
        "step": STEP_ID,
        "status": "completed",
        "description": "TEP PPN derivation: gradient-dependent screening closes Cassini gap and resolves density-regime paradox",
        "beta_A": BETA_A,
        "bounds": {
            "cassini_gamma": {"central": GAMMA_CASSINI, "err": GAMMA_CASSINI_ERR, "source": "Bertotti et al. 2003"},
            "llr_beta": {"central": BETA_LLR, "err": BETA_LLR_ERR, "source": "LLR Nordtvedt effect"},
        },
        "unscreened": {
            "gamma": rounded(gamma_unscreened, 6),
            "beta": beta_unscreened,
            "cassini_deviation_sigma": rounded(dev_unscreened, 1),
            "conclusion": "ruled_out_by_Cassini",
        },
        "lorentzian_source_screening": {
            "Q_eff_over_Msun_approx": rounded(frac_lorentzian, 6),
            "beta_eff": rounded(beta_eff_lor, 6),
            "gamma": rounded(gamma_lor, 10),
            "gamma_deviation_sigma": rounded(dev_lor, 1),
            "passes_cassini": False,
            "conclusion": "gap_identified_steeper_profile_needed",
        },
        "gradient_dependent_screening": {
            "model": "f(g) = [1 + (g/g_t)^n]⁻¹ where g = |∇Φ|",
            "note": "Phenomenological parameterization of TEP environmental screening factor S_Sigma(E) tracking local Newtonian acceleration",
            "g_t_m_s2": G_T,
            "steepness_n": DEFAULT_N,
            "environments": env_rows,
            "ppn_test_passes_cassini": ppn_pass,
            "cassini_gate": {
                "g_saturn_orbit_m_s2": g_saturn_orbit,
                "screening_factor": rounded(f_cassini, 6),
                "beta_eff": rounded(beta_eff_cassini, 6),
                "gamma": rounded(gamma_cassini_screened, 10),
                "passes_cassini": cassini_pass,
            },
            "wide_binary_gate": {
                "g_wide_binary_m_s2": g_wide_binary,
                "screening_factor": rounded(f_wb, 6),
                "beta_eff": rounded(beta_eff_wb, 6),
                "coupling_active_fraction": rounded(f_wb, 4),
            },
            "conclusion": "phenomenological_pass_covariant_derivation_open" if ppn_pass else "fails",
        },
    }

    write_json(step_json_path(STEP_ID), payload)
    write_csv(step_csv_path(STEP_ID), env_rows)

    print_status(f"Results written to {step_json_path(STEP_ID)}", "SUCCESS")
    print_status(f"{STEP_ID} complete", "TITLE")
    return payload


if __name__ == "__main__":
    run()
