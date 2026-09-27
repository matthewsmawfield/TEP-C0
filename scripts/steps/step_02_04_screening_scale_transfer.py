#!/usr/bin/env python3
"""TEP-C0 Step 02-04: Screening Scale Transfer — two-sector evaluation
=====================================================================

Evaluates the single covariant shear-sector operator

    S_Sigma(E) = [1 + (sqrt(Sigma^2)/g_t)^n + (rho_amb/rho_half)^2]^-1

and the clock-amplitude-sector geometric saturation factor

    S_geom = (rho_bar / rho_T)^(1/3) = R_T(M) / R_phys

on the same corpus environments, demonstrating that the two scales
rho_T = 20 g/cm^3 (material saturation, S_A sector) and
rho_half = 0.5 M_sun/pc^3 (ambient half-suppression, S_Sigma density
channel) parameterize different response projections of one scalar
configuration — not two values of a single screening variable.

Key points computed here:

1. The density argument rho in S_Sigma is the *ambient* coarse-grained
   environmental density of the region (galactic habitat ~0.5 M_sun/pc^3,
   near-Earth interplanetary medium ~1e-23 g/cm^3), not the material
   density of an embedded body.  Evaluating the ambient channel at a
   material density is the cross-domain substitution excluded by the
   corpus projection ontology (Paper 0 Section 7); the out-of-domain
   evaluation is computed explicitly and flagged.

2. At Earth's surface the correct environmental inputs give
   S_Sigma ~ 1e-21 (shear sector screened, as Cassini/Eoetvoes/LLR
   require) while the clock-amplitude sector retains S_geom ~ 0.65 —
   the transitional regime in which the GNSS clock-covariance signal
   (S_A channel, Papers 1/6/14/33) lives.  The two coexist because they
   are different functionals of the same field configuration.

3. The quartic density-dependent-mass completion V = lambda phi^4/4
   (candidate realization, Paper 0 Section 2.2) gives m_eff ∝ rho^(1/3),
   i.e. Compton length lambda_C ∝ rho^(-1/3): a single density-dependent
   mass spanning the ~24 orders of magnitude between ambient galactic
   and solid-matter densities.

4. R_T(M_earth) = (3 M_earth / 4 pi rho_T)^(1/3) ~ 4140 km matches the
   GNSS covariance scale lambda_T = 3330-4549 km (Papers 1/6/14).

Outputs:
- environments            : per-environment S_Sigma with correct inputs
- geometric_S_A_factors   : per-body S_geom = (rho_bar/rho_T)^(1/3)
- category_error_check    : ambient channel at material density (out of domain)
- compton_bridge          : lambda_C ∝ rho^(-1/3) across the density range
- R_T_earth_km            : geometric saturation radius vs measured lambda_T
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np

from c0_common import (
    RESULTS_DIR,
    TEPLogger,
    ensure_dirs,
    print_status,
    set_step_logger,
    step_json_path,
    write_json,
)

STEP_ID = "step_02_04_screening_scale_transfer"

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
G_SI = 6.67430e-11      # m^3 kg^-1 s^-2
C_SI = 2.99792458e8     # m/s
M_SUN_KG = 1.98847e30   # kg
M_SUN_G = 1.98847e33    # g
PC_CM = 3.086e18        # cm
AU_M = 1.495978707e11   # m
BETA_A = -1.0           # locked TEP conformal coupling

# ---------------------------------------------------------------------------
# The two scales — different sectors, one scalar configuration
# ---------------------------------------------------------------------------
RHO_T = 20.0                                   # g/cm^3  material saturation (S_A sector)
RHO_HALF_MSUN_PC3 = 0.5                        # M_sun/pc^3  ambient half-suppression (S_Sigma)
RHO_HALF_G_CM3 = RHO_HALF_MSUN_PC3 * M_SUN_G / PC_CM**3   # ~3.39e-23 g/cm^3
A_T = 3.4e-10                                  # m/s^2  shear acceleration scale ~ cH0/2
N_SCREEN = 2.0                                 # transition steepness

# Cassini requirement on the single-body solar source charge (Paper 0 Section 7)
S_SIGMA_CASSINI_MAX = 5.8e-6

# Measured GNSS covariance scale (Papers 1/14, Paper 6 calibration)
LAMBDA_T_MEASURED_KM = (3330.0, 4549.0)


def s_sigma(rho_amb_g_cm3: float, g_m_s2: float, n: float = N_SCREEN) -> dict:
    """Covariant shear-sector factor S_Sigma(E) with correct environmental inputs.

    In the Newtonian limit Sigma^2 ~ (beta_A g / c^2)^2 and g_t = a_t/c^2, so
    sqrt(Sigma^2)/g_t = |beta_A| g / a_t.  rho_amb is the *ambient*
    coarse-grained density of the region (not the body's material density).
    """
    shear_ratio = abs(BETA_A) * g_m_s2 / A_T
    dens_ratio = rho_amb_g_cm3 / RHO_HALF_G_CM3
    kinetic_term = shear_ratio ** n
    density_term = dens_ratio ** 2
    s = 1.0 / (1.0 + kinetic_term + density_term)
    return {
        "rho_ambient_g_cm3": rho_amb_g_cm3,
        "g_m_s2": g_m_s2,
        "shear_ratio_g_over_a_t": shear_ratio,
        "kinetic_term": kinetic_term,
        "density_term": density_term,
        "S_Sigma": s,
    }


def s_geom(rho_bar_g_cm3: float) -> float:
    """Clock-amplitude-sector geometric saturation factor S = (rho_bar/rho_T)^(1/3).

    Equals R_T(M)/R_phys for a uniform body; S < 1 places the surface on the
    unscreened side of the transition (S_A channel active).
    """
    return (rho_bar_g_cm3 / RHO_T) ** (1.0 / 3.0)


def r_t_km(mass_kg: float) -> float:
    """Geometric saturation radius R_T(M) = (3M / 4 pi rho_T)^(1/3) in km."""
    mass_g = mass_kg * 1.0e3
    return (3.0 * mass_g / (4.0 * np.pi * RHO_T)) ** (1.0 / 3.0) / 1.0e5


def run() -> dict:
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    print_status(f"Material saturation scale rho_T  = {RHO_T} g/cm^3  (S_A geometric sector)", "INFO")
    print_status(f"Ambient half-suppression rho_half = {RHO_HALF_G_CM3:.3e} g/cm^3  (S_Sigma density channel)", "INFO")
    print_status(f"Shear acceleration scale a_t      = {A_T:.2e} m/s^2  (~cH0/2)", "INFO")
    print_status(f"Hierarchy rho_T / rho_half        = {RHO_T / RHO_HALF_G_CM3:.2e}", "INFO")

    # ------------------------------------------------------------------
    # 1. Shear-sector evaluation with CORRECT environmental inputs
    #    rho_amb = ambient coarse-grained density of the region
    # ------------------------------------------------------------------
    # Ambient densities (g/cm^3): local interplanetary medium ~ solar wind,
    # ISM, and the galactic-habitat mean used for the cosmological channel.
    RHO_AMB_INTERPLANETARY = 1.0e-23   # ~5-10 protons/cm^3 at 1 AU
    RHO_AMB_ISM = 1.7e-24              # ~1 proton/cm^3
    RHO_AMB_GALACTIC = RHO_HALF_G_CM3  # 0.5 M_sun/pc^3 habitat mean
    RHO_AMB_VOID = 3.4e-26             # deep void
    RHO_AMB_COSMIC_MEAN = 4.0e-30      # critical-density mean today

    environments = {
        "earth_surface": dict(rho_amb=RHO_AMB_INTERPLANETARY, g=9.81,
                              note="body surface; density channel ~unsuppressed, gradient channel dominates"),
        "earth_orbit_1AU": dict(rho_amb=RHO_AMB_INTERPLANETARY, g=5.93e-3,
                                note="Solar acceleration at 1 AU"),
        "saturn_orbit_9p5AU": dict(rho_amb=RHO_AMB_INTERPLANETARY, g=6.5e-5,
                                 note="Cassini conjunction environment"),
        "solar_surface": dict(rho_amb=RHO_AMB_INTERPLANETARY, g=274.0,
                              note="deep heliospheric gradient"),
        "galactic_disk": dict(rho_amb=RHO_AMB_GALACTIC, g=1.2e-10,
                              note="halo/disk habitat: BOTH channels transitional"),
        "wide_binary": dict(rho_amb=RHO_AMB_GALACTIC, g=1.2e-10,
                            note="deep-MOND acceleration regime"),
        "cosmic_void": dict(rho_amb=RHO_AMB_VOID, g=1.0e-11,
                            note="ambient density and shear both below thresholds"),
        "cosmic_mean": dict(rho_amb=RHO_AMB_COSMIC_MEAN, g=1.0e-11,
                            note="homogeneous background: S_Sigma -> 1"),
    }

    env_results = {}
    for name, env in environments.items():
        r = s_sigma(env["rho_amb"], env["g"])
        r["note"] = env["note"]
        env_results[name] = r
        print_status(f"  {name:22s}: S_Sigma = {r['S_Sigma']:.3e}  "
                     f"(kinetic {r['kinetic_term']:.2e}, density {r['density_term']:.2e})", "INFO")

    cassini = env_results["saturn_orbit_9p5AU"]
    cassini_pass = cassini["S_Sigma"] < S_SIGMA_CASSINI_MAX
    print_status(f"Cassini gate: S_Sigma(Saturn) = {cassini['S_Sigma']:.2e} "
                 f"< {S_SIGMA_CASSINI_MAX:.1e} required -> {'PASS' if cassini_pass else 'FAIL'}", "SUCCESS")

    # ------------------------------------------------------------------
    # 2. Clock-amplitude-sector geometric factors S = (rho_bar/rho_T)^(1/3)
    # ------------------------------------------------------------------
    bodies = {
        "moon": 3.34,
        "sun": 1.41,
        "earth": 5.51,
        "white_dwarf": 1.0e6,
        "ucd_compact": 1.0e4,
        "neutron_star": 4.0e14,
    }
    geom_results = {}
    for name, rho_bar in bodies.items():
        s = s_geom(rho_bar)
        geom_results[name] = {
            "rho_bar_g_cm3": rho_bar,
            "S_geom": s,
            "regime": ("transitional (S_A clock channel active)" if s < 1.0
                       else "saturated interior (S_A suppressed)"),
        }
        print_status(f"  {name:14s}: rho_bar = {rho_bar:.2e} g/cm^3 -> S_geom = {s:.3f}", "INFO")

    # ------------------------------------------------------------------
    # 3. The reviewer's category error, computed and flagged
    #    Evaluating the AMBIENT channel at a MATERIAL density is out of domain.
    # ------------------------------------------------------------------
    rho_material_earth = 5.51  # g/cm^3 — Earth's mean material density
    s_out_of_domain = 1.0 / (1.0 + (rho_material_earth / RHO_HALF_G_CM3) ** 2)
    category_error = {
        "evaluation": "ambient-density channel at Earth's MATERIAL density 5.51 g/cm^3",
        "result_S": s_out_of_domain,
        "verdict": (
            "out-of-domain substitution: the density argument of S_Sigma is the "
            "coarse-grained environmental mean of the region, not the material "
            "density of an embedded body (Paper 0 Section 7: projections are not "
            "interchangeable). The correct ambient inputs give S_Sigma ~ 1e-21 "
            "via the gradient channel; the body's own material density enters the "
            "S_A sector through S_geom = 0.65 instead."
        ),
    }
    print_status(f"Out-of-domain check: ambient channel at 5.51 g/cm^3 -> S = {s_out_of_domain:.1e} "
                 f"(flagged, not the operator's prediction for Earth)", "WARNING")

    # ------------------------------------------------------------------
    # 4. R_T(M_earth) vs measured GNSS covariance scale
    # ------------------------------------------------------------------
    M_EARTH_KG = 5.9722e24
    r_t_earth = r_t_km(M_EARTH_KG)
    in_band = LAMBDA_T_MEASURED_KM[0] <= r_t_earth <= LAMBDA_T_MEASURED_KM[1]
    print_status(f"R_T(M_earth) = {r_t_earth:.0f} km vs measured lambda_T = "
                 f"{LAMBDA_T_MEASURED_KM[0]:.0f}-{LAMBDA_T_MEASURED_KM[1]:.0f} km "
                 f"-> {'within' if in_band else 'outside'} band", "SUCCESS" if in_band else "WARNING")

    # ------------------------------------------------------------------
    # 5. Quartic-completion bridge: m_eff ∝ rho^(1/3) -> lambda_C ∝ rho^(-1/3)
    #    One density-dependent mass spanning material to ambient densities.
    # ------------------------------------------------------------------
    RHO_ISM_REF = 1.7e-24          # g/cm^3
    LAMBDA_C_ISM_CM = 1.0 * PC_CM  # ~1 pc Compton length in the ISM (Paper 0 text)
    bridge = {}
    for label, rho in [
        ("cosmic_mean", RHO_AMB_COSMIC_MEAN),
        ("void", RHO_AMB_VOID),
        ("ISM", RHO_ISM_REF),
        ("rho_half_galactic", RHO_HALF_G_CM3),
        ("solar_interior", 150.0),
        ("rho_T_material", RHO_T),
        ("white_dwarf", 1.0e6),
    ]:
        lam_cm = LAMBDA_C_ISM_CM * (rho / RHO_ISM_REF) ** (-1.0 / 3.0)
        bridge[label] = {
            "rho_g_cm3": rho,
            "lambda_C_cm": lam_cm,
            "lambda_C_km": lam_cm / 1.0e5,
            "lambda_C_pc": lam_cm / PC_CM,
        }
        print_status(f"  lambda_C({label:18s}, rho={rho:.2e}) = {lam_cm:.3e} cm", "INFO")

    # ------------------------------------------------------------------
    # 6. Two-sector statement for Earth — the requested explicit S-values
    # ------------------------------------------------------------------
    earth = {
        "S_Sigma_shear_sector": env_results["earth_surface"]["S_Sigma"],
        "S_A_geometric_sector": geom_results["earth"]["S_geom"],
        "reading": (
            "S_Sigma(Earth) ~ 1e-21: shear/source-charge sector screened, as "
            "Cassini, LLR and Eoetvoes require. S_A(Earth) ~ 0.65: clock-amplitude "
            "sector at the geometric transition — the channel in which the GNSS "
            "clock-covariance signal is detected (Papers 1/14/33). The two are "
            "distinct projections of one scalar configuration; the ~24-order "
            "separation of rho_T and rho_half is the physical hierarchy between "
            "solid-matter and cosmic ambient densities, not an inconsistency."
        ),
    }

    payload = {
        "step": STEP_ID,
        "description": (
            "Two-sector screening evaluation: rho_T (material saturation, S_A "
            "geometric sector) and rho_half (ambient half-suppression, S_Sigma "
            "density channel) parameterize different response projections of one "
            "scalar configuration"
        ),
        "status": "completed",
        "scales": {
            "rho_T_g_cm3": RHO_T,
            "rho_half_Msun_pc3": RHO_HALF_MSUN_PC3,
            "rho_half_g_cm3": RHO_HALF_G_CM3,
            "a_t_m_s2": A_T,
            "n": N_SCREEN,
            "hierarchy_rho_T_over_rho_half": RHO_T / RHO_HALF_G_CM3,
        },
        "environments_S_Sigma": env_results,
        "geometric_S_A_factors": geom_results,
        "category_error_check": category_error,
        "earth_two_sector": earth,
        "cassini": {
            "S_Sigma_saturn": cassini["S_Sigma"],
            "required_max": S_SIGMA_CASSINI_MAX,
            "passes": cassini_pass,
        },
        "R_T_earth_km": r_t_earth,
        "lambda_T_measured_km": list(LAMBDA_T_MEASURED_KM),
        "R_T_within_measured_band": in_band,
        "compton_bridge_quartic": bridge,
        "interpretation": (
            "The GNSS clock-covariance signal lives in the S_A (clock-amplitude) "
            "sector at the geometric transition S ~ 0.65, while the shear sector "
            "is screened to S_Sigma ~ 1e-21 at Earth's surface by the gradient "
            "channel of the same covariant operator. No contradiction: different "
            "response functionals, one field configuration."
        ),
        "timestamp": int(time.time()),
    }

    write_json(step_json_path(STEP_ID), payload)
    print_status(f"Step {STEP_ID} completed", "SUCCESS")
    return payload


if __name__ == "__main__":
    print(run())
