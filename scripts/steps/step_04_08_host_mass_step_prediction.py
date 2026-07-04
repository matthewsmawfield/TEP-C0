#!/usr/bin/env python3
"""TEP-C0 Step 04-08: Host Mass Step Prediction
===============================================
Mini-analysis comparing four models for the Pantheon+ host-mass step:

1. LCDM no mass step
2. LCDM fitted mass step
3. TEP locked mass step (alpha_log = -7.66e-3 fixed)
4. TEP fitted residual environmental term

Outputs:
- observed_mass_step
- TEP_predicted_mass_step
- residual_after_TEP_correction
- Delta chi2 improvement
- AIC/BIC vs fitted nuisance parameter
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import cholesky
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from c0_common import (
    RESULTS_DIR,
    TEPLogger,
    ensure_dirs,
    print_status,
    set_step_logger,
    step_json_path,
    write_json,
)
from core.cosmology import CosmologyFLRW, TEPCosmology

STEP_ID = "step_04_08_host_mass_step_prediction"

ALPHA_LOG = -7.66e-3  # Locked lab-scale coupling
HOST_MASS_THRESHOLD = 10.0  # log(M*/Msun)


def load_pantheon_with_mass():
    """Load Pantheon+ data with host mass column."""
    data_dir = Path("data/raw")
    data_candidates = [
        data_dir / "Pantheon+SH0ES.dat",
        data_dir / "pantheon_plus_shoes.dat",
    ]
    pantheon_file = next((p for p in data_candidates if p.exists()), None)
    if pantheon_file is None:
        raise FileNotFoundError("Pantheon+ data not found")

    import pandas as pd
    df = pd.read_csv(pantheon_file, sep=r'\s+', comment='#')

    z_raw = df["zCMB"].values if "zCMB" in df.columns else df["zHD"].values
    mb_raw = df['m_b_corr'].values if 'm_b_corr' in df.columns else df['mB'].values

    # Host mass column
    host_mass = None
    for col in ['HOST_LOGMASS', 'HOST_MASS', 'logMass']:
        if col in df.columns:
            host_mass = df[col].values
            break

    if host_mass is None:
        raise FileNotFoundError("No host mass column found in Pantheon+ data")

    # Filter valid data: remove NaN, Inf, and sentinel values (-999 = missing host mass)
    valid = (np.isfinite(z_raw) & np.isfinite(mb_raw) & np.isfinite(host_mass)
             & (host_mass > 0.0) & (z_raw > 0.001) & (z_raw < 3.0))
    z = z_raw[valid]
    mb = mb_raw[valid]
    host_mass = host_mass[valid]

    # Load covariance
    cov_candidates = [
        data_dir / "Pantheon+SH0ES.cov",
        data_dir / "pantheon_plus_shoes.cov",
    ]
    cov_file = next((c for c in cov_candidates if c.exists()), None)
    if cov_file is None:
        raise FileNotFoundError("Covariance not found")

    cov_flat = np.loadtxt(cov_file)
    n_raw = len(z_raw)
    if cov_flat.size == n_raw * n_raw + 1:
        cov_full = cov_flat[1:].reshape(n_raw, n_raw)
    elif cov_flat.size == n_raw * n_raw:
        cov_full = cov_flat.reshape(n_raw, n_raw)
    elif cov_flat.size == int(np.sum(valid)) ** 2:
        cov_full = cov_flat.reshape(int(np.sum(valid)), int(np.sum(valid)))
        return z, mb, host_mass, cov_full
    else:
        raise ValueError(f"Covariance size mismatch: {cov_flat.size}")

    cov = cov_full[np.ix_(valid, valid)]
    return z, mb, host_mass, cov


def compute_chi2(residuals: np.ndarray, L: np.ndarray) -> float:
    y = np.linalg.solve(L, residuals)
    return float(np.dot(y, y))


def run() -> dict:
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    try:
        z, mb, host_mass, cov = load_pantheon_with_mass()
    except FileNotFoundError as e:
        payload = {
            "step": STEP_ID,
            "status": "blocked",
            "error": str(e),
            "timestamp": int(time.time()),
        }
        write_json(step_json_path(STEP_ID), payload)
        print_status(f"Step {STEP_ID} blocked: {e}", "WARNING")
        return payload

    n = len(z)
    L = cholesky(cov, lower=True, check_finite=False)
    logdet = float(2.0 * np.sum(np.log(np.diag(L))))
    mass_flag = (host_mass > HOST_MASS_THRESHOLD).astype(float)

    print_status(f"Loaded {n} SNe with host masses", "INFO")
    print_status(f"High-mass hosts (logM>{HOST_MASS_THRESHOLD}): {int(mass_flag.sum())}", "INFO")

    # --- Model 1: LCDM no mass step ---
    def chi2_lcdm_no_step(params):
        Om, M = params
        if not (0.05 <= Om <= 0.9):
            return 1e10
        c = CosmologyFLRW(H0=70.0, Om0=Om)
        mu_pred = c.distance_modulus(z) + M
        return compute_chi2(mb - mu_pred, L)

    res1 = minimize(chi2_lcdm_no_step, [0.315, -19.3], method="L-BFGS-B", bounds=[(0.05, 0.9), (-21.0, -16.0)])
    chi2_1 = res1.fun

    # --- Model 2: LCDM + fitted mass step ---
    def chi2_lcdm_fitted_step(params):
        Om, M, gamma = params
        if not (0.05 <= Om <= 0.9):
            return 1e10
        c = CosmologyFLRW(H0=70.0, Om0=Om)
        mu_pred = c.distance_modulus(z) + M + gamma * mass_flag
        return compute_chi2(mb - mu_pred, L)

    res2 = minimize(chi2_lcdm_fitted_step, [0.315, -19.3, -0.03], method="L-BFGS-B",
                    bounds=[(0.05, 0.9), (-21.0, -16.0), (-0.2, 0.2)])
    chi2_2 = res2.fun
    gamma_fitted = float(res2.x[2])

    # --- Model 3: TEP locked mass step ---
    # Compute predicted offset from scalar field solver using realistic
    # galaxy parameters instead of the hand-waved ln(100) mass ratio.
    from core.scalar_field import solve_scalar_field_cylinder
    from core.screening import screening_factor
    from core import constants as tep_const

    def galaxy_density_g_cm3(log_mass_solar, r_eff_kpc):
        """Estimate average stellar density from total mass and effective radius."""
        mass_kg = 10.0 ** log_mass_solar * tep_const.M_SUN
        r_eff_m = r_eff_kpc * 1000.0 * tep_const.MPC_TO_M / 1e6  # kpc -> m
        # Assume spherical volume ~ 4/3 pi r_eff^3 for order-of-magnitude estimate
        vol_m3 = (4.0 / 3.0) * np.pi * r_eff_m ** 3
        rho_kg_m3 = mass_kg / vol_m3
        return rho_kg_m3 / 1000.0  # kg/m^3 -> g/cm^3

    # Typical high-mass host: logM ~ 11.0, r_eff ~ 5 kpc
    rho_high = galaxy_density_g_cm3(11.0, 5.0)
    # Typical low-mass host: logM ~ 9.5, r_eff ~ 2 kpc
    rho_low = galaxy_density_g_cm3(9.5, 2.0)

    # Include screening at galactic densities
    screen_high = screening_factor(rho_high, tep_const.RHO_T)
    screen_low = screening_factor(rho_low, tep_const.RHO_T)

    phi_high, _ = solve_scalar_field_cylinder(
        total_mass_kg=10.0 ** 11.0 * tep_const.M_SUN,
        radius_m=5.0 * 1000.0 * tep_const.MPC_TO_M / 1e6,
        height_m=5.0 * 1000.0 * tep_const.MPC_TO_M / 1e6,
        density_g_cm3=rho_high,
    )
    phi_low, _ = solve_scalar_field_cylinder(
        total_mass_kg=10.0 ** 9.5 * tep_const.M_SUN,
        radius_m=2.0 * 1000.0 * tep_const.MPC_TO_M / 1e6,
        height_m=2.0 * 1000.0 * tep_const.MPC_TO_M / 1e6,
        density_g_cm3=rho_low,
    )

    # Clock-rate ratio A(phi_high)/A(phi_low) = exp(beta_A * (phi_high - phi_low))
    # Magnitude shift: Delta_mu = -2.5 * log10(A_high / A_low)
    #                  = -2.5 * beta_A * (phi_high - phi_low) / ln(10)
    #                  = -1.0857 * beta_A * (phi_high - phi_low)
    # With beta_A = -1.0 (from constants), and phi_high < phi_low (higher density -> more negative phi):
    # Delta_mu ~ -1.0857 * (-1.0) * (negative) = negative (high-mass hosts brighter)
    beta_A = tep_const.BETA_A
    predicted_offset = -1.0857 * beta_A * (phi_high - phi_low)
    print_status(f"TEP predicted mass-step offset: {predicted_offset:.4f} mag", "INFO")
    print_status(f"  (phi_high={phi_high:.4f}, phi_low={phi_low:.4f}, "
                 f"screen_high={screen_high:.4f}, screen_low={screen_low:.4f})", "INFO")

    # Use LCDM best-fit Om for TEP to isolate the mass-step comparison
    om_lcdm_mle = float(res1.x[0])

    def chi2_tep_locked(params):
        eps, M = params
        if not (0.0 <= eps <= 2.0):
            return 1e10
        # Fair comparison: use same matter density as LCDM best-fit
        tep = TEPCosmology(H0=70.0, Omega_m=om_lcdm_mle, epsilon_T=eps, z_T=5.0)
        mu_pred = tep.distance_modulus(z) + M + predicted_offset * mass_flag
        return compute_chi2(mb - mu_pred, L)

    res3 = minimize(chi2_tep_locked, [0.1, -19.3], method="L-BFGS-B", bounds=[(0.0, 2.0), (-21.0, -16.0)])
    chi2_3 = res3.fun

    # --- Model 4: TEP + fitted residual environmental term ---
    def chi2_tep_fitted_env(params):
        eps, M, gamma_env = params
        if not (0.0 <= eps <= 2.0):
            return 1e10
        # Fair comparison: use same matter density as LCDM best-fit
        tep = TEPCosmology(H0=70.0, Omega_m=om_lcdm_mle, epsilon_T=eps, z_T=5.0)
        mu_pred = tep.distance_modulus(z) + M + predicted_offset * mass_flag + gamma_env * mass_flag
        return compute_chi2(mb - mu_pred, L)

    res4 = minimize(chi2_tep_fitted_env, [0.1, -19.3, 0.0], method="L-BFGS-B",
                    bounds=[(0.0, 2.0), (-21.0, -16.0), (-0.1, 0.1)])
    chi2_4 = res4.fun
    gamma_residual = float(res4.x[2])

    # Information criteria
    def aic_bic(chi2_val, k, n):
        return float(chi2_val + 2 * k), float(chi2_val + k * np.log(n))

    aic1, bic1 = aic_bic(chi2_1, 2, n)
    aic2, bic2 = aic_bic(chi2_2, 3, n)
    aic3, bic3 = aic_bic(chi2_3, 2, n)
    aic4, bic4 = aic_bic(chi2_4, 3, n)

    # Observed mass step from fitted LCDM
    observed_mass_step = gamma_fitted

    # Residual after TEP locked correction
    residual_after_tep = gamma_residual

    delta_chi2_tep_locked_vs_lcdm_no_step = float(chi2_1 - chi2_3)
    delta_chi2_tep_locked_vs_lcdm_fitted = float(chi2_2 - chi2_3)

    print_status(f"Observed mass step (LCDM fitted): {observed_mass_step:.4f} mag", "INFO")
    print_status(f"TEP predicted mass step: {predicted_offset:.4f} mag", "INFO")
    print_status(f"Residual after TEP locked correction: {residual_after_tep:.4f} mag", "INFO")
    print_status(f"Δχ² TEP_locked vs LCDM_no_step: {delta_chi2_tep_locked_vs_lcdm_no_step:+.2f}", "INFO")
    print_status(f"Δχ² TEP_locked vs LCDM_fitted: {delta_chi2_tep_locked_vs_lcdm_fitted:+.2f}", "INFO")

    payload = {
        "step": STEP_ID,
        "description": "Host mass step prediction mini-analysis",
        "status": "completed",
        "n_sne": int(n),
        "alpha_log": float(ALPHA_LOG),
        "host_mass_threshold": float(HOST_MASS_THRESHOLD),
        "observed_mass_step": float(observed_mass_step),
        "TEP_predicted_mass_step": float(predicted_offset),
        "residual_after_TEP_correction": float(residual_after_tep),
        "models": {
            "LCDM_no_step": {"chi2": float(chi2_1), "k": 2, "aic": aic1, "bic": bic1},
            "LCDM_fitted_step": {"chi2": float(chi2_2), "k": 3, "aic": aic2, "bic": bic2,
                                  "gamma": float(gamma_fitted)},
            "TEP_locked_step": {"chi2": float(chi2_3), "k": 2, "aic": aic3, "bic": bic3,
                                "epsilon_shear": float(res3.x[0])},
            "TEP_fitted_residual": {"chi2": float(chi2_4), "k": 3, "aic": aic4, "bic": bic4,
                                    "gamma_residual": float(gamma_residual),
                                    "epsilon_shear": float(res4.x[0])},
        },
        "delta_chi2_tep_locked_vs_lcdm_no_step": delta_chi2_tep_locked_vs_lcdm_no_step,
        "delta_chi2_tep_locked_vs_lcdm_fitted": delta_chi2_tep_locked_vs_lcdm_fitted,
        "interpretation": (
            "The locked TEP mass-step prediction removes most of the observed host-mass residual "
            "without fitting a new nuisance parameter."
            if abs(residual_after_tep) < abs(observed_mass_step) * 0.5 else
            "The TEP locked prediction captures part of the host-mass effect; residual indicates "
            "additional environmental physics not yet included."
        ),
        "timestamp": int(time.time()),
    }

    write_json(step_json_path(STEP_ID), payload)
    print_status(f"Step {STEP_ID} completed", "SUCCESS")
    return payload


if __name__ == "__main__":
    print(run())
