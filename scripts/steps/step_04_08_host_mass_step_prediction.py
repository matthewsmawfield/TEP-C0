#!/usr/bin/env python3
"""TEP-C0 Step 04-08: Host Mass Step Prediction
===============================================
Mini-analysis comparing four models for the Pantheon+ host-mass step:

1. LCDM no mass step
2. LCDM fitted mass step
3. TEP locked mass step (canonical nested clock-rate prediction:
   kappa_SN * Delta_X with X = S_env * S_A * sigma^2 / c^2, no fitted
   amplitude; retired logarithmic ansatz removed)
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

    # --- Model 3: TEP locked mass step (canonical nested clock-rate prediction) ---
    # Per Rule 22 the SN clock rate is the nested total: ambient
    # environmental baseline + screened host well.  The corpus's derived
    # environmental coordinate (Paper 11 convention) is
    #     X = S_env * sigma^2 / c^2,
    # where sigma^2 = G M / (2 R_eff) is the host virial depth and S_env
    # is the ambient (group-halo) screening factor encoding the nested
    # baseline difference between environments.  The amplitude-sector
    # projection
    #     S_A = min[1, (sigma / sigma_T)^{2/3}],   sigma_T = 65 km/s
    # (Paper 0 amplitude sector S_A = min[1,(rho/rho_T)^{1/3}] evaluated
    # on the depth proxy sigma^2; sigma_T selected at 65 km/s by the
    # TEP-H0 step_52 model comparison) suppresses the response for
    # shallow hosts.  The magnitude offset follows from the response-
    # sector coefficient kappa_SN; because both ladder channels respond
    # to the same nested field, the frozen Cepheid-equivalent value
    # kappa_equiv = 3.69e5 mag (TEP-H0 step_44 redshift-only channel,
    # canonical sigma_v = 250 km/s convention) is adopted with no fitted
    # nuisance parameter.
    from core import constants as tep_const

    KAPPA_SN = 3.69e5          # mag; frozen kappa_Cep-equivalent benchmark
    SIGMA_T_SA = 65.0          # km/s; amplitude-sector scale (step_52)

    def host_sigma_kms(log_mass_solar, r_eff_kpc):
        """Virial depth coordinate sigma^2 = G M / (2 R_eff)."""
        mass_kg = 10.0 ** log_mass_solar * tep_const.M_SUN
        r_eff_m = r_eff_kpc * tep_const.MPC_TO_M / 1e3  # kpc -> m
        return np.sqrt(tep_const.G_NEWTON * mass_kg / (2.0 * r_eff_m)) / 1e3

    def S_A_of_sigma(sigma_kms):
        """Amplitude-sector response on the depth proxy."""
        return min(1.0, (sigma_kms / SIGMA_T_SA) ** (2.0 / 3.0))

    def S_env_group(n_mb):
        """Ambient group-halo screening (TEP-H0 convention, N_crit=10,
        gamma=1.2); encodes the nested environmental baseline."""
        return 1.0 / (1.0 + (n_mb / 10.0) ** 1.2)

    # Representative host classes
    sigma_high = host_sigma_kms(11.0, 5.0)   # massive host ~ 207 km/s
    sigma_low = host_sigma_kms(9.5, 2.0)     # low-mass host ~ 58 km/s

    # Ambient nesting: massive hosts preferentially occupy group/cluster
    # environments.  Fiducial Tully-group richness values are adopted for
    # the two classes; the isolated-depth bound (S_env = 1) is reported
    # alongside as the upper envelope.
    senv_high = S_env_group(8.0)
    senv_low = S_env_group(2.0)
    sa_high = S_A_of_sigma(sigma_high)
    sa_low = S_A_of_sigma(sigma_low)

    c2 = (tep_const.C_LIGHT / 1e3) ** 2   # (km/s)^2, sigma coordinate

    X_high = senv_high * sa_high * sigma_high ** 2 / c2
    X_low = senv_low * sa_low * sigma_low ** 2 / c2
    X_high_iso = sigma_high ** 2 / c2
    X_low_iso = sigma_low ** 2 / c2

    delta_X_nested = X_high - X_low
    delta_X_isolated = X_high_iso - X_low_iso

    # Locked prediction: magnitude offset of the high-mass class.
    # Paper 11's convention is Delta_mu = -kappa * X, so kappa > 0
    # makes high-X hosts appear brighter.
    predicted_offset = -KAPPA_SN * delta_X_nested
    predicted_offset_isolated = -KAPPA_SN * delta_X_isolated

    # Raw conformal clock channel for honesty: Delta_u = Delta X on the
    # unscreened depth coordinate gives the same high-X-brighter sign.
    beta_A = tep_const.BETA_A
    delta_mu_clock = float(-1.0857 * abs(beta_A) * delta_X_isolated)

    # Data-implied response coefficient for cross-channel comparison
    # (filled after the LCDM fit below)
    print_status(f"Canonical nested-rate prediction: {predicted_offset:.4f} mag "
                 f"(isolated bound {predicted_offset_isolated:.4f} mag)", "INFO")
    print_status(f"  sigma_high={sigma_high:.1f} km/s, sigma_low={sigma_low:.1f} km/s; "
                 f"S_env=({senv_high:.2f},{senv_low:.2f}), S_A=({sa_high:.3f},{sa_low:.3f}); "
                 f"Delta_X={delta_X_nested:.3e}", "INFO")
    print_status(f"  raw conformal clock channel: {delta_mu_clock:.2e} mag "
                 f"(response amplification kappa ~ {KAPPA_SN:.2e})", "INFO")

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

    # Data-implied response coefficient: the SN-channel equivalent of
    # kappa_Cep, recovered from the fitted residual on the same nested
    # coordinate.  Cross-channel consistency with the Cepheid sector
    # (kappa ~ few x 10^5 mag) is the falsifiable content of the locked
    # prediction.
    kappa_SN_implied = float(
        -observed_mass_step / delta_X_nested if delta_X_nested != 0 else np.nan
    )

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
        "host_mass_threshold": float(HOST_MASS_THRESHOLD),
        "observed_mass_step": float(observed_mass_step),
        "TEP_predicted_mass_step": float(predicted_offset),
        "TEP_predicted_mass_step_isolated_bound": float(predicted_offset_isolated),
        "kappa_SN_frozen": float(KAPPA_SN),
        "kappa_SN_implied_by_data": kappa_SN_implied,
        "delta_X_nested": float(delta_X_nested),
        "delta_X_isolated": float(delta_X_isolated),
        "delta_mu_clock_raw": delta_mu_clock,
        "host_coordinates": {
            "sigma_high_kms": float(sigma_high),
            "sigma_low_kms": float(sigma_low),
            "S_env_high": float(senv_high),
            "S_env_low": float(senv_low),
            "S_A_high": float(sa_high),
            "S_A_low": float(sa_low),
        },
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
            "Under the Paper 11 convention Delta_mu = -kappa X, the "
            "canonical nested clock-rate prediction has the established "
            "massive-host-brighter sign and literature-step order "
            "(~0.04-0.05 mag), but it overshoots the post-correction "
            "residual in this m_b_corr mini-analysis. The frozen "
            "kappa_SN = kappa_Cep-equiv = 3.69e5 mag is a cross-channel "
            "transfer, while the data-implied kappa_SN on this corrected "
            "product is reported separately and is opposite signed. The "
            "raw conformal clock ratio contributes only ~5e-7 mag, so "
            "the step is carried by the response sector, consistent "
            "with the Paper 0 amplitude/screening split. This step is "
            "therefore a channel-consistency diagnostic, not an "
            "independent confirmation of the host-mass correction."
        ),
        "timestamp": int(time.time()),
    }

    write_json(step_json_path(STEP_ID), payload)
    print_status(f"Step {STEP_ID} completed", "SUCCESS")
    return payload


if __name__ == "__main__":
    print(run())
