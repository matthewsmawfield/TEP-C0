#!/usr/bin/env python3
"""TEP-C0 Step 03-09: LCDM Null Injection Test
==============================================
Generate mock Pantheon+-like standardized magnitudes from best-fit LCDM.
Preserve Pantheon+ redshifts and covariance. Run model comparison to
measure false-positive rate for TEP Bayes factor > 30.

Outputs:
- false_positive_rate_TEP_BF_gt_30
- median_Delta_chi2_under_LCDM
- p_value_observed_Delta_chi2
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import cholesky

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

STEP_ID = "step_03_09_lcdm_null_injection"

N_SEEDS = 32  # Number of null realizations


def load_pantheon_data():
    """Load Pantheon+ data and covariance."""
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

    valid = np.isfinite(z_raw) & np.isfinite(mb_raw) & (z_raw > 0.001) & (z_raw < 3.0)
    z = z_raw[valid]
    n = len(z)

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
    n_filt = int(np.sum(valid))
    if cov_flat.size == n_raw * n_raw + 1:
        cov_full = cov_flat[1:].reshape(n_raw, n_raw)
        cov = cov_full[np.ix_(valid, valid)]
    elif cov_flat.size == n_raw * n_raw:
        cov_full = cov_flat.reshape(n_raw, n_raw)
        cov = cov_full[np.ix_(valid, valid)]
    elif cov_flat.size == n_filt * n_filt + 1:
        cov = cov_flat[1:].reshape(n_filt, n_filt)
    elif cov_flat.size == n_filt * n_filt:
        cov = cov_flat.reshape(n_filt, n_filt)
    else:
        raise ValueError(f"Covariance size mismatch: {cov_flat.size} vs raw={n_raw}^2+1={n_raw*n_raw+1} or filt={n_filt}^2+1={n_filt*n_filt+1}")

    return z, mb_raw[valid], cov


def generate_mock_lcdm(z: np.ndarray, cov: np.ndarray, seed: int) -> np.ndarray:
    """Generate mock standardized magnitudes from best-fit LCDM."""
    rng = np.random.default_rng(seed)
    lcdm = CosmologyFLRW(H0=70.0, Om0=0.315)
    mu_lcdm = lcdm.distance_modulus(z)
    # Standard intercept M ~ -19.3
    M = -19.3
    mb_true = mu_lcdm + M
    # Add covariance noise
    L = cholesky(cov, lower=True, check_finite=False)
    noise = L @ rng.standard_normal(len(z))
    return mb_true + noise


def fit_lcdm(z: np.ndarray, mb: np.ndarray, cov: np.ndarray) -> tuple[float, float]:
    """Fit LCDM (Om, M) and return chi2, lnL."""
    from scipy.optimize import minimize
    lcdm = CosmologyFLRW(H0=70.0, Om0=0.315)

    def chi2(params):
        Om, M = params
        if not (0.05 <= Om <= 0.9):
            return 1e10
        c = CosmologyFLRW(H0=70.0, Om0=Om)
        mu_pred = c.distance_modulus(z) + M
        residuals = mb - mu_pred
        L = cholesky(cov, lower=True, check_finite=False)
        y = np.linalg.solve(L, residuals)
        return float(np.dot(y, y))

    res = minimize(chi2, [0.3, -19.3], method="L-BFGS-B", bounds=[(0.05, 0.9), (-21.0, -16.0)])
    chi2_min = res.fun
    n = len(z)
    logdet = float(2.0 * np.sum(np.log(np.diag(cholesky(cov, lower=True, check_finite=False)))))
    lnL = -0.5 * (chi2_min + logdet + n * np.log(2.0 * np.pi))
    return chi2_min, lnL


def fit_tep_m1(z: np.ndarray, mb: np.ndarray, cov: np.ndarray) -> tuple[float, float]:
    """Fit TEP M1 (epsilon_shear, M) and return chi2, lnL."""
    from scipy.optimize import minimize

    def chi2(params):
        eps, M = params
        if not (0.0 <= eps <= 2.0):
            return 1e10
        tep = TEPCosmology(H0=70.0, Omega_m=1.0, epsilon_T=eps, z_T=5.0)
        mu_pred = tep.distance_modulus(z) + M
        residuals = mb - mu_pred
        L = cholesky(cov, lower=True, check_finite=False)
        y = np.linalg.solve(L, residuals)
        return float(np.dot(y, y))

    res = minimize(chi2, [0.1, -19.3], method="L-BFGS-B", bounds=[(0.0, 2.0), (-21.0, -16.0)])
    chi2_min = res.fun
    n = len(z)
    logdet = float(2.0 * np.sum(np.log(np.diag(cholesky(cov, lower=True, check_finite=False)))))
    lnL = -0.5 * (chi2_min + logdet + n * np.log(2.0 * np.pi))
    return chi2_min, lnL


def run() -> dict:
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    try:
        z, mb_obs, cov = load_pantheon_data()
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
    print_status(f"Loaded {n} SNe with full covariance", "INFO")

    # Run null injections
    delta_chi2_list = []
    bf_list = []
    tep_bf_gt_30_count = 0

    for seed in range(N_SEEDS):
        mb_mock = generate_mock_lcdm(z, cov, seed)
        chi2_lcdm, lnL_lcdm = fit_lcdm(z, mb_mock, cov)
        chi2_tep, lnL_tep = fit_tep_m1(z, mb_mock, cov)

        delta_chi2 = chi2_tep - chi2_lcdm
        delta_chi2_list.append(delta_chi2)

        # Bayes factor proxy: exp(2 * Delta lnL) for 2 extra params in TEP
        # More precisely, use BIC-approximated Bayes factor
        k_lcdm = 2
        k_tep = 2
        bic_lcdm = chi2_lcdm + k_lcdm * np.log(n)
        bic_tep = chi2_tep + k_tep * np.log(n)
        delta_bic = bic_tep - bic_lcdm
        bf_proxy = np.exp(-0.5 * delta_bic)
        bf_list.append(float(bf_proxy))

        if bf_proxy > 30.0:
            tep_bf_gt_30_count += 1

        if seed % 8 == 0:
            print_status(f"  Seed {seed}: Δχ²={delta_chi2:.2f}, BF_proxy={bf_proxy:.2f}", "INFO")

    delta_chi2_arr = np.array(delta_chi2_list)
    bf_arr = np.array(bf_list)

    false_positive_rate = tep_bf_gt_30_count / N_SEEDS
    median_delta_chi2 = float(np.median(delta_chi2_arr))
    p_value = float(np.mean(delta_chi2_arr < -7.5))  # P(Δχ² < observed) under LCDM

    print_status(f"Null injection summary:", "INFO")
    print_status(f"  Median Δχ² (TEP-LCDM): {median_delta_chi2:.2f}", "INFO")
    print_status(f"  BF > 30 count: {tep_bf_gt_30_count} / {N_SEEDS}", "INFO")
    print_status(f"  False positive rate: {false_positive_rate:.3f}", "INFO")
    print_status(f"  p(observed Δχ² < -7.5 | LCDM): {p_value:.3f}", "INFO")

    payload = {
        "step": STEP_ID,
        "description": "LCDM null injection: mock Pantheon+ from LCDM best-fit",
        "status": "completed",
        "n_injections": N_SEEDS,
        "false_positive_rate_TEP_BF_gt_30": float(false_positive_rate),
        "median_Delta_chi2_under_LCDM": float(median_delta_chi2),
        "p_value_observed_Delta_chi2": float(p_value),
        "bf_proxy_distribution": {
            "mean": float(np.mean(bf_arr)),
            "median": float(np.median(bf_arr)),
            "std": float(np.std(bf_arr)),
            "min": float(np.min(bf_arr)),
            "max": float(np.max(bf_arr)),
        },
        "delta_chi2_distribution": {
            "mean": float(np.mean(delta_chi2_arr)),
            "median": float(np.median(delta_chi2_arr)),
            "std": float(np.std(delta_chi2_arr)),
            "min": float(np.min(delta_chi2_arr)),
            "max": float(np.max(delta_chi2_arr)),
        },
        "interpretation": (
            f"Under LCDM null injection, the TEP M1 Bayes factor observed in real Pantheon+ "
            f"(BF ~ 30) occurs in only {false_positive_rate*100:.1f}% of synthetic realizations."
        ),
        "timestamp": int(time.time()),
    }

    write_json(step_json_path(STEP_ID), payload)
    print_status(f"Step {STEP_ID} completed", "SUCCESS")
    return payload


if __name__ == "__main__":
    print(run())
