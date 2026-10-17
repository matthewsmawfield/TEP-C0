#!/usr/bin/env python3
"""TEP-C0 Step 03-10: Pantheon+ Subset Robustness
=================================================
Leave-one-survey-out and redshift-window robustness tests.

Runs:
- all Pantheon+
- low-z only
- z > 0.01, z > 0.023, z > 0.05
- remove SH0ES anchors
- remove each major survey one at a time
- high-z only

For each, report Delta chi2 M1 vs LCDM, best epsilon_shear_los, best z_T.
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

STEP_ID = "step_03_10_pantheon_subset_robustness"


def load_pantheon_data():
    """Load Pantheon+ data, covariance, and survey identifiers."""
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

    # zHD is the corpus-standard Pantheon+ Hubble-diagram redshift
    # (peculiar-velocity corrected); the covariance's VPEC terms are
    # defined for that convention.
    z_raw = df["zHD"].values if "zHD" in df.columns else df["zCMB"].values
    mb_raw = df['m_b_corr'].values if 'm_b_corr' in df.columns else df['mB'].values

    # Survey identifier: use 'SURVEY' or 'survey' column if present
    survey_col = None
    for col in ['SURVEY', 'survey', 'IDSURVEY']:
        if col in df.columns:
            survey_col = df[col].values
            break

    calibrator_raw = (
        df['IS_CALIBRATOR'].values if 'IS_CALIBRATOR' in df.columns else None
    )

    valid = np.isfinite(z_raw) & np.isfinite(mb_raw) & (z_raw > 0.001) & (z_raw < 3.0)
    z = z_raw[valid]
    zhel = df["zHEL"].values[valid]
    if not np.all(np.isfinite(zhel)) or np.any(zhel <= -1):
        raise ValueError("Invalid observer redshift")
    # Fit brightness in the distance-redshift frame; released data unchanged.
    mb = mb_raw[valid] - 5.0 * np.log10((1.0 + zhel) / (1.0 + z))
    survey = survey_col[valid] if survey_col is not None else np.full(len(z), "UNKNOWN")
    is_calibrator = calibrator_raw[valid].astype(bool) if calibrator_raw is not None else None

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
        return z, mb, survey, is_calibrator, cov_full
    else:
        raise ValueError(f"Covariance size mismatch: {cov_flat.size}")

    cov = cov_full[np.ix_(valid, valid)]
    return z, mb, survey, is_calibrator, cov


def fit_models(z: np.ndarray, mb: np.ndarray, cov: np.ndarray) -> dict:
    """Fit LCDM and TEP M1, return chi2 and best-fit params."""
    n = len(z)
    L = cholesky(cov, lower=True, check_finite=False)
    logdet = float(2.0 * np.sum(np.log(np.diag(L))))

    def chi2_lcdm(params):
        Om, M = params
        if not (0.05 <= Om <= 0.9):
            return 1e10
        c = CosmologyFLRW(H0=70.0, Om0=Om)
        mu_pred = c.distance_modulus(z) + M
        residuals = mb - mu_pred
        y = np.linalg.solve(L, residuals)
        return float(np.dot(y, y))

    def chi2_tep(params):
        eps, M = params
        if not (0.0 <= eps <= 2.0):
            return 1e10
        tep = TEPCosmology(H0=70.0, Omega_m=1.0, epsilon_T=eps, z_T=5.0)
        mu_pred = tep.distance_modulus(z) + M
        residuals = mb - mu_pred
        y = np.linalg.solve(L, residuals)
        return float(np.dot(y, y))

    res_lcdm = minimize(chi2_lcdm, [0.315, -19.3], method="L-BFGS-B", bounds=[(0.05, 0.9), (-21.0, -16.0)])
    res_tep = minimize(chi2_tep, [0.1, -19.3], method="L-BFGS-B", bounds=[(0.0, 2.0), (-21.0, -16.0)])

    chi2_lcdm_val = res_lcdm.fun
    chi2_tep_val = res_tep.fun
    delta_chi2 = chi2_tep_val - chi2_lcdm_val
    lnL_lcdm = -0.5 * (chi2_lcdm_val + logdet + n * np.log(2.0 * np.pi))
    lnL_tep = -0.5 * (chi2_tep_val + logdet + n * np.log(2.0 * np.pi))

    return {
        "chi2_lcdm": float(chi2_lcdm_val),
        "chi2_tep": float(chi2_tep_val),
        "delta_chi2": float(delta_chi2),
        "delta_lnL": float(lnL_tep - lnL_lcdm),
        "bestfit_lcdm": {"Om": float(res_lcdm.x[0]), "M": float(res_lcdm.x[1])},
        "bestfit_tep": {"epsilon_shear": float(res_tep.x[0]), "M": float(res_tep.x[1])},
        "n_sne": int(n),
    }


def run() -> dict:
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    try:
        z_all, mb_all, survey_all, is_calibrator, cov_all = load_pantheon_data()
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

    # Define subsets
    subsets = {}

    # Redshift cuts
    subsets["all"] = np.ones(len(z_all), dtype=bool)
    subsets["low_z"] = z_all < 0.1
    subsets["z_gt_0.01"] = z_all > 0.01
    subsets["z_gt_0.023"] = z_all > 0.023
    subsets["z_gt_0.05"] = z_all > 0.05
    subsets["high_z"] = z_all > 0.5

    # Survey cuts (if survey info available)
    unique_surveys = np.unique(survey_all)
    if len(unique_surveys) > 1 and unique_surveys[0] != "UNKNOWN":
        for surv in unique_surveys:
            subsets[f"no_{surv}"] = survey_all != surv
        # SH0ES anchors: the z>0.023 proxy removes the 77 Cepheid
        # calibrators plus ~250 additional low-z Hubble-flow SNe.
        subsets["no_SH0ES"] = z_all > 0.023  # proxy: remove very local anchors
        if is_calibrator is not None:
            # Exact anchor removal: drops only the 77 IS_CALIBRATOR
            # Cepheid-anchored SNe, retaining the rest of the low-z flow.
            subsets["no_calibrators"] = ~is_calibrator

    results = {}
    for name, mask in subsets.items():
        n_sel = int(np.sum(mask))
        if n_sel < 10:
            print_status(f"  Skipping subset '{name}': only {n_sel} SNe", "WARNING")
            continue

        z_sub = z_all[mask]
        mb_sub = mb_all[mask]
        cov_sub = cov_all[np.ix_(mask, mask)]

        print_status(f"  Fitting subset '{name}': {n_sel} SNe", "INFO")
        fit = fit_models(z_sub, mb_sub, cov_sub)
        results[name] = fit
        print_status(
            f"    Δχ²(TEP-LCDM)={fit['delta_chi2']:+.2f}, eps={fit['bestfit_tep']['epsilon_shear']:.3f}",
            "INFO",
        )

    # Robustness summary. Sign alone is not a robustness metric: a
    # subset with |delta_chi2| ~ 0.1 is statistically consistent with no
    # preference. The classification therefore combines sign
    # persistence with the magnitude of the residual improvement and
    # reports where the full-sample signal is carried.
    delta_chi2_values = [r["delta_chi2"] for r in results.values() if "delta_chi2" in r]
    all_negative = all(d < 0 for d in delta_chi2_values) if delta_chi2_values else False
    majority_negative = (sum(d < 0 for d in delta_chi2_values) / len(delta_chi2_values)) > 0.5 if delta_chi2_values else False

    # Per-subset improvement large enough to matter at MLE level.
    # With one additional fitted parameter, |Delta chi2| ~ 1 is within
    # fit-level noise; a nontrivial residual preference requires
    # delta_chi2 < -2.
    nontrivial = {name: r["delta_chi2"] < -2.0 for name, r in results.items() if "delta_chi2" in r}
    n_nontrivial = sum(nontrivial.values())
    total = len(nontrivial)
    full_dchi2 = results.get("all", {}).get("delta_chi2", None)

    if all_negative and n_nontrivial >= max(1, int(0.75 * total)):
        robustness = "strong"
    elif all_negative:
        robustness = "marginal"
    elif majority_negative:
        robustness = "mixed"
    else:
        robustness = "weak"

    payload = {
        "step": STEP_ID,
        "description": "Pantheon+ subset and redshift-window robustness",
        "status": "completed",
        "subsets_tested": list(results.keys()),
        "results_per_subset": results,
        "robustness_assessment": robustness,
        "all_subsets_prefer_tep": bool(all_negative),
        "majority_prefer_tep": bool(majority_negative),
        "nontrivial_subsets": nontrivial,
        "n_nontrivial": int(n_nontrivial),
        "n_subsets": int(total),
        "full_sample_delta_chi2": full_dchi2,
        "interpretation": (
            f"The M1-vs-LCDM improvement is nontrivial (Delta chi2 < -2) in "
            f"{n_nontrivial} of {total} subsets (full sample Delta chi2 = "
            f"{full_dchi2:+.2f}). Redshift-window and leave-one-out tests "
            f"localise where the full-sample preference is carried; subsets "
            f"with |Delta chi2| ~ 0 are consistent with no preference rather "
            f"than independent confirmation."
            if majority_negative else
            "Subset tests reveal heterogeneity; several subsets favor LCDM, "
            "indicating where the preference originates."
        ),
        "timestamp": int(time.time()),
    }

    write_json(step_json_path(STEP_ID), payload)
    print_status(f"Step {STEP_ID} completed", "SUCCESS")
    return payload


if __name__ == "__main__":
    print(run())
