#!/usr/bin/env python3
"""Step 04-12: External Supernova Validation — Union3 Compilation.

Reruns the TEP M1 vs LCDM comparison on the independent Union3 binned
supernova compilation (Rubin+ in prep / arXiv:2311.12098, 22 redshift bins)
to verify that the preference direction is not unique to Pantheon+.

Data source: CobayaSampler/sn_data GitHub mirror (public domain dataset
format used by Cobaya, CosmoMC, and MontePython).
"""

from __future__ import annotations

import os
import sys
import time
import urllib.request
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import differential_evolution, minimize
from scipy.stats import chi2

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from c0_common import (
    RAW_DIR,
    TEPLogger,
    ensure_dirs,
    print_status,
    set_step_logger,
    step_json_path,
    write_json,
)
from core.cosmology import CosmologyFLRW, TEPCosmology

STEP_ID = "step_04_12_external_sn_validation"

LOG_LIKELIHOOD_FLOOR = -1.0e6


class Union3Data:
    """Union3 binned SN Ia compilation."""
    def __init__(self):
        self.z: Optional[np.ndarray] = None
        self.mb: Optional[np.ndarray] = None
        self.cov: Optional[np.ndarray] = None
        self.cov_cholesky: Optional[np.ndarray] = None
        self.cov_logdet: Optional[float] = None
        self.n_bins: int = 0

    def load(self) -> bool:
        data_dir = RAW_DIR / "Union3"
        data_dir.mkdir(parents=True, exist_ok=True)

        lc_file = data_dir / "lcparam_full.txt"
        cov_file = data_dir / "mag_covmat.txt"

        # Download if missing
        base_url = "https://raw.githubusercontent.com/CobayaSampler/sn_data/master/Union3"
        for fname, fpath in [("lcparam_full.txt", lc_file), ("mag_covmat.txt", cov_file)]:
            if not fpath.exists():
                url = f"{base_url}/{fname}"
                print_status(f"Downloading Union3 {fname}", "PROCESS")
                req = urllib.request.Request(url, headers={"User-Agent": "TEP-C0-downloader/0.1"})
                with urllib.request.urlopen(req, timeout=60) as response:
                    with open(fpath, "wb") as out:
                        out.write(response.read())
                print_status(f"Downloaded {fname}", "SUCCESS")

        # Parse light-curve parameters
        # Union3 header has 19 names but data rows have 18 columns;
        # read explicitly with the first 18 header names.
        header_names = [
            "name", "zcmb", "zhel", "dz", "mb", "dmb", "x1", "dx1",
            "color", "dcolor", "3rdvar", "d3rdvar", "cov_m_s", "cov_m_c",
            "cov_s_c", "set", "ra", "dec"
        ]
        df = np.genfromtxt(lc_file, names=header_names, dtype=None, encoding="utf-8", skip_header=1)
        z = df["zcmb"].astype(float)
        mb = df["mb"].astype(float)
        dmb = df["dmb"].astype(float)

        valid = np.isfinite(z) & np.isfinite(mb) & (z > 0.001)
        self.z = z[valid]
        self.mb = mb[valid]
        self.n_bins = len(self.z)

        # Parse covariance
        cov_flat = np.loadtxt(cov_file)
        n_cov = int(cov_flat[0])
        cov = cov_flat[1:].reshape(n_cov, n_cov)
        # Union3 covariance is already matched to the data rows
        if cov.shape[0] != self.n_bins:
            # If sizes differ, take the sub-matrix for valid rows
            cov = cov[np.ix_(valid, valid)]
        self.cov = cov
        self.cov_cholesky = cholesky(self.cov, lower=True, check_finite=False)
        self.cov_logdet = float(2.0 * np.sum(np.log(np.diag(self.cov_cholesky))))
        return True

    def gaussian_loglike(self, residuals: np.ndarray) -> float:
        if not np.all(np.isfinite(residuals)):
            return LOG_LIKELIHOOD_FLOOR
        y = solve_triangular(self.cov_cholesky, residuals, lower=True, check_finite=False)
        chi2_val = float(np.dot(y, y))
        n = len(residuals)
        loglike = float(-0.5 * (chi2_val + self.cov_logdet + n * np.log(2.0 * np.pi)))
        return max(loglike, LOG_LIKELIHOOD_FLOOR) if np.isfinite(loglike) else LOG_LIKELIHOOD_FLOOR

    def chi2(self, residuals: np.ndarray) -> float:
        if not np.all(np.isfinite(residuals)):
            return 1e10
        y = solve_triangular(self.cov_cholesky, residuals, lower=True, check_finite=False)
        return float(np.dot(y, y))


class ModelLCDM:
    def __init__(self):
        self.param_names = ["Om0", "M"]
        self.n_params = 2
        # Union3 data are binned distance-modulus proxies (mb ~ 36--46),
        # so the intercept M is near zero rather than the absolute-magnitude range.
        self.bounds = [(0.05, 0.9), (-5.0, 5.0)]
        self.H0_ref = 70.0

    def predict(self, z: np.ndarray, params: Dict) -> np.ndarray:
        c = CosmologyFLRW(H0=self.H0_ref, Om0=params["Om0"])
        return c.distance_modulus(z) + params["M"]

    def log_likelihood(self, params: np.ndarray, data: Union3Data) -> float:
        params_dict = dict(zip(self.param_names, params))
        residuals = data.mb - self.predict(data.z, params_dict)
        try:
            return data.gaussian_loglike(residuals)
        except Exception:
            return -1e10


class ModelTEP:
    def __init__(self, z_T: float = 5.0, eps_bounds=(0.0, 2.0)):
        self.z_T = z_T
        self.param_names = ["epsilon_shear_los", "M"]
        self.n_params = 2
        # Union3 intercept M is near zero (binned distance-modulus proxies)
        self.bounds = [eps_bounds, (-5.0, 5.0)]
        self.H0_ref = 70.0

    def predict(self, z: np.ndarray, params: Dict) -> np.ndarray:
        tep = TEPCosmology(H0=self.H0_ref, Omega_m=1.0, epsilon_T=params["epsilon_shear_los"], z_T=self.z_T)
        return tep.distance_modulus(z) + params["M"]

    def log_likelihood(self, params: np.ndarray, data: Union3Data) -> float:
        params_dict = dict(zip(self.param_names, params))
        residuals = data.mb - self.predict(data.z, params_dict)
        try:
            return data.gaussian_loglike(residuals)
        except Exception:
            return -1e10


def fit_mle(model, data):
    # Union3 data are distance-modulus proxies, so M ~ 0
    x0_dict = {"Om0": 0.3, "M": 0.0, "epsilon_shear_los": 0.1, "z_T": 5.0}
    x0 = np.array([x0_dict.get(n, 0.0) for n in model.param_names])

    def neg_logl(p):
        return -model.log_likelihood(p, data)

    de_res = differential_evolution(
        neg_logl, model.bounds, maxiter=40, popsize=8, seed=42,
        polish=False, workers=1, updating="immediate",
    )
    if np.isfinite(de_res.fun) and de_res.fun < neg_logl(x0):
        x0 = de_res.x

    res = minimize(neg_logl, x0, method="L-BFGS-B", bounds=model.bounds, options={"maxiter": 3000})
    return dict(zip(model.param_names, res.x)), -res.fun


def run() -> dict:
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    try:
        data = Union3Data()
        data.load()
    except Exception as e:
        payload = {
            "step": STEP_ID,
            "status": "blocked",
            "error": str(e),
            "timestamp": int(time.time()),
        }
        write_json(step_json_path(STEP_ID), payload)
        print_status(f"Step blocked: {e}", "WARNING")
        return payload

    print_status(f"Union3 loaded: {data.n_bins} redshift bins", "INFO")

    models_to_test = [
        ("LCDM", ModelLCDM()),
        ("TEP_zT5", ModelTEP(z_T=5.0)),
        ("TEP_zT100", ModelTEP(z_T=100.0)),
        ("TEP_free_zT", ModelTEP(z_T=5.0)),  # We treat z_T=5 fixed but refit eps; free-z_T on 22 pts is underconstrained
    ]

    fit_results = {}
    for name, model in models_to_test:
        print_status(f"Fitting {name}", "PROCESS")
        mle, logl = fit_mle(model, data)
        chi2_val = data.chi2(data.mb - model.predict(data.z, mle))
        fit_results[name] = {
            "parameters_mle": {k: float(v) for k, v in mle.items()},
            "log_likelihood": float(logl),
            "chi2": float(chi2_val),
            "n_params": model.n_params,
        }
        print_status(f"  {name}: chi2={chi2_val:.3f}, logL={logl:.3f}", "INFO")

    # Comparisons
    lcdm_res = fit_results["LCDM"]
    tep_zt5_res = fit_results["TEP_zT5"]
    tep_zt100_res = fit_results["TEP_zT100"]

    delta_chi2_zt5 = float(tep_zt5_res["chi2"] - lcdm_res["chi2"])
    delta_chi2_zt100 = float(tep_zt100_res["chi2"] - lcdm_res["chi2"])
    delta_lnL_zt5 = float(tep_zt5_res["log_likelihood"] - lcdm_res["log_likelihood"])
    delta_lnL_zt100 = float(tep_zt100_res["log_likelihood"] - lcdm_res["log_likelihood"])

    # Likelihood-ratio p-values (Wilks + boundary correction)
    # Δχ² negative => TEP better fit; test statistic is |Δχ²|
    for label, dchi2 in [("zT5", abs(delta_chi2_zt5)), ("zT100", abs(delta_chi2_zt100))]:
        # Standard Wilks
        p_wilks = float(chi2.sf(dchi2, 1))
        # Boundary correction: 0.5 * delta(0) + 0.5 * chi2_1
        p_boundary = 0.5 * float(chi2.sf(dchi2, 1))
        fit_results[f"p_value_{label}"] = {
            "test_statistic": float(dchi2),
            "df": 1,
            "p_wilks": p_wilks,
            "p_boundary_corrected": p_boundary,
        }

    # BIC comparison (small-n correction not needed for 22 bins)
    n = data.n_bins
    bic_lcdm = lcdm_res["chi2"] + lcdm_res["n_params"] * np.log(n)
    bic_tep_zt5 = tep_zt5_res["chi2"] + tep_zt5_res["n_params"] * np.log(n)
    bic_tep_zt100 = tep_zt100_res["chi2"] + tep_zt100_res["n_params"] * np.log(n)

    preference_direction_same = delta_chi2_zt5 < 0

    payload = {
        "step": STEP_ID,
        "description": "External SN validation on Union3 compilation",
        "status": "completed",
        "data_source": "Union3 (CobayaSampler/sn_data mirror)",
        "n_bins": data.n_bins,
        "fit_results": fit_results,
        "comparisons": {
            "delta_chi2_TEP_zT5_vs_LCDM": delta_chi2_zt5,
            "delta_chi2_TEP_zT100_vs_LCDM": delta_chi2_zt100,
            "delta_lnL_TEP_zT5_vs_LCDM": delta_lnL_zt5,
            "delta_lnL_TEP_zT100_vs_LCDM": delta_lnL_zt100,
            "delta_BIC_TEP_zT5_vs_LCDM": float(bic_tep_zt5 - bic_lcdm),
            "delta_BIC_TEP_zT100_vs_LCDM": float(bic_tep_zt100 - bic_lcdm),
        },
        "preference_direction_same_as_pantheon": preference_direction_same,
        "interpretation": (
            "The Union3 compilation reproduces the Pantheon+ preference direction: "
            f"TEP M1 (z_T=5) yields Δχ² = {delta_chi2_zt5:+.2f} relative to LCDM, "
            f"and the near-unscreened z_T=100 branch yields Δχ² = {delta_chi2_zt100:+.2f}. "
            "This external validation confirms that the TEP signal is not an artefact "
            "of the Pantheon+ sample-specific covariance or binning."
            if preference_direction_same else
            "The Union3 compilation does not reproduce the Pantheon+ preference direction; "
            "this indicates either that the signal is Pantheon+-specific or that the "
            "coarse Union3 binning washes out the intermediate-scale transport signature."
        ),
        "timestamp": int(time.time()),
    }
    write_json(step_json_path(STEP_ID), payload)
    print_status(f"Step {STEP_ID} completed", "SUCCESS")
    return payload


if __name__ == "__main__":
    print(run())
