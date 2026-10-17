#!/usr/bin/env python3
"""Step 03_11: Pantheon+ Joint Profile Likelihood for Chandrasekhar Luminosity Drift.

Tier: RESEARCH GRADE

Scientific Purpose:
Tests whether a flexible TEP distance-redshift relation d_L(z) can absorb
the intrinsic Chandrasekhar luminosity drift predicted by k-mouflage tracking:
    dm(z) = -7.5 * eps * log10(1+z)
where eps is the interior tracking fraction (eps = 1 implies G_loc ∝ (1+z)^-2
inside source galaxies, as would occur if interior conformal factors tracked
the cosmological clock drift without pinning).

This step evaluates the joint profile likelihood chi^2(eps) over the full
Pantheon+ sample (1,701 SNe Ia) using the complete Cholesky-decomposed
covariance matrix from step 03_01, re-optimizing the TEP cosmological parameters
(line-of-sight shear epsilon_shear_los and magnitude normalization M) at every eps.

Key Results:
- On the fixed-z_T=5 M1 branch, full tracking (eps = 1.0) incurs
  Delta chi2 about 801 under the stated Chandrasekhar luminosity prescription.
- The interpolated 3-sigma upper bound is eps < 0.155 on this branch.
- This SN Ia progenitor constraint does not by itself constrain BBH source
  evolution or constitute a selection-corrected siren-population likelihood.

Dependencies:
- data/pantheon_plus/Pantheon+SH0ES.dat
- data/pantheon_plus/Pantheon+SH0ES.cov
- scripts/steps/step_03_01_three_model_comparison.py (PantheonData, TEPCosmology)

Outputs:
- results/step_03_11_sn_luminosity_drift_profile.json
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Dict

import numpy as np
from scipy.optimize import minimize

# Ensure c0_common and sister steps are importable
STEP_DIR = Path(__file__).resolve().parent
REPO_ROOT = STEP_DIR.parent.parent
sys.path.insert(0, str(STEP_DIR))
sys.path.insert(0, str(REPO_ROOT))

from c0_common import (
    TEPLogger, ensure_dirs, print_status, set_step_logger, write_json
)
from step_03_01_three_model_comparison import PantheonData, TEPCosmology
from scripts.utils.static_metric import StaticCosmology

STEP_ID = "step_03_11_sn_luminosity_drift_profile"


def run() -> None:
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}: Pantheon+ Drift Profile Likelihood", "TITLE")
    ensure_dirs()

    # Load Pantheon+ data
    data = PantheonData()
    if not data.load():
        raise RuntimeError("Failed to load Pantheon+ dataset")

    print_status(f"Loaded {len(data.z)} SNe Ia with covariance {data.cov.shape}", "INFO")

    # 1. Profile likelihood on M1 NoLambda (matter-only, flexible shear, z_T = 5)
    eps_grid = np.linspace(-0.25, 0.40, 14)
    chi2_prof = []
    best_params_prof = []

    print_status("Computing profile likelihood across eps in [-0.25, 0.40]...", "PROCESS")
    for eps_val in eps_grid:
        def obj(params):
            eps_shear, M = params
            if eps_shear < 0:
                return 1e10
            cosmo = TEPCosmology(H0=70.0, Omega_m=1.0, epsilon_T=eps_shear, z_T=5.0)
            mu = cosmo.distance_modulus(data.z) + M - 7.5 * eps_val * np.log10(1.0 + data.z)
            res = data.mb - mu - data.observer_correction
            return -data.gaussian_loglike(res)

        res = minimize(obj, [0.8, -19.32], method="Nelder-Mead")
        cosmo = TEPCosmology(H0=70.0, Omega_m=1.0, epsilon_T=res.x[0], z_T=5.0)
        mu = cosmo.distance_modulus(data.z) + res.x[1] - 7.5 * eps_val * np.log10(1.0 + data.z)
        c2 = float(data.chi2(data.mb - mu - data.observer_correction))
        chi2_prof.append(c2)
        best_params_prof.append({"eps_shear": float(res.x[0]), "M": float(res.x[1])})

    chi2_prof = np.array(chi2_prof)
    dchi2 = chi2_prof - np.min(chi2_prof)
    best_idx = int(np.argmin(chi2_prof))
    eps_best = float(eps_grid[best_idx])

    # Interpolate bounds
    eps_dense = np.linspace(-0.25, 0.40, 1000)
    dchi2_dense = np.interp(eps_dense, eps_grid, dchi2)
    e_1sig_lo = float(eps_dense[np.where((dchi2_dense < 1.0) & (eps_dense <= eps_best))[0][0]])
    e_1sig_hi = float(eps_dense[np.where((dchi2_dense < 1.0) & (eps_dense >= eps_best))[0][-1]])
    e_3sig_hi = float(eps_dense[np.where((dchi2_dense < 9.0) & (eps_dense >= eps_best))[0][-1]])

    print_status(f"Best-fit eps = {eps_best:.3f}, 1-sigma: [{e_1sig_lo:.3f}, {e_1sig_hi:.3f}], 3-sigma: eps < {e_3sig_hi:.3f}", "SUCCESS")

    # 2. Test full tracking eps = 1.0 on M1 NoLambda
    def obj_m1_eps1(params):
        eps_shear, M = params
        if eps_shear < 0:
            return 1e10
        cosmo = TEPCosmology(H0=70.0, Omega_m=1.0, epsilon_T=eps_shear, z_T=5.0)
        mu = cosmo.distance_modulus(data.z) + M - 7.5 * 1.0 * np.log10(1.0 + data.z)
        return -data.gaussian_loglike(data.mb - mu - data.observer_correction)

    res_m1_1 = minimize(obj_m1_eps1, [0.89, -19.32], method="Nelder-Mead")
    cosmo_eps1 = TEPCosmology(H0=70.0, Omega_m=1.0, epsilon_T=res_m1_1.x[0], z_T=5.0)
    mu_eps1 = cosmo_eps1.distance_modulus(data.z) + res_m1_1.x[1] - 7.5 * np.log10(1.0 + data.z)
    chi2_m1_eps1 = float(data.chi2(data.mb - mu_eps1 - data.observer_correction))
    delta_chi2_eps1 = float(chi2_m1_eps1 - np.min(chi2_prof))

    print_status(f"eps = 1.0 test on M1 NoLambda: Delta chi2 = {delta_chi2_eps1:.1f}", "INFO")

    out: Dict = {
        "step": STEP_ID,
        "dataset": "Pantheon+ (1,701 SNe Ia, full covariance matrix)",
        "observer_redshift": "zHEL",
        "distance_redshift": data.redshift_column_used,
        "question": "Does a flexible distance-redshift relation absorb the Chandrasekhar luminosity drift dm(z) = -7.5 eps log10(1+z)?",
        "profile_likelihood_M1": {
            "eps_grid": [float(x) for x in eps_grid],
            "chi2": [float(x) for x in chi2_prof],
            "delta_chi2": [float(x) for x in dchi2],
            "best_fit_eps": eps_best,
            "one_sigma_interval": [e_1sig_lo, e_1sig_hi],
            "three_sigma_upper_bound": e_3sig_hi,
        },
        "full_tracking_test_eps_1": {
            "M1_NoLambda": {
                "chi2": chi2_m1_eps1,
                "delta_chi2_vs_best": delta_chi2_eps1,
                "best_shear_param": float(res_m1_1.x[0]),
                "best_M": float(res_m1_1.x[1]),
                "verdict": f"EXCLUDED at Delta chi2 = {delta_chi2_eps1:.1f} (>> 9, the 3-sigma threshold)"
            }
        },
        "verdict": (
            "Conditional on the fixed-z_T=5 M1 NoLambda distance relation and the "
            "Chandrasekhar luminosity prescription, Pantheon+ constrains the SN Ia "
            "progenitor tracking fraction. Full tracking (eps = 1.0) incurs "
            f"Delta chi2 = {delta_chi2_eps1:.1f} after re-optimizing the "
            "nonnegative line-of-sight shear and magnitude intercept; the "
            f"interpolated upper bound is eps < {e_3sig_hi:.3f} at Delta chi2 = 9. "
            "No BBH population evolution, host-to-source transfer, or GW selection "
            "function is fitted, so this is not a model-independent bound on "
            "standard-siren source masses."
        )
    }

    out_file = REPO_ROOT / "results" / f"{STEP_ID}.json"
    with open(out_file, "w") as f:
        json.dump(out, f, indent=2)
    print_status(f"Saved results to {out_file}", "SUCCESS")


if __name__ == "__main__":
    run()
