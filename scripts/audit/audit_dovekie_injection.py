#!/usr/bin/env python3
"""Injection calibration for the DES-Dovekie environmental-kernel result.

The free logistic depth gate gains dchi2 = -16.3 vs LCDM on the real
catalog. This audit asks whether that value is achievable under:

  A) GR-true mocks  -- inject the LCDM best-fit mu plus covariance draw,
     refit all models. Distribution of dchi2(E3-LCDM) measures the
     false-positive rate of the gate on GR-true data.
  B) Env-true mocks -- inject the E3 best-fit model mu plus covariance
     draw, refit. Measures whether the machinery recovers the injected
     gate and what dchi2 a real environmental signal produces.
  C) Generic mass-correlation control -- inject LCDM mu plus a smooth
     mass-dependent offset (power law in logM, amplitude matched to the
     observed residual span). If E3 soaks up an arbitrary smooth mass
     systematic as easily, the -16.3 is not kernel-specific.

Output: results/outputs/audit_dovekie_injection.json
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "steps"))
sys.path.insert(0, str(ROOT / "scripts"))

from audit_dovekie_env_kernel import load_dovekie
from step_03_11_env_modulated_kernel import EnvKernelFit

OUT = ROOT / "results" / "outputs" / "audit_dovekie_injection.json"
N_MOCKS = 24


def fit_suite(fit):
    """Refit LCDM and all env kernels on a mocked mu; returns dchi2 map."""
    r_l = minimize(lambda t: fit.chi2_lcdm(t), [0.35, 0.0],
                   method="Nelder-Mead",
                   options={"xatol": 1e-8, "fatol": 1e-8})
    chi2_l = r_l.fun
    out = {"chi2_lcdm": float(chi2_l)}

    def fit_w(wf, bounds_env):
        all_bounds = [(0.0, 2.0)] + list(bounds_env) + [(-5.0, 5.0)]

        def objective(theta):
            w = wf(tuple(theta[1:1 + len(bounds_env)]))
            if not np.all(np.isfinite(w)) or np.any(w <= 0):
                return 1e10
            return fit.chi2_given_w(theta[0], theta[-1], w)

        rng = np.random.default_rng(7)
        best = None
        for _ in range(10):
            x0 = [rng.uniform(0.0, 2.0)] + \
                 [rng.uniform(*b) for b in bounds_env] + \
                 [rng.uniform(-5, 5)]
            r = minimize(objective, x0, method="L-BFGS-B",
                         bounds=all_bounds)
            if best is None or r.fun < best.fun:
                best = r
        return best

    for name, wf, b in [
        ("E0", fit.w_uniform, []),
        ("E2", fit.w_mass, [(-1.0, 1.0)]),
        ("E3", fit.w_logistic, [(0.0, 3.0), (8.5, 11.5)]),
        ("E4", fit.w_locked, []),
    ]:
        r = fit_w(wf, b)
        out[name] = {"dchi2": float(r.fun - chi2_l),
                     "theta": [float(t) for t in r.x]}
    return out


def main():
    z, mu, mass, is_cal, cov, _ = load_dovekie()
    rng = np.random.default_rng(2026)
    L = np.linalg.cholesky(cov)

    # Real-data reference fits for injection truths
    fit0 = EnvKernelFit(z, mu, mass, is_cal, cov)
    r_l0 = minimize(lambda t: fit0.chi2_lcdm(t), [0.35, 0.0],
                    method="Nelder-Mead",
                    options={"xatol": 1e-8, "fatol": 1e-8})
    Om_l, M_l = r_l0.x
    from core.cosmology import _integrated_distance
    chi_l = _integrated_distance(
        z, lambda zz: 299792.458 / (70 * np.sqrt(
            Om_l * (1 + zz) ** 3 + (1 - Om_l))))
    mu_lcdm = 5 * np.log10(np.maximum(chi_l * (1 + z), 1e-30)) + 25.0 + M_l

    # E3 truth: use the real-data fitted parameters (eps=0.998, a=0.462,
    # m_t=11.5, M) from the previous audit
    eps_e3, a_e3, mt_e3, M_e3 = 0.998, 0.462, 11.5, 0.017
    w_e3 = fit0.w_logistic((a_e3, mt_e3))
    mu_env = fit0.mu_model(eps_e3, w_e3) + M_e3

    # Generic smooth mass-offset truth: LCDM + c*(logM-10.5) offset with
    # amplitude spanning the observed tercile residual gap (~0.06 mag)
    mass_c = mass - 10.5
    mu_massoff = mu_lcdm + 0.03 * mass_c

    results = {"n": int(len(z)), "n_mocks": N_MOCKS,
               "observed_dchi2_E3": -16.31, "truths": {}}

    for label, mu_true in [("GR_true", mu_lcdm),
                           ("env_true_E3", mu_env),
                           ("mass_offset", mu_massoff)]:
        d_e3, d_e2, d_e0, d_e4 = [], [], [], []
        thetas = []
        for i in range(N_MOCKS):
            mu_mock = mu_true + L @ rng.standard_normal(len(z))
            fit = EnvKernelFit(z, mu_mock, mass, is_cal, cov)
            s = fit_suite(fit)
            d_e3.append(s["E3"]["dchi2"]); d_e2.append(s["E2"]["dchi2"])
            d_e0.append(s["E0"]["dchi2"]); d_e4.append(s["E4"]["dchi2"])
            thetas.append(s["E3"]["theta"])
            if (i + 1) % 4 == 0:
                print(f"  {label}: {i+1}/{N_MOCKS} "
                      f"E3 median {np.median(d_e3):+.2f}", flush=True)
        arr = np.array(d_e3)
        results["truths"][label] = {
            "dchi2_E3_median": float(np.median(arr)),
            "dchi2_E3_q16_84": [float(np.percentile(arr, 16)),
                                float(np.percentile(arr, 84))],
            "frac_E3_better_than_obs": float(np.mean(arr <= -16.31)),
            "dchi2_E2_median": float(np.median(d_e2)),
            "dchi2_E0_median": float(np.median(d_e0)),
            "dchi2_E4_median": float(np.median(d_e4)),
            "E3_thetas": thetas,
        }
        print(f"{label}: E3 dchi2 median {np.median(arr):+.2f} "
              f"[{np.percentile(arr,16):+.2f},{np.percentile(arr,84):+.2f}] "
              f"P(<=-16.3)={np.mean(arr<=-16.31):.2f}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Saved {OUT}")


if __name__ == "__main__":
    main()
