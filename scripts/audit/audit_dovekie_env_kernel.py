#!/usr/bin/env python3
"""Audit: environmental transport kernel on DES-Dovekie -- independent
catalog WITH host masses.

DES-SN5YR's Hubble diagram alone cannot test the screened kernel (no
host-mass column in the HD file), and the Pantheon+ z>0.023 subset
actively disfavoured the locked screening weights (dchi2 = +41.6).
The DES-Dovekie release (arXiv:2511.07517 -- the DES-SN5YR+
Foundation+Low-z re-calibration with its own standardization) ships
HOST_LOGMASS per object plus the official STAT+SYS inverse covariance,
so the environmental prediction can be scored on a genuinely
independent catalog.

Models (identical EnvKernelFit machinery and constants to step_03_11):
  LCDM        (Om0, M)
  E0 uniform  w_i = 1                              (eps, M)
  E1 binary   w_i = exp(a h_i), h=+/-0.5 @ logM=10 (eps, a, M)
  E3 logistic w_i = exp[-a L_i]                    (eps, a, m_t, M)
  E4 locked   w_i = S_A(sigma_i)/norm              (eps, M)

Also scored: per-host-mass-tercile residual mean under the uniform-kernel
MLE -- does the wing residual structure track host depth on an
independent catalog?

Output: results/outputs/audit_dovekie_env_kernel.json
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "steps"))

from step_03_11_env_modulated_kernel import (
    EnvKernelFit, MASS_SPLIT, MASS_SENTINEL,
    SIGMA_REF, M_REF, SIGMA_T_SA, FJ_SLOPE,
)

OUT = ROOT / "results" / "outputs" / "audit_dovekie_env_kernel.json"


def load_dovekie():
    hd = pd.read_csv(ROOT / "data/raw/DES_Dovekie/DES-Dovekie_HD.csv",
                     comment="#", sep=r"\s+")
    md = pd.read_csv(ROOT / "data/raw/DES_Dovekie/DES-Dovekie_Metadata.csv",
                     comment="#", sep=r"\s+", low_memory=False)
    m = hd.merge(md[["CID", "HOST_LOGMASS"]], on="CID", how="left")
    z = m["zHD"].values.astype(float)
    mu = m["MU"].values.astype(float)
    mass = m["HOST_LOGMASS"].values.astype(float)
    ok = (np.isfinite(z) & np.isfinite(mu) & np.isfinite(mass)
          & (mass > 5.0) & (z > 0.001) & (z < 3.0))
    # official STAT+SYS inverse covariance (upper-tri packed)
    d = np.load(ROOT / "data/raw/DES_Dovekie/STAT+SYS.npz")
    n = int(d[d.files[0]][0])
    inv = np.zeros((n, n))
    inv[np.triu_indices(n)] = d[d.files[1]]
    il = np.tril_indices(n, -1)
    inv[il] = inv.T[il]
    assert n == len(z), (n, len(z))
    cov = np.linalg.inv(inv)
    return z[ok], mu[ok], mass[ok], np.zeros(ok.sum(), bool), \
        cov[np.ix_(ok, ok)], int((~ok).sum())


def main():
    z, mu, mass, is_cal, cov, ndrop = load_dovekie()
    print(f"DES-Dovekie: n={len(z)} (dropped {ndrop}), "
          f"z {z.min():.3f}-{z.max():.3f}, "
          f"mass {mass.min():.1f}-{mass.max():.1f}")
    fit = EnvKernelFit(z, mu, mass, is_cal, cov)

    res = {"n": int(len(z)), "n_dropped": ndrop,
           "n_lowmass": int((mass < MASS_SPLIT).sum()),
           "w_sa_range": [float(fit.w_sa.min()), float(fit.w_sa.max())],
           "frac_w_lt_0p5": float((fit.w_sa < 0.5).mean()),
           "models": {}}

    # Dovekie MU is a distance modulus (mag ~33-45), not an apparent
    # magnitude; the intercept M must sit near 0, not at M_B ~ -19.3.
    # Local multistart with the correct M bounds.
    def fit_mu(w_func, bounds_env, eps_bounds=(0.0, 2.0),
               m_bounds=(-5.0, 5.0)):
        all_bounds = [eps_bounds] + list(bounds_env) + [m_bounds]

        def objective(theta):
            eps = theta[0]
            ep = tuple(theta[1:1 + len(bounds_env)])
            M = theta[-1]
            w = w_func(ep)
            if not np.all(np.isfinite(w)) or np.any(w <= 0):
                return 1e10
            return fit.chi2_given_w(eps, M, w)

        best = None
        rng = np.random.default_rng(42)
        starts = [[0.87] + [np.mean(b) for b in bounds_env] + [0.0]]
        for _ in range(15):
            starts.append([
                rng.uniform(*eps_bounds),
                *[rng.uniform(*b) for b in bounds_env],
                rng.uniform(*m_bounds),
            ])
        for x0 in starts:
            r = minimize(objective, x0, method="L-BFGS-B",
                         bounds=all_bounds)
            if best is None or r.fun < best.fun:
                best = r
        return best

    # LCDM reference (chi2 MLE; M free, lands near 0 for mu data)
    r_l = minimize(lambda t: fit.chi2_lcdm(t), [0.35, 0.0],
                   method="Nelder-Mead",
                   options={"xatol": 1e-8, "fatol": 1e-8})
    chi2_lcdm = r_l.fun
    res["models"]["LCDM"] = {"chi2": float(chi2_lcdm),
                             "mle": {"Om0": float(r_l.x[0]),
                                     "M": float(r_l.x[1])}}

    for name, wf, b in [
        ("E0_uniform", fit.w_uniform, []),
        ("E1_binary_mass_gate", fit.w_binary, [(-2.0, 2.0)]),
        ("E2_continuous_mass_gate", fit.w_mass, [(-1.0, 1.0)]),
        ("E3_logistic_depth_gate", fit.w_logistic, [(0.0, 3.0), (8.5, 11.5)]),
        ("E4_locked_SA_prediction", fit.w_locked, []),
    ]:
        r = fit_mu(wf, b)
        res["models"][name] = {"chi2": float(r.fun),
                               "dchi2_vs_LCDM": float(r.fun - chi2_lcdm),
                               "theta": [float(t) for t in r.x]}
        print(f"  {name}: chi2={r.fun:.2f} d={r.fun - chi2_lcdm:+.2f} "
              f"theta={np.round(r.x,3)}")

    # residual-vs-host-mass structure under the uniform MLE
    eps_hat = res["models"]["E0_uniform"]["theta"][0]
    M_hat = res["models"]["E0_uniform"]["theta"][-1]
    resid = mu - (fit.mu_model(eps_hat, np.ones(fit.n)) + M_hat)
    from core.cosmology import _integrated_distance
    Om, Ml = r_l.x[0], r_l.x[1]
    chi_l = _integrated_distance(
        z, lambda zz: 299792.458 / (70 * np.sqrt(
            Om * (1 + zz) ** 3 + (1 - Om))))
    resid_l = mu - (5 * np.log10(np.maximum(
        chi_l * (1 + z), 1e-30)) + 25.0 + Ml)
    q = np.quantile(mass, [0, 1 / 3, 2 / 3, 1])
    terc = {}
    for i, lab in enumerate(["low", "mid", "high"]):
        mk = (mass >= q[i]) & (mass < q[i + 1] if i < 2 else mass <= q[i + 1])
        terc[lab] = {"n": int(mk.sum()),
                     "mass_range": [float(q[i]), float(q[i + 1])],
                     "mean_resid_TEP": float(resid[mk].mean()),
                     "mean_resid_LCDM": float(resid_l[mk].mean())}
    res["resid_by_mass_tercile"] = terc
    print("  resid by mass tercile:",
          {k: round(v["mean_resid_TEP"], 4) for k, v in terc.items()})

    # Weight-vs-mass correlations backing the manuscript's sign-mismatch
    # statement: the locked S_A(sigma_FJ) map rises with host mass while
    # the fitted E3 logistic gate falls with it.
    w_e3_fit = fit.w_logistic(tuple(res["models"]["E3_logistic_depth_gate"]["theta"][1:3]))
    res["weight_mass_correlations"] = {
        "corr_w_locked_sa_vs_mass": float(np.corrcoef(fit.w_sa, mass)[0, 1]),
        "corr_w_e3_fitted_vs_mass": float(np.corrcoef(w_e3_fit, mass)[0, 1]),
        "frac_w_sa_gt_0p9": float((fit.w_sa > 0.9).mean()),
        "frac_w_sa_gt_0p99": float((fit.w_sa > 0.99).mean()),
    }
    print("  corr w_sa-mass: %.3f | corr w_E3-mass: %.3f | frac w>0.9: %.3f"
          % (res["weight_mass_correlations"]["corr_w_locked_sa_vs_mass"],
             res["weight_mass_correlations"]["corr_w_e3_fitted_vs_mass"],
             res["weight_mass_correlations"]["frac_w_sa_gt_0p9"]))

    res["interpretation"] = (
        "Independent-catalog environmental test. If the locked screening "
        "weights again lose to the uniform kernel, the Pantheon+ wing "
        "result is confirmed as catalog-independent structure; if the "
        "mass gates gain here, the wing's adverse reading is "
        "sample-dependent. Residual-by-mass-tercile records the "
        "environmental direction directly.")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=2)
    print(f"Saved {OUT}")


if __name__ == "__main__":
    main()
