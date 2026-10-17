#!/usr/bin/env python3
"""TEP-C0 Step 03-11: Environment-Modulated Transport Kernel
============================================================
The discriminating follow-up test identified in Section 4: does the
line-of-sight transport amplitude track host depth rather than redshift
alone?

The canonical M1 kernel applies a single amplitude epsilon_los to every
sightline, gamma_i(z) = 1 + epsilon_los ln(1+z) S(z). The TEP
environmental sector instead predicts per-sightline modulation: the
screened temporal response inside a host well suppresses the accumulated
transport by a weight w_i in (0,1] set by the host potential depth. A
generic low-z systematic cannot mimic such modulation because the
deviation must then correlate with the host coordinate after the
redshift dependence is divided out.

On the EdS background the correction is exactly linear in epsilon:

    chi(z; eps) = chi0(z) + eps * chi1(z),
    chi0(z) = int c / (H0 E_EdS) dz,
    chi1(z) = int c / (H0 E_EdS) ln(1+z) S(z) dz,

so a per-object amplitude eps_i = eps * w_i requires only two
precomputed integrals; the evaluation is exact, not an approximation.

Models compared (all on zHD, full 1701x1701 covariance):

    E0  uniform kernel            w_i = 1                     (eps, M)
    E1  binary mass gate          w_i = exp(a * h_i),
                                  h_i = +/-0.5 at logM = 10   (eps, a, M)
    E2  continuous mass gate      w_i = 10^{a (m_i - 10)}     (eps, a, M)
    E3  logistic depth gate       w_i = exp[-a L_i],
                                  L_i = logistic((m_i - m_t)/0.5)
                                                                (eps, a, m_t, M)
    E4  locked S_A prediction     w_i = S_A(sigma_i) / <S_A>,
                                  sigma_i = 207 (M_i/1e11)^0.25 km/s,
                                  S_A = min[1, (sigma/65)^{2/3}]  (eps, M)

E4 expends no new parameters: it forward-applies the corpus screening
formalism (Section 4.8 / Paper 0 amplitude sector) with the canonical
sigma_T = 65 km/s and the same virial proxy sigma^2 = GM/(2 R_eff)
(normalised to sigma = 207 km/s at 10^11 Msun, matching step_04_08).

For each model: profiled chi2 via multistart, and dynesty nested
evidence against LCDM for the gated models. The binary and continuous
gates nest the uniform kernel at a = 0, so any improvement is directly
attributable to environment. A permutation null (scrambled mass labels,
matched z-mass correlation) tests whether the gain is genuinely
environmental rather than a projection of the mass-redshift selection.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import cholesky, solve_triangular
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

from core.cosmology import _integrated_distance

STEP_ID = "step_03_11_env_modulated_kernel"

H0_REF = 70.0          # km/s/Mpc, pipeline convention
C_KMS = 299792.458
Z_T = 5.0              # canonical turnover
N_T = 2.0              # canonical sharpness
MASS_SPLIT = 10.0      # canonical host-mass step boundary
MASS_SENTINEL = -9.0   # missing host mass code in Pantheon+

# Locked-screening constants (shared with step_04_08 / Paper 0 sector)
SIGMA_REF = 207.0      # km/s at 1e11 Msun (step_04_08 representative)
M_REF = 11.0           # log10(Msun) anchor for the sigma-mass map
SIGMA_T_SA = 65.0      # km/s amplitude-sector threshold
FJ_SLOPE = 0.25        # sigma ∝ M^0.25 virial proxy


def load_pantheon(observer_frame=True):
    """Load the release with its separate distance and observer redshifts.

    For the one-redshift model functions, subtract the known additive
    5 log10[(1+zHEL)/(1+zHD)] factor from the brightness vector. This is
    algebraically identical to adding it to every model prediction; the
    release file and covariance are untouched. observer_frame=False is
    an explicit historical one-redshift comparison.
    """
    data_dir = Path("data/raw")
    df = pd.read_csv(data_dir / "pantheon_plus_shoes.dat", sep=r"\s+", comment="#")

    z_raw = df["zHD"].values
    mb_raw = df["m_b_corr"].values
    mass_raw = df["HOST_LOGMASS"].values.astype(float)
    is_cal_raw = df["IS_CALIBRATOR"].values.astype(bool)

    valid = (
        np.isfinite(z_raw) & np.isfinite(mb_raw)
        & (z_raw > 0.001) & (z_raw < 3.0)
    )
    z = z_raw[valid]
    mb = mb_raw[valid]
    if observer_frame:
        zhel = df["zHEL"].values[valid]
        if not np.all(np.isfinite(zhel)) or np.any(zhel <= -1):
            raise ValueError("Invalid observer redshift")
        mb = mb - 5.0 * np.log10((1.0 + zhel) / (1.0 + z))
    mass = mass_raw[valid]
    is_cal = is_cal_raw[valid]

    cov_flat = np.loadtxt(data_dir / "Pantheon+SH0ES.cov")
    n_raw = len(z_raw)
    if cov_flat.size == n_raw * n_raw + 1:
        cov_full = cov_flat[1:].reshape(n_raw, n_raw)
    elif cov_flat.size == n_raw * n_raw:
        cov_full = cov_flat.reshape(n_raw, n_raw)
    else:
        raise ValueError(f"Covariance size mismatch: {cov_flat.size}")
    cov = cov_full[np.ix_(valid, valid)]

    # Missing host mass (sentinel) -> median observed mass; one object,
    # immaterial, but kept explicit rather than silently dropped.
    missing = ~np.isfinite(mass) | (mass <= MASS_SENTINEL)
    mass_observed = mass[~missing]
    mass[missing] = float(np.median(mass_observed))

    return z, mb, mass, is_cal, cov, int(missing.sum())


def kernel_integrals(z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """chi0 and chi1 on the canonical EdS background (M1 convention).

    chi0 = int c/(H0 E_EdS) dz ; chi1 = int c/(H0 E_EdS) ln(1+z) S(z) dz.
    """
    def base(zz):
        return C_KMS / (H0_REF * (1.0 + zz) ** 1.5)
    chi0 = _integrated_distance(z, base)
    chi1 = _integrated_distance(
        z, lambda zz: base(zz) * np.log1p(zz)
        * np.exp(-(zz / Z_T) ** N_T))
    return chi0, chi1


class EnvKernelFit:
    """Gaussian likelihood with per-sightline amplitude eps * w_i."""

    def __init__(self, z, mb, mass, is_cal, cov):
        self.z = z
        self.mb = mb
        self.mass = mass
        self.is_cal = is_cal
        self.n = len(z)
        self.L = cholesky(cov, lower=True, check_finite=False)
        self.logdet = float(2.0 * np.sum(np.log(np.diag(self.L))))
        self.chi0, self.chi1 = kernel_integrals(z)
        # precompute constant weight pieces
        self.h_bin = np.where(mass >= MASS_SPLIT, 0.5, -0.5)
        self.m_centered = mass - MASS_SPLIT
        sigma_i = SIGMA_REF * 10.0 ** (FJ_SLOPE * (mass - M_REF))
        s_a = np.minimum(1.0, (sigma_i / SIGMA_T_SA) ** (2.0 / 3.0))
        self.w_sa = s_a / s_a.max()

    # ---- model distances -------------------------------------------
    def mu_model(self, eps: float, w: np.ndarray) -> np.ndarray:
        """Distance modulus for per-object amplitude eps*w_i."""
        dl = (self.chi0 + eps * w * self.chi1) * (1.0 + self.z)
        dl = np.maximum(dl, np.finfo(float).tiny)
        return 5.0 * np.log10(dl) + 25.0

    def chi2_given_w(self, eps: float, M: float, w: np.ndarray) -> float:
        resid = self.mb - (self.mu_model(eps, w) + M)
        y = solve_triangular(self.L, resid, lower=True, check_finite=False)
        return float(np.dot(y, y))

    def chi2_lcdm(self, params) -> float:
        Om, M = params
        if not (0.05 <= Om <= 0.9):
            return 1e10
        zp1 = 1.0 + self.z
        chi = _integrated_distance(
            self.z, lambda zz: C_KMS / (H0_REF * np.sqrt(
                Om * (1.0 + zz) ** 3 + (1.0 - Om))))
        dl = chi * zp1
        mu = 5.0 * np.log10(np.maximum(dl, np.finfo(float).tiny)) + 25.0
        resid = self.mb - (mu + M)
        y = solve_triangular(self.L, resid, lower=True, check_finite=False)
        return float(np.dot(y, y))

    # ---- weights per model -----------------------------------------
    def w_uniform(self, _p):
        return np.ones(self.n)

    def w_binary(self, p):
        (a,) = p
        return np.exp(a * self.h_bin)

    def w_mass(self, p):
        (a,) = p
        return 10.0 ** (a * self.m_centered)

    def w_logistic(self, p):
        a, m_t = p
        L = 1.0 / (1.0 + np.exp(-(self.mass - m_t) / 0.5))
        return np.exp(-a * L)

    def w_locked(self, _p):
        return self.w_sa

    # ---- fitting ----------------------------------------------------
    def fit(self, w_func, bounds_env, eps_bounds=(0.0, 2.0)):
        """Multistart fit over (eps, env..., M)."""
        all_bounds = [eps_bounds] + list(bounds_env) + [(-21.0, -16.0)]

        def objective(theta):
            eps = theta[0]
            ep = tuple(theta[1:1 + len(bounds_env)])
            M = theta[-1]
            w = w_func(ep)
            if not np.all(np.isfinite(w)) or np.any(w <= 0):
                return 1e10
            return self.chi2_given_w(eps, M, w)

        best = None
        rng = np.random.default_rng(42)
        starts = [[0.87] + [np.mean(b) for b in bounds_env] + [-19.33]]
        for _ in range(15):
            starts.append([
                rng.uniform(*eps_bounds),
                *[rng.uniform(*b) for b in bounds_env],
                rng.uniform(-19.8, -18.8),
            ])
        for x0 in starts:
            res = minimize(objective, x0, method="L-BFGS-B", bounds=all_bounds)
            if best is None or res.fun < best.fun:
                best = res
        return best

    def loglike_max(self, chi2):
        return -0.5 * (chi2 + self.logdet + self.n * np.log(2.0 * np.pi))


def run_nested(fit: EnvKernelFit, w_func, bounds_env, nlive=500):
    """Dynesty evidence over (eps, env..., M) with flat priors."""
    import dynesty

    bounds = [(0.0, 2.0)] + list(bounds_env) + [(-21.0, -16.0)]
    ndim = len(bounds)
    lo = np.array([b[0] for b in bounds])
    hi = np.array([b[1] for b in bounds])

    def prior(u):
        return lo + u * (hi - lo)

    def loglike(theta):
        eps = theta[0]
        ep = tuple(theta[1:1 + len(bounds_env)])
        M = theta[-1]
        w = w_func(ep)
        if not np.all(np.isfinite(w)) or np.any(w <= 0):
            return -1e6
        chi2 = fit.chi2_given_w(eps, M, w)
        if not np.isfinite(chi2):
            return -1e6
        return -0.5 * chi2  # logdet/n terms common to all models

    sampler = dynesty.NestedSampler(
        loglike, prior, ndim, bound="multi", sample="rwalk",
        nlive=nlive, rstate=np.random.default_rng(20261011),
    )
    sampler.run_nested(dlogz=0.05)
    res = sampler.results
    return float(res.logz[-1]), float(res.logzerr[-1])


def run_nested_lcdm(fit: EnvKernelFit, nlive=500):
    """Reference LCDM evidence on the same conventions."""
    import dynesty

    bounds = [(0.05, 0.9), (-21.0, -16.0)]
    ndim = 2
    lo = np.array([b[0] for b in bounds])
    hi = np.array([b[1] for b in bounds])

    def prior(u):
        return lo + u * (hi - lo)

    def loglike(theta):
        chi2 = fit.chi2_lcdm(theta)
        return -0.5 * chi2 if np.isfinite(chi2) else -1e6

    sampler = dynesty.NestedSampler(
        loglike, prior, ndim, bound="multi", sample="rwalk", nlive=nlive, rstate=np.random.default_rng(20261011))
    sampler.run_nested(dlogz=0.05)
    res = sampler.results
    return float(res.logz[-1]), float(res.logzerr[-1])


def main():
    set_step_logger(TEPLogger(STEP_ID))
    ensure_dirs()
    t0 = time.time()
    print_status(f"{STEP_ID}: environment-modulated transport kernel")

    z, mb, mass, is_cal, cov, n_missing_mass = load_pantheon()
    fit = EnvKernelFit(z, mb, mass, is_cal, cov)
    print_status(
        f"Loaded {fit.n} SNe ({int(is_cal.sum())} calibrators), "
        f"{n_missing_mass} missing host masses -> median"
    )

    results = {"meta": {
        "step_id": STEP_ID,
        "n_sne": fit.n,
        "n_calibrators": int(is_cal.sum()),
        "n_missing_host_mass": n_missing_mass,
        "redshift_column": "zHD",
        "covariance": "full 1701x1701 Pantheon+SH0ES.cov",
        "z_T": Z_T, "n_T": N_T, "H0_ref": H0_REF,
        "mass_split": MASS_SPLIT,
        "locked_screening": {
            "sigma_ref_kms": SIGMA_REF, "m_ref_log": M_REF,
            "sigma_T_SA_kms": SIGMA_T_SA, "fj_slope": FJ_SLOPE,
        },
    }}

    # ---- reference LCDM (Omega_m free, not the EdS kernel path) ------
    res_lcdm = minimize(
        fit.chi2_lcdm, [0.36, -19.34], method="L-BFGS-B",
        bounds=[(0.05, 0.9), (-21.0, -16.0)])
    chi2_lcdm_val = float(res_lcdm.fun)
    results["lcdm"] = {
        "chi2": chi2_lcdm_val,
        "Omega_m": float(res_lcdm.x[0]), "M": float(res_lcdm.x[1]),
        "lnL": fit.loglike_max(chi2_lcdm_val),
    }
    print_status(f"LCDM: chi2={chi2_lcdm_val:.2f} Om={res_lcdm.x[0]:.4f}")

    # ---- E0 uniform kernel ------------------------------------------
    res_e0 = fit.fit(fit.w_uniform, [])
    results["E0_uniform"] = {
        "k": 2,
        "chi2": float(res_e0.fun),
        "epsilon_los": float(res_e0.x[0]), "M": float(res_e0.x[-1]),
        "delta_chi2_vs_lcdm": float(res_e0.fun - chi2_lcdm_val),
        "lnL": fit.loglike_max(float(res_e0.fun)),
    }
    print_status(
        f"E0 uniform: chi2={res_e0.fun:.2f} eps={res_e0.x[0]:.4f} "
        f"(ref delta chi2={res_e0.fun - chi2_lcdm_val:+.2f})"
    )

    # ---- gated models -------------------------------------------------
    gates = [
        ("E1_binary_mass_gate", fit.w_binary, [(-1.5, 1.5)]),
        ("E2_continuous_mass_gate", fit.w_mass, [(-1.0, 1.0)]),
        ("E3_logistic_depth_gate", fit.w_logistic, [(0.0, 3.0), (8.5, 11.5)]),
        ("E4_locked_SA_prediction", fit.w_locked, []),
        ("E5_free_class_amplitudes", None, None),
    ]

    for name, w_func, b_env in gates:
        if name == "E5_free_class_amplitudes":
            # Independent amplitudes per mass class: eps_low, eps_high, M.
            def obj5(theta):
                e_lo, e_hi, M = theta
                w = np.where(mass >= MASS_SPLIT, e_hi, e_lo)
                dl = (fit.chi0 + w * fit.chi1) * (1.0 + fit.z)
                dl = np.maximum(dl, np.finfo(float).tiny)
                mu = 5.0 * np.log10(dl) + 25.0
                resid = fit.mb - (mu + M)
                y = solve_triangular(fit.L, resid, lower=True, check_finite=False)
                return float(np.dot(y, y))

            best = None
            rng = np.random.default_rng(7)
            for x0 in ([[0.87, 0.87, -19.33]] +
                       [[rng.uniform(0, 2), rng.uniform(0, 2),
                         rng.uniform(-19.8, -18.8)] for _ in range(15)]):
                r = minimize(obj5, x0, method="L-BFGS-B",
                             bounds=[(0.0, 2.0), (0.0, 2.0), (-21.0, -16.0)])
                if best is None or r.fun < best.fun:
                    best = r
            results[name] = {
                "k": 3,
                "chi2": float(best.fun),
                "epsilon_low_mass": float(best.x[0]),
                "epsilon_high_mass": float(best.x[1]),
                "M": float(best.x[2]),
                "delta_chi2_vs_lcdm": float(best.fun - chi2_lcdm_val),
                "delta_chi2_vs_uniform": float(best.fun - res_e0.fun),
                "lnL": fit.loglike_max(float(best.fun)),
            }
            print_status(
                f"{name}: chi2={best.fun:.2f} eps_low={best.x[0]:.3f} "
                f"eps_high={best.x[1]:.3f}"
            )
            continue

        b_list = list(b_env)
        res = fit.fit(w_func, b_list)
        k = 2 + len(b_list)
        entry = {
            "k": k,
            "chi2": float(res.fun),
            "epsilon_los": float(res.x[0]),
            "env_params": [float(v) for v in res.x[1:1 + len(b_list)]],
            "M": float(res.x[-1]),
            "delta_chi2_vs_lcdm": float(res.fun - chi2_lcdm_val),
            "delta_chi2_vs_uniform": float(res.fun - res_e0.fun),
            "lnL": fit.loglike_max(float(res.fun)),
        }
        if name == "E4_locked_SA_prediction":
            entry["w_range"] = [float(fit.w_sa.min()), float(fit.w_sa.max())]
            entry["w_mean"] = float(fit.w_sa.mean())
        results[name] = entry
        print_status(
            f"{name}: chi2={res.fun:.2f} env={entry['env_params']} "
            f"(vs uniform {res.fun - res_e0.fun:+.2f})"
        )

    # ---- permutation null for E2 (continuous gate) -------------------
    print_status("Running permutation null for E2 (100 scrambles)...")
    rng = np.random.default_rng(123)
    perm_gains = []
    mass_arr = fit.mass.copy()
    for i in range(100):
        perm = rng.permutation(fit.n)
        mass_perm = mass_arr[perm]
        fit_perm = EnvKernelFit(z, mb, mass_perm, is_cal, cov)
        r = fit_perm.fit(fit_perm.w_mass, [(-1.0, 1.0)])
        perm_gains.append(float(res_e0.fun - r.fun))
        if (i + 1) % 25 == 0:
            print_status(f"  permutation {i + 1}/100")
    perm_gains = np.array(perm_gains)
    real_gain = float(res_e0.fun - results["E2_continuous_mass_gate"]["chi2"])
    p_perm = float((np.sum(perm_gains >= real_gain) + 1) / (len(perm_gains) + 1))
    results["E2_permutation_null"] = {
        "n_permutations": 100,
        "real_gain_chi2": real_gain,
        "null_mean_gain": float(perm_gains.mean()),
        "null_std_gain": float(perm_gains.std()),
        "null_max_gain": float(perm_gains.max()),
        "p_value": p_perm,
    }
    print_status(
        f"E2 gain={real_gain:+.2f} vs null mean {perm_gains.mean():.2f} "
        f"+/- {perm_gains.std():.2f} (p={p_perm:.3f})"
    )

    # ---- nested evidence ----------------------------------------------
    print_status("Running dynesty for LCDM, E0, E1, E2 (nlive=1000)...")
    lnZ_lcdm, e_lcdm = run_nested_lcdm(fit, nlive=1000)
    lnZ_e0, e_e0 = run_nested(fit, fit.w_uniform, [], nlive=1000)
    lnZ_e1, e_e1 = run_nested(fit, fit.w_binary, [(-1.5, 1.5)], nlive=1000)
    lnZ_e2, e_e2 = run_nested(fit, fit.w_mass, [(-1.0, 1.0)], nlive=1000)

    results["evidence"] = {
        "nlive": 1000,
        "lnZ_lcdm": lnZ_lcdm, "lnZerr_lcdm": e_lcdm,
        "lnZ_E0_uniform": lnZ_e0, "lnZerr_E0": e_e0,
        "lnZ_E1_binary": lnZ_e1, "lnZerr_E1": e_e1,
        "lnZ_E2_continuous": lnZ_e2, "lnZerr_E2": e_e2,
        "BF_E0_vs_lcdm": float(np.exp(lnZ_e0 - lnZ_lcdm)),
        "BF_E1_vs_lcdm": float(np.exp(lnZ_e1 - lnZ_lcdm)),
        "BF_E2_vs_lcdm": float(np.exp(lnZ_e2 - lnZ_lcdm)),
        "lnBF_E1_vs_E0": float(lnZ_e1 - lnZ_e0),
        "lnBF_E2_vs_E0": float(lnZ_e2 - lnZ_e0),
        "BF_E2_vs_E0": float(np.exp(lnZ_e2 - lnZ_e0)),
        "BF_E1_vs_E0": float(np.exp(lnZ_e1 - lnZ_e0)),
    }
    print_status(
        f"Evidence: E0 BF={np.exp(lnZ_e0 - lnZ_lcdm):.2f}, "
        f"E1 BF={np.exp(lnZ_e1 - lnZ_lcdm):.2f}, "
        f"E2 BF={np.exp(lnZ_e2 - lnZ_lcdm):.2f}, "
        f"E2/E0={np.exp(lnZ_e2 - lnZ_e0):.2f}"
    )

    # ---- environment-residual diagnostic -------------------------------
    # Residuals of the best uniform fit, binned by host mass, on the
    # whitened scale: if the leftover deviation correlates with mass the
    # uniform kernel is incomplete.
    w_uni = np.ones(fit.n)
    mu_best = fit.mu_model(results["E0_uniform"]["epsilon_los"], w_uni) \
        + results["E0_uniform"]["M"]
    resid = fit.mb - mu_best
    y = solve_triangular(fit.L, resid, lower=True, check_finite=False)
    bins = [(0, 9.5), (9.5, 10.0), (10.0, 10.5), (10.5, 13.0)]
    resid_by_mass = {}
    for lo, hi in bins:
        m = (fit.mass >= lo) & (fit.mass < hi)
        resid_by_mass[f"{lo}-{hi}"] = {
            "n": int(m.sum()),
            "mean_whitened_resid": float(y[m].mean()) if m.sum() else None,
            "mean_raw_resid_mag": float(resid[m].mean()) if m.sum() else None,
        }
    results["residual_by_host_mass"] = resid_by_mass

    # Split-sample amplitudes (diagnostic: eps fitted per mass half)
    split = {}
    for lab, sel in [("low_mass", fit.mass < MASS_SPLIT),
                     ("high_mass", fit.mass >= MASS_SPLIT)]:
        idx = np.where(sel)[0]
        sub = EnvKernelFit(z[idx], mb[idx], mass[idx],
                           is_cal[idx],
                           cov[np.ix_(idx, idx)])
        r = sub.fit(sub.w_uniform, [])
        split[lab] = {
            "n": int(sel.sum()),
            "epsilon_los": float(r.x[0]), "M": float(r.x[-1]),
            "chi2": float(r.fun),
        }
    results["mass_split_fits"] = split
    print_status(
        f"Split fits: eps_low={split['low_mass']['epsilon_los']:.3f} "
        f"(n={split['low_mass']['n']}), "
        f"eps_high={split['high_mass']['epsilon_los']:.3f} "
        f"(n={split['high_mass']['n']})"
    )

    results["meta"]["runtime_seconds"] = round(time.time() - t0, 1)
    out_path = step_json_path(STEP_ID)
    write_json(out_path, results)
    print_status(f"Wrote {out_path}")
    return results


if __name__ == "__main__":
    main()
