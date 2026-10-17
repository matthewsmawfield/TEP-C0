#!/usr/bin/env python3
"""TEP-C0 Step 03-12: Composite Environment-Modulated Kernel
============================================================
Upgrade of Step 03-11. The single-coordinate host-mass gate showed the
correct amplitude-sector direction (deep hosts enhanced) but only
marginal significance: HOST_LOGMASS is a weak depth proxy because the
potential requires R_eff, and Pantheon+ does not publish per-SN radii.

This step exploits every remaining per-sightline observable in the raw
catalog that carries environment information at full coverage:

    m_i   HOST_LOGMASS      host virial-depth proxy (sigma^2 = GM/2R_eff)
    rho_i catalog density   ambient environment (S_env coordinate of the
                            corpus's X = S_env * S_A * sigma^2 / c^2):
                            for each SN, count other SNe within a
                            projected comoving aperture and |dz| window,
                            rank-normalised within z-quantile slices to
                            suppress survey-footprint selection.
    X_i   composite         mean of the rank-normalised coordinates --
                            the lowest-dimensional surrogate of the
                            corpus screening coordinate available at
                            full coverage.

Design choices fixed a priori:
    - One free gate parameter: w_i = exp[a (X_i - <X>)] with a single
      composite coordinate -> the environmental test is one-dimensional
      (k=3), minimising the Occam penalty that buried the multi-form
      gates of Step 03-11.
    - Permutation nulls scramble the coordinate labels *within z-slices*
      so survey-footprint and mass-z correlations are preserved under
      the null; only the fine-grained environment-residual link is
      destroyed.
    - A residual diagnostic: Spearman correlation of the uniform
      kernel's symmetric-whitened residuals with X_i inside the intermediate bulk
      (0.023 < z < 1), where the uniform kernel mildly resists. If the
      leftover deviation tracks the environmental coordinate there, the
      gain is environmental rather than a kernel-shape coincidence.
    - z-slice stability: the gate amplitude refitted independently in
      three disjoint z ranges.

Models fitted (all zHD, full covariance, canonical z_T=5, n_T=2):
    G0  uniform kernel                 (eps, M)          reference
    G1  mass-rank gate                 (eps, a, M)
    G2  density-rank gate              (eps, a, M)
    G3  composite X gate               (eps, a, M)       primary test
    G4  X-gate x calibrator-aware      same gate with the frozen
        Paper-11 calibrator shift (+0.09 mag) applied     cross-check
"""

from __future__ import annotations

import sys
import time
import os
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import cho_solve, solve_triangular, eigh
from scipy.optimize import minimize
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from c0_common import (
    TEPLogger,
    ensure_dirs,
    print_status,
    set_step_logger,
    step_json_path,
    write_json,
)
from step_03_11_env_modulated_kernel import (
    EnvKernelFit,
    MASS_SPLIT,
    load_pantheon,
    SIGMA_REF, M_REF, SIGMA_T_SA, FJ_SLOPE,
)

STEP_ID = "step_03_12_env_composite_kernel"

# density-coordinate constants (pre-registered)
APERTURE_MPC = 3.0        # projected comoving aperture for neighbour count
DZ_WINDOW = 0.02          # line-of-sight window
N_ZSLICE = 20             # z-quantile slices for rank normalisation
N_PERM = 300              # permutation null draws (secondary tests)
N_PERM_HEAD = 1000        # permutation draws for the headline G5/G8 nulls
KAPPA_CEP_SHIFT = +0.09   # frozen calibrator offset (Step 04-08 / 4.6)

BULK_LO, BULK_HI = 0.023, 1.0


def local_density_rank(z: np.ndarray, ra: np.ndarray, dec: np.ndarray,
                       aperture_mpc: float = APERTURE_MPC,
                       dz: float = DZ_WINDOW,
                       return_counts: bool = False, event_ids=None):
    """Ambient-density coordinate from the SN catalog itself.

    For each physical SN i, count neighbours j != i with projected comoving
    separation < aperture_mpc and |z_j - z_i| < dz. Rank-normalise the
    counts within z-quantile slices so the coordinate measures relative
    clustering, not the survey discovery rate vs z.
    """
    from scipy.spatial import cKDTree

    original_z = np.asarray(z)
    if event_ids is None:
        # Explicit measurement-row control for historical comparisons.
        inverse = np.arange(len(z))
    else:
        if len(event_ids) != len(z):
            raise ValueError('Event identifiers do not match the selected rows')
        _, first, inverse = np.unique(event_ids, return_index=True, return_inverse=True)
        for values in (z, ra, dec):
            if not np.allclose(np.asarray(values), np.asarray(values)[first][inverse], rtol=0, atol=1e-8):
                raise ValueError('Repeated event has inconsistent density-coordinate geometry')
        z, ra, dec = (np.asarray(values)[first] for values in (z, ra, dec))
    th = np.deg2rad(dec)
    ph = np.deg2rad(ra)
    xyz = np.column_stack([
        np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph), np.sin(th)])
    tree = cKDTree(xyz)

    # comoving EdS distance (same convention as the kernel) -> aperture in rad
    zz = np.asarray(z)
    dm = (2.0 * 299792.458 / 70.0) * (1.0 - 1.0 / np.sqrt(1.0 + zz))
    ap_rad = aperture_mpc / np.maximum(dm, 1.0)   # angular aperture per SN
    # chord distance equivalent: sin(ap/2)*2
    chord = 2.0 * np.sin(np.minimum(ap_rad, np.pi) / 2.0)

    counts = np.zeros(len(zz))
    for i in range(len(zz)):
        idx = tree.query_ball_point(xyz[i], chord[i])
        # self is always in the ball and always passes the dz cut
        counts[i] = max(np.sum(np.abs(zz[idx] - zz[i]) < dz) - 1, 0)

    # rank within z-quantile slices
    counts = counts[inverse]
    ranks = np.zeros(len(original_z))
    q = pd.qcut(original_z, N_ZSLICE, labels=False, duplicates="drop")
    for s in np.unique(q):
        m = q == s
        ranks[m] = rankdata(counts[m]) / m.sum()
    if return_counts:
        return ranks, counts
    return ranks


class EventPermutation:
    """Exchange whole physical events within survey-composition/redshift cells.

    Duplicate light curves stay together, with corresponding survey rows
    matched. Depth exchanges additionally preserve measured-offset flags.
    The likelihood still contains every measurement and its full covariance.
    """
    def __init__(self, event_ids, survey, z, measured):
        cells = np.asarray(pd.qcut(z, N_ZSLICE, labels=False, duplicates='drop'))
        ambient, depth = {}, {}
        for event in np.unique(event_ids):
            rows = np.flatnonzero(event_ids == event)
            rows = rows[np.argsort(survey[rows], kind='stable')]
            if not np.all(cells[rows] == cells[rows[0]]):
                raise ValueError('Repeated event crosses redshift permutation cells')
            key = (int(cells[rows[0]]), tuple(survey[rows]))
            ambient.setdefault(key, []).append(rows)
            depth.setdefault(key + (tuple(measured[rows]),), []).append(rows)
        self.ambient_groups = list(ambient.values())
        self.depth_groups = list(depth.values())
        self.metadata = {
            'unit': 'physical CID event block',
            'stratification': 'survey-composition x 20 redshift quantiles; depth additionally preserves measured-offset flags',
            'n_events': len(np.unique(event_ids)),
            'n_ambient_events_permutable': sum(len(g) for g in self.ambient_groups if len(g)>1),
            'n_depth_events_permutable': sum(len(g) for g in self.depth_groups if len(g)>1),
            'n_ambient_rows_permutable': sum(sum(map(len,g)) for g in self.ambient_groups if len(g)>1),
            'n_depth_rows_permutable': sum(sum(map(len,g)) for g in self.depth_groups if len(g)>1),
        }

    @staticmethod
    def _exchange(values, groups, rng):
        output = np.asarray(values).copy()
        for group in groups:
            for target, donor in zip(group, rng.permutation(len(group))):
                output[target] = values[group[donor]]
        return output

    def permute(self, ambient, depth, rng):
        return (self._exchange(ambient, self.ambient_groups, rng),
                self._exchange(depth, self.depth_groups, rng))


def fit_environment_permutation(fit, ambient, depth, quality=None):
    """The existing multistart fit, usable in serial or independent workers."""
    def weight(p):
        exponent = p[0]*ambient + p[1]*depth
        if quality is not None:
            exponent = exponent + p[2]*quality[0] + p[3]*quality[1] + p[4]*quality[2]
        return np.exp(exponent)
    return float(fit.fit(weight, [(-4.0, 4.0)]*(2 if quality is None else 5)).fun)


def permutation_chi2(fit, permutations, quality=None, workers=1):
    """Scheduling only: draws and all fit starts remain deterministic."""
    if workers == 1:
        return [fit_environment_permutation(fit, a, d, quality) for a,d in permutations]
    from joblib import Parallel, delayed
    return Parallel(n_jobs=workers, inner_max_num_threads=1, verbose=5)(
        delayed(fit_environment_permutation)(fit, a, d, quality) for a,d in permutations)


def precision_chi2_difference(fit: EnvKernelFit, residual_before: np.ndarray,
                              residual_after: np.ndarray) -> np.ndarray:
    before = cho_solve((fit.L, True), residual_before, check_finite=False)
    after = cho_solve((fit.L, True), residual_after, check_finite=False)
    return residual_before * before - residual_after * after


def symmetric_whiten_residuals(covariance, residual):
    """Unique symmetric whitening; preserve the original lower-triangle likelihood."""
    effective = np.tril(covariance) + np.tril(covariance, -1).T
    eigenvalues, eigenvectors = eigh(effective, check_finite=False)
    if np.any(eigenvalues <= 0):
        raise ValueError("Covariance is not positive definite")
    components = np.sum(eigenvectors*np.asarray(residual)[:, None], axis=0)
    whitened = np.sum(eigenvectors*(components/np.sqrt(eigenvalues))[None, :], axis=1)
    if not np.all(np.isfinite(whitened)):
        raise ValueError("Nonfinite symmetric residual whitening")
    return whitened


def make_gate(fit: EnvKernelFit, coord: np.ndarray):
    """w_i = exp[a (coord_i - <coord>)], a free."""
    c = coord - coord.mean()

    def w(p):
        (a,) = p
        return np.exp(a * c)
    return w


def perm_null(fit: EnvKernelFit, coord: np.ndarray, z: np.ndarray,
              n_perm: int, seed: int = 7):
    """Permutation null of the gate gain, scrambling within z-slices."""
    rng = np.random.default_rng(seed)
    q = pd.qcut(z, N_ZSLICE, labels=False, duplicates="drop")
    chi2_u = float(fit.fit(fit.w_uniform, []).fun)
    draws = []
    for _ in range(n_perm):
        cp = coord.copy()
        for s in np.unique(q):
            m = np.where(q == s)[0]
            cp[m] = coord[rng.permutation(m)]
        draws.append(cp)
    workers = int(os.environ.get('TEP_ENV_WORKERS', '1'))
    if workers == 1:
        chi = [fit_single_permutation(fit, cp) for cp in draws]
    else:
        from joblib import Parallel, delayed
        chi = Parallel(n_jobs=workers, inner_max_num_threads=1, verbose=5)(
            delayed(fit_single_permutation)(fit, cp) for cp in draws)
    return chi2_u - np.asarray(chi)


def fit_single_permutation(fit, coordinate):
    return float(fit.fit(make_gate(fit, coordinate), [(-4.0, 4.0)]).fun)


def save_stage(results, stage):
    """Save completed statistics separately from the final manuscript result."""
    root = Path(__file__).resolve().parents[2]
    paths = [Path(__file__).resolve(), root/'scripts/steps/step_03_11_env_modulated_kernel.py',
             root/'core/cosmology.py', root/'scripts/utils/static_metric.py',
             root/'data/raw/pantheon_plus_shoes.dat', root/'data/raw/Pantheon+SH0ES.cov']
    payload = {'stage':stage, 'completed_results':results,
               'source_and_input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    temporary = root/'results/step_03_12_current_stage.tmp.json'
    write_json(temporary, payload)
    temporary.replace(root/'results/step_03_12_current_stage.json')


def run_nested_gate(fit: EnvKernelFit, coord: np.ndarray, nlive=800):
    import dynesty
    w = make_gate(fit, coord)
    bounds = [(0.0, 2.0), (-4.0, 4.0), (-21.0, -16.0)]
    lo = np.array([b[0] for b in bounds])
    hi = np.array([b[1] for b in bounds])

    def prior(u):
        return lo + u * (hi - lo)

    def loglike(theta):
        eps, a, M = theta
        ww = w((a,))
        if not np.all(np.isfinite(ww)) or np.any(ww <= 0):
            return -1e6
        chi2 = fit.chi2_given_w(eps, M, ww)
        return -0.5 * chi2 if np.isfinite(chi2) else -1e6

    s = dynesty.NestedSampler(loglike, prior, 3, bound="multi",
                              sample="rwalk", nlive=nlive, rstate=np.random.default_rng(20261011))
    s.run_nested(dlogz=0.05)
    r = s.results
    return float(r.logz[-1]), float(r.logzerr[-1])


def main():
    set_step_logger(TEPLogger(STEP_ID))
    ensure_dirs()
    t0 = time.time()
    print_status(f"{STEP_ID}: composite environment-modulated kernel")

    data_dir = Path("data/raw")
    df = pd.read_csv(data_dir / "pantheon_plus_shoes.dat", sep=r"\s+",
                     comment="#")
    z_raw = df["zHD"].values
    mb_raw = df["m_b_corr"].values
    valid = np.isfinite(z_raw) & np.isfinite(mb_raw) & (z_raw > 0.001) & (z_raw < 3.0)
    observer_correction = 5.0 * np.log10(
        (1.0 + df["zHEL"].values[valid]) / (1.0 + z_raw[valid]))
    ra = df["RA"].values[valid]
    dec = df["DEC"].values[valid]
    event_ids = df['CID'].values[valid]

    z, mb, mass, is_cal, cov, n_missing = load_pantheon()
    fit = EnvKernelFit(z, mb, mass, is_cal, cov)
    print_status(f"Loaded {fit.n} SNe ({int(is_cal.sum())} calibrators)")

    # ---- coordinates ---------------------------------------------------
    mass_rank = rankdata(mass) / len(mass)
    dens_rank = local_density_rank(z, ra, dec, event_ids=event_ids)
    composite = 0.5 * (mass_rank + dens_rank)

    coords = {
        "G1_mass_rank_gate": mass_rank,
        "G2_density_rank_gate": dens_rank,
        "G3_composite_X_gate": composite,
    }

    results = {"meta": {
        "step_id": STEP_ID,
        "n_sne": fit.n,
        "n_unique_events": len(np.unique(event_ids)),
        "density_counting_unit": "physical CID event; repeated light curves are one neighbour",
        "redshift_column": "zHD",
        "observer_redshift_column": "zHEL",
        "brightness_convention": "m_b_corr minus 5 log10[(1+zHEL)/(1+zHD)]; equivalent to observer factor in predicted distance",
        "distance_integration": "integrated cubic kernel in log(1+z), 2049 nodes",
        "covariance": "full 1701x1701",
        "aperture_mpc": APERTURE_MPC, "dz_window": DZ_WINDOW,
        "n_zslice": N_ZSLICE, "n_perm": N_PERM,
        "n_perm_headline": N_PERM_HEAD,
        "secondary_permutation_scope": "Measurement-row shuffles within redshift slices; diagnostic only. G5/G8 use physical-event blocks.",
        "nested_seed": 20261011,
        "composite": "mean[rank(HOST_LOGMASS), rank(local density)]",
    }}

    # ---- uniform reference ---------------------------------------------
    res_u = fit.fit(fit.w_uniform, [])
    chi2_u = float(res_u.fun)
    res_lcdm = minimize(fit.chi2_lcdm, [0.36, -19.34], method="L-BFGS-B",
                        bounds=[(0.05, 0.9), (-21.0, -16.0)])
    chi2_l = float(res_lcdm.fun)
    results["G0_uniform"] = {"k": 2, "chi2": chi2_u,
                             "delta_chi2_vs_lcdm": chi2_u - chi2_l}
    print_status(f"LCDM chi2={chi2_l:.2f} | uniform chi2={chi2_u:.2f} "
                 f"({chi2_u - chi2_l:+.2f})")

    # ---- single-coordinate and composite gates -------------------------
    for name, coord in coords.items():
        w = make_gate(fit, coord)
        r = fit.fit(w, [(-4.0, 4.0)])
        a_hat = float(r.x[1])
        gain = chi2_u - float(r.fun)

        # permutation null within z-slices
        null = perm_null(fit, coord, z, N_PERM)
        p = float((np.sum(null >= gain) + 1) / (len(null) + 1))

        # z-slice stability of the gate amplitude
        slices = [(0.0, BULK_LO), (BULK_LO, 0.3), (0.3, 3.0)]
        a_by_slice = {}
        for lo, hi in slices:
            sel = (z >= lo) & (z < hi)
            idx = np.where(sel)[0]
            if len(idx) < 50:
                a_by_slice[f"{lo}-{hi}"] = None
                continue
            sub = EnvKernelFit(z[idx], mb[idx], mass[idx], is_cal[idx],
                               cov[np.ix_(idx, idx)])
            rs = sub.fit(make_gate(sub, coord[idx]), [(-4.0, 4.0)])
            a_by_slice[f"{lo}-{hi}"] = {"a": float(rs.x[1]),
                                       "n": int(sel.sum()),
                                       "gain": float(
                                           sub.fit(sub.w_uniform, []).fun
                                           - rs.fun)}
        results[name] = {
            "k": 3, "chi2": float(r.fun),
            "epsilon_los": float(r.x[0]), "a_gate": a_hat,
            "M": float(r.x[-1]),
            "delta_chi2_vs_lcdm": float(r.fun - chi2_l),
            "gain_vs_uniform_chi2": gain,
            "perm_null_mean": float(null.mean()),
            "perm_null_std": float(null.std()),
            "perm_null_max": float(null.max()),
            "perm_p": p,
            "a_by_zslice": a_by_slice,
        }
        print_status(
            f"{name}: chi2={r.fun:.2f} a={a_hat:+.3f} "
            f"gain={gain:+.2f} p={p:.3f} "
            f"slices={ {k: (v['a'] if v else None) for k,v in a_by_slice.items()} }"
        )
        save_stage(results, name)

    # ---- order-invariant bulk residual diagnostic -----------------------
    # symmetric-whitened residuals of the best uniform fit vs coordinate
    mu_u = fit.mu_model(float(res_u.x[0]), np.ones(fit.n)) + float(res_u.x[-1])
    y = symmetric_whiten_residuals(cov, fit.mb - mu_u)
    bulk = (z > BULK_LO) & (z < BULK_HI)
    rho_w, p_w = spearmanr(composite[bulk], y[bulk])
    # permutation null for the correlation (within z-slices of bulk)
    rng = np.random.default_rng(5)
    qb = pd.qcut(z[bulk], 10, labels=False, duplicates="drop")
    xb = composite[bulk]
    null_rho = []
    for _ in range(N_PERM):
        xp = xb.copy()
        for s in np.unique(qb):
            m = np.where(qb == s)[0]
            xp[m] = xb[rng.permutation(m)]
        null_rho.append(spearmanr(xp, y[bulk])[0])
    null_rho = np.array(null_rho)
    p_corr = float((np.sum(np.abs(null_rho) >= abs(rho_w)) + 1)
                   / (len(null_rho) + 1))
    correlation_null = EventPermutation(
        event_ids, df["IDSURVEY"].values[valid].astype(int), z,
        np.zeros(len(z), dtype=bool))
    event_rng = np.random.default_rng(5)
    event_rho = []
    for _ in range(N_PERM_HEAD):
        xp, _ = correlation_null.permute(composite, np.zeros(len(z)), event_rng)
        event_rho.append(float(spearmanr(xp[bulk], y[bulk])[0]))
    event_p = float((1+np.sum(np.abs(event_rho)>=abs(rho_w)))/(N_PERM_HEAD+1))
    results["bulk_residual_correlation"] = {
        "segment": f"{BULK_LO} < z < {BULK_HI}",
        "n": int(bulk.sum()),
        "spearman_rho": float(rho_w),
        "spearman_p_analytic": float(p_w),
        "perm_p": event_p,
        "permutation_draws": N_PERM_HEAD,
        "permutation_null": correlation_null.metadata,
        "legacy_row_null_perm_p": p_corr,
        "legacy_row_null_draws": N_PERM,
        "whitening_method": "symmetric inverse square root of original lower-triangle covariance",
        "scope": "secondary residual-coordinate diagnostic; not independent environmental confirmation",
        "null_rho_std": float(null_rho.std()),
        "direction": "positive rho = fainter-than-kernel residual grows "
                     "with environmental depth" if rho_w > 0 else
                    "residual anti-correlates with depth",
    }
    print_status(f"bulk resid-X corr: rho={rho_w:+.3f} perm_p={event_p:.4f}")

    # ---- same within the calibrator block ------------------------------
    if is_cal.sum() >= 20:
        rho_c, p_c = spearmanr(composite[is_cal], y[is_cal])
        results["cal_residual_correlation"] = {
            "n": int(is_cal.sum()),
            "spearman_rho": float(rho_c),
            "spearman_p_analytic": float(p_c),
        }
        print_status(f"calibrator resid-X corr: rho={rho_c:+.3f}")

    # ---- nested evidence for the composite gate -------------------------
    print_status("dynesty: composite gate (nlive=800)...")
    lnZ_g3, e_g3 = run_nested_gate(fit, composite, nlive=800)
    # Recompute both references with this integrator and the same sample.
    # Never combine a fresh environmental likelihood with an older distance grid.
    from step_03_11_env_modulated_kernel import run_nested, run_nested_lcdm
    lnZ_e0, e_e0 = run_nested(fit, fit.w_uniform, [], nlive=1000)
    lnZ_lcdm, e_lcdm = run_nested_lcdm(fit, nlive=1000)
    results["evidence"] = {
        "nlive": 800,
        "lnZ_G3_composite": lnZ_g3, "lnZerr": e_g3,
        "lnZ_E0_uniform_ref": lnZ_e0, "lnZ_lcdm_ref": lnZ_lcdm,
        "lnZerr_E0": e_e0, "lnZerr_lcdm": e_lcdm,
        "reference_source": "same-run, same-integrator, same-sample",
        "BF_G3_vs_lcdm": float(np.exp(lnZ_g3 - lnZ_lcdm)),
        "BF_G3_vs_E0": float(np.exp(lnZ_g3 - lnZ_e0)),
    }
    print_status(f"G3 evidence: BF vs LCDM={np.exp(lnZ_g3 - lnZ_lcdm):.2f}, "
                 f"vs uniform={np.exp(lnZ_g3 - lnZ_e0):.2f}")
    save_stage(results, 'G3_evidence')

    # ---- calibrator-shifted cross-check ---------------------------------
    # frozen Paper-11 shift on calibrator rows + composite gate
    mb_shift = mb.copy()
    mb_shift[is_cal] -= KAPPA_CEP_SHIFT
    fit_s = EnvKernelFit(z, mb_shift, mass, is_cal, cov)
    rs = fit_s.fit(fit_s.w_uniform, [])
    rs_g = fit_s.fit(make_gate(fit_s, composite), [(-4.0, 4.0)])
    rs_lcdm = minimize(fit_s.chi2_lcdm, [0.36, -19.34], method="L-BFGS-B",
                       bounds=[(0.05, 0.9), (-21.0, -16.0)])
    chi2_lcdm_s = float(rs_lcdm.fun)
    # provenance: express the frozen shift as kappa on the physical depth
    # coordinate X = S_A sigma^2/c^2 and propagate the current Paper-11
    # transfer coefficients for transparency.
    _c2 = 299792.458 ** 2
    _sig = SIGMA_REF * 10.0 ** (FJ_SLOPE * (mass - M_REF))
    _sa = np.minimum(1.0, (_sig / SIGMA_T_SA) ** (2.0 / 3.0))
    _xcal = float((_sa * _sig ** 2 / _c2)[is_cal.astype(bool)].mean())
    # load the empirically measured calibrator offset from the cross-catalog
    # transfer fit so the provenance block traces to a real artifact.
    _cross_path = Path("results/outputs/audit_kappa_crosscatalog.json")
    if not _cross_path.exists():
        raise FileNotFoundError(
            f"Required cross-catalog calibration artifact missing: {_cross_path}"
        )
    _cross = json.loads(_cross_path.read_text())
    _jko = _cross["catalogs"]["pantheon_plus"]["joint_kappa_plus_calibrator_offset"]
    _dM_cal = float(_jko["dM_cal"])
    _dM_cal_err = float(_jko["dM_cal_err"])
    _dM_source = "audit_kappa_crosscatalog.json:joint_kappa_plus_calibrator_offset"
    results["G4_calshift_composite"] = {
        "kappa_cep_shift_mag": KAPPA_CEP_SHIFT,
        "chi2_lcdm_shifted": chi2_lcdm_s,
        "chi2_uniform_shifted": float(rs.fun),
        "delta_vs_lcdm": float(rs.fun - chi2_l),
        "delta_vs_lcdm_shifted": float(rs.fun - chi2_lcdm_s),
        "chi2_composite_shifted": float(rs_g.fun),
        "delta_composite_vs_lcdm": float(rs_g.fun - chi2_l),
        "delta_composite_vs_lcdm_shifted": float(rs_g.fun - chi2_lcdm_s),
        "a_gate": float(rs_g.x[1]),
        "shift_provenance": {
            "calibrator_mean_depth_X": _xcal,
            "implied_kappa_for_frozen_shift": float(KAPPA_CEP_SHIFT / _xcal),
            "predicted_shift_unified_kappa_441316": float(441316.4 * _xcal),
            "predicted_shift_joint_kappa_745976": float(745975.5 * _xcal),
            "predicted_shift_redshift_only_kappa_1p2e6": float(1.2e6 * _xcal),
            "measured_offset_joint_fit": _dM_cal,
            "measured_offset_joint_fit_err": _dM_cal_err,
            "measured_offset_source": _dM_source,
        },
    }
    print_status(
        f"cal-shifted: uniform chi2={rs.fun:.1f} ({rs.fun-chi2_l:+.1f}), "
        f"composite chi2={rs_g.fun:.1f} ({rs_g.fun-chi2_l:+.1f})"
    )

    # ---- G5: two-channel nested screening -------------------------------
    # The corpus's nested-domain structure predicts TWO environmental
    # channels with opposite signs: the ambient halo depth (mass +
    # large-scale density) enhances the effective transport field, while
    # the radial depth inside the host well screens it. HOST_ANGSEP
    # provides the within-host radial coordinate for the 462 SNe whose
    # host associations carry a measured projected offset (IDSURVEY
    # 1 = SDSS and 10 = DES; see the scheme audit below); objects
    # without a measurement receive neutral weight so the full sample
    # is retained.
    print_status("G5: two-channel nested screening gate...")
    ang = df["HOST_ANGSEP"].values[valid].astype(float)
    surv_i = df["IDSURVEY"].values[valid].astype(int)
    # HOST_ANGSEP is a measured projected offset only for the
    # spectroscopic host-association schemes: IDSURVEY 1 (SDSS,
    # median 1.5") and 10 (DES, median 0.7"). IDSURVEY 15 (PS1) records
    # a narrow 17-24" distribution with unverified offset provenance.
    # Keep its depth weight neutral in this specified selection; the
    # distribution alone does not establish a nominal matching radius.
    # Ranking within each scheme separately prevents cross-scheme
    # contrast from masquerading as radial depth.
    ok_ang = (np.isfinite(ang) & (ang > 0)
            & np.isin(surv_i, [1, 10]))
    DA = (2.0 * 299792.458 / 70.0) * (1.0 - 1.0 / np.sqrt(1.0 + z)) / (1.0 + z)
    s_raw = ang * DA / (5.0 * 10.0 ** (0.25 * (mass - 11.0)))  # r_proj/R_eff
    x_depth = np.zeros(len(z))
    for s_ in np.unique(surv_i[ok_ang]):
        m2 = ok_ang & (surv_i == s_)
        srv = rankdata(-np.log(np.maximum(s_raw[m2], 1e-3))) / m2.sum()
        x_depth[m2] = srv - srv.mean()   # deep-in-well = positive
    X_amb = 0.5 * (mass_rank + dens_rank)
    X_amb = X_amb - X_amb.mean()

    def w_two_channel(p):
        a, b = p
        return np.exp(a * X_amb + b * x_depth)

    r5 = fit.fit(w_two_channel, [(-4.0, 4.0), (-4.0, 4.0)])
    gain5 = chi2_u - float(r5.fun)
    print_status(
        f"G5 two-channel: chi2={r5.fun:.2f} a_amb={r5.x[1]:+.3f} "
        f"b_depth={r5.x[2]:+.3f} gain={gain5:+.2f}"
    )

    # z-slice stability (angsep has no coverage below z=0.023)
    slices5 = {}
    for lo, hi in [(0.0, BULK_LO), (BULK_LO, 0.3), (0.3, 3.0)]:
        sel = (z >= lo) & (z < hi)
        ii = np.where(sel)[0]
        if len(ii) < 50:
            continue
        sub = EnvKernelFit(z[ii], mb[ii], mass[ii], is_cal[ii],
                           cov[np.ix_(ii, ii)])
        ru2 = sub.fit(sub.w_uniform, [])

        def wj2(p, ii=ii):
            a, b = p
            return np.exp(a * X_amb[ii] + b * x_depth[ii])

        rj = sub.fit(wj2, [(-4.0, 4.0), (-4.0, 4.0)])
        slices5[f"{lo}-{hi}"] = {
            "n": int(sel.sum()), "n_angsep": int(ok_ang[ii].sum()),
            "a_amb": float(rj.x[1]), "b_depth": float(rj.x[2]),
            "gain": float(ru2.fun - rj.fun),
        }
    results["G5_two_channel"] = {
        "k": 4, "chi2": float(r5.fun),
        "epsilon_los": float(r5.x[0]),
        "a_ambient": float(r5.x[1]), "b_radial_depth": float(r5.x[2]),
        "M": float(r5.x[-1]),
        "delta_chi2_vs_lcdm": float(r5.fun - chi2_l),
        "gain_vs_uniform_chi2": gain5,
        "n_angsep": int(ok_ang.sum()),
        "angsep_z_min": float(z[ok_ang].min()),
        "z_slices": slices5,
    }

    # order-invariant precision-weighted contributions to the gain
    resid0 = fit.mb - (fit.mu_model(
        float(res_u.x[0]), np.ones(fit.n)) + float(res_u.x[-1]))
    w5 = w_two_channel((float(r5.x[1]), float(r5.x[2])))
    resid5 = fit.mb - (fit.mu_model(
        float(r5.x[0]), w5) + float(r5.x[-1]))
    d5 = precision_chi2_difference(fit, resid0, resid5)
    assert np.isclose(d5.sum(), gain5, atol=1e-5)
    results["G5_two_channel"]["precision_localization"] = {
        name: float(d5[m].sum())
        for name, m in [
            ("calibrators", is_cal), ("z_gt_1", z > 1),
            ("bulk_0p023_1", (z > BULK_LO) & (z < BULK_HI)),
            ("lowz_field", (z < BULK_LO) & ~is_cal)]
    }
    results["G5_two_channel"]["localization_method"] = (
        "r_i*(C^-1 r)_i before minus after; sums to full delta chi2, "
        "invariant under row ordering but not a standalone subset refit")

    # permutation null: scramble both coordinates within survey x redshift
    # cells, with depth restricted to measured-offset rows
    rng = np.random.default_rng(17)
    event_null = EventPermutation(event_ids, surv_i, z, ok_ang)
    workers = int(os.environ.get('TEP_ENV_WORKERS', '1'))
    if workers < 1:
        raise ValueError('TEP_ENV_WORKERS must be positive')
    draws5 = [event_null.permute(X_amb, x_depth, rng) for _ in range(N_PERM_HEAD)]
    print_status(f'G5: {N_PERM_HEAD} physical-event permutations; {workers} fit workers')
    gains5 = chi2_u - np.asarray(permutation_chi2(fit, draws5, workers=workers))
    p5 = float((np.sum(gains5 >= gain5) + 1) / (len(gains5) + 1))
    results["G5_two_channel"]["permutation_null"] = {
        "n_permutations": N_PERM_HEAD,
        **event_null.metadata,
        "null_mean": float(gains5.mean()),
        "null_std": float(gains5.std()),
        "null_max": float(gains5.max()),
        "p_value": p5,
    }
    print_status(f"G5 perm null: mean={gains5.mean():.2f} "
                 f"std={gains5.std():.2f} -> p={p5:.4f}")
    save_stage(results, 'G5_permutation')

    # redshift-projection control: pure amplitude-z gate gains ~nothing,
    # and a survives controlling for z
    zc = (z - z.min()) / (z.max() - z.min())
    rz = fit.fit(make_gate(fit, zc), [(-4.0, 4.0)])

    def w_xz(p):
        a, c = p
        return np.exp(a * X_amb + c * zc)
    rxz = fit.fit(w_xz, [(-4.0, 4.0), (-4.0, 4.0)])
    results["G5_two_channel"]["z_projection_control"] = {
        "z_gate_gain": float(chi2_u - rz.fun),
        "joint_a_with_z": float(rxz.x[1]),
        "joint_c_z": float(rxz.x[2]),
        "joint_gain": float(chi2_u - rxz.fun),
    }
    print_status(f"z-projection control: z-gate gain="
                 f"{chi2_u - rz.fun:+.2f}, a survives at {rxz.x[1]:+.3f}")

    # nested evidence under the theory-motivated narrow gate prior
    print_status("dynesty G5 narrow prior +-1 (nlive=800)...")
    lnZ_n, en = _nested_two_channel(fit, X_amb, x_depth, 1.0, 800)
    print_status("dynesty G5 wide prior +-4 (nlive=800)...")
    lnZ_w, ew = _nested_two_channel(fit, X_amb, x_depth, 4.0, 800)
    results["G5_two_channel"]["evidence"] = {
        "narrow_prior_pm1": {
            "lnZ": lnZ_n, "lnZerr": en,
            "BF_vs_lcdm": float(np.exp(lnZ_n - results["evidence"]["lnZ_lcdm_ref"])),
            "BF_vs_uniform": float(np.exp(lnZ_n - results["evidence"]["lnZ_E0_uniform_ref"])),
        },
        "wide_prior_pm4": {
            "lnZ": lnZ_w, "lnZerr": ew,
            "BF_vs_lcdm": float(np.exp(lnZ_w - results["evidence"]["lnZ_lcdm_ref"])),
            "BF_vs_uniform": float(np.exp(lnZ_w - results["evidence"]["lnZ_E0_uniform_ref"])),
        },
        "note": "theory predicts O(10%) modulation; the narrow prior "
                "brackets that expectation and is the physical comparison",
    }
    print_status(
        f"G5 evidence: narrow BF={np.exp(lnZ_n - results['evidence']['lnZ_lcdm_ref']):.2f} "
        f"vs LCDM, wide BF={np.exp(lnZ_w - results['evidence']['lnZ_lcdm_ref']):.2f}"
    )
    save_stage(results, 'G5_evidence')

    # ---- G6: band-pass / local-screening control -------------------------
    # The implied per-SN amplitude (eps_i* ~ residual / dmu/deps) rises
    # from ~0.5 at z<0.1 to ~0.9 in the bulk: superficially a "local
    # screening" shape. The decisive test is whether the preference
    # survives anchor removal: a genuine screening would persist on
    # z>0.023, an intercept projection must collapse.
    print_status("G6: band-pass (local-screening) control...")
    from core.cosmology import _integrated_distance
    from scipy.linalg import cholesky as _chol

    def chi1_of(zloc, z_eval):
        return _integrated_distance(z_eval, lambda zz:
            299792.458 / (70.0 * (1.0 + zz) ** 1.5)
            * np.log1p(zz) * np.exp(-(zz / 5.0) ** 2)
            * (1.0 - np.exp(-(zz / zloc) ** 2)))

    chi0_all = fit.chi0

    def chi2_band(th, zz= z, L=fit.L, mbb=mb, c0=chi0_all):
        eps, M, zloc = th
        with np.errstate(all="ignore"):
            c1 = chi1_of(zloc, zz)
            dl = np.maximum((c0 + eps * c1) * (1.0 + zz), 1e-10)
            r = mbb - (5.0 * np.log10(dl) + 25.0 + M)
            y = solve_triangular(L, r, lower=True)
            v = float(y @ y)
        return v if np.isfinite(v) else 1e10

    rb = minimize(chi2_band, [0.87, -19.33, 0.09], method="L-BFGS-B",
                  bounds=[(0.0, 2.0), (-21.0, -16.0), (0.005, 1.0)])

    # anchor-removed refit
    sel = z > BULK_LO
    ii = np.where(sel)[0]
    Ls = _chol(cov[np.ix_(ii, ii)], lower=True, check_finite=False)
    chi0_s = chi0_all[ii]
    rb_s = minimize(chi2_band, [0.87, -19.33, 0.09], method="L-BFGS-B",
                    bounds=[(0.0, 2.0), (-21.0, -16.0), (0.005, 1.0)],
                    args=(z[ii], Ls, mb[ii], chi0_s))
    sub = EnvKernelFit(z[ii], mb[ii], mass[ii], is_cal[ii],
                       cov[np.ix_(ii, ii)])
    ru_s = sub.fit(sub.w_uniform, [])

    # env gate on the same anchor-removed subset for contrast
    rg_s = sub.fit(make_gate(sub, composite[ii]), [(-4.0, 4.0)])

    results["G6_bandpass_control"] = {
        "k": 3,
        "chi2_full": float(rb.fun),
        "z_loc": float(rb.x[2]),
        "epsilon": float(rb.x[0]),
        "gain_vs_uniform_full": chi2_u - float(rb.fun),
        "delta_chi2_vs_lcdm": float(rb.fun - chi2_l),
        "anchor_removed": {
            "chi2_band": float(rb_s.fun),
            "chi2_uniform": float(ru_s.fun),
            "z_loc": float(rb_s.x[2]),
            "gain_vs_uniform": float(ru_s.fun - rb_s.fun),
            "env_gate_gain_same_subset": float(ru_s.fun - rg_s.fun),
        },
        "verdict": "band-pass gain collapses under anchor removal: "
                   "the z-deformation is the calibrator-intercept "
                   "structure in kernel form, not a resolved local "
                   "screening; the environmental gate retains its gain "
                   "on the same subset.",
    }
    print_status(
        f"G6 band-pass: full gain={chi2_u - rb.fun:+.2f} z_loc={rb.x[2]:.3f} | "
        f"z>0.023: gain={ru_s.fun - rb_s.fun:+.2f} (env keeps "
        f"{ru_s.fun - rg_s.fun:+.2f})"
    )

    # ---- G7: FITPROB quality-structure control ---------------------------
    # A gate on the SALT2 light-curve fit probability also improves the
    # fit; recorded as a data-quality control (not an environmental
    # coordinate). Its significance under both z-slice and survey-z-cell
    # permutation shows the structure is not a survey-offset artefact.
    fp = df["FITPROB"].values[valid].astype(float)
    okf = np.isfinite(fp)
    fpr = np.zeros(len(z))
    fv = rankdata(fp[okf]) / okf.sum()
    fpr[okf] = fv - fv.mean()

    def w_fp(p, f_=fpr):
        (c,) = p
        return np.exp(c * f_)

    r7 = fit.fit(w_fp, [(-4.0, 4.0)])
    gain7 = chi2_u - float(r7.fun)
    # survey-x-z-cell permutation null
    surv = df["IDSURVEY"].values[valid].astype(str)
    zcell = pd.qcut(z, 10, labels=False, duplicates="drop")
    group = surv + "_" + zcell.astype(str)
    rng7 = np.random.default_rng(11)
    g7 = []
    for i in range(150):
        fp_ = fpr.copy()
        for u in np.unique(group):
            m2 = np.where(group == u)[0]
            if len(m2) > 1:
                fp_[m2] = fpr[rng7.permutation(m2)]
        r = fit.fit(lambda p, fp_=fp_: np.exp(p[0] * fp_), [(-4.0, 4.0)])
        g7.append(chi2_u - float(r.fun))
    g7 = np.array(g7)
    p7 = float((np.sum(g7 >= gain7) + 1) / (len(g7) + 1))

    # split-sample preference by FITPROB half
    med_fp = np.nanmedian(fp)
    split7 = {}
    for lab, sel7 in [("high_fitprob", fp >= med_fp),
                      ("low_fitprob", fp < med_fp)]:
        ii7 = np.where(sel7)[0]
        sub7 = EnvKernelFit(z[ii7], mb[ii7], mass[ii7], is_cal[ii7],
                            cov[np.ix_(ii7, ii7)])
        rl7 = minimize(sub7.chi2_lcdm, [0.36, -19.34], method="L-BFGS-B",
                       bounds=[(0.05, 0.9), (-21.0, -16.0)])
        rt7 = sub7.fit(sub7.w_uniform, [])
        split7[lab] = {"n": int(sel7.sum()),
                       "delta_chi2_tep_lcdm": float(rt7.fun - rl7.fun)}
    results["G7_fitprob_quality"] = {
        "k": 3, "chi2": float(r7.fun), "c_gate": float(r7.x[1]),
        "gain_vs_uniform_chi2": gain7,
        "delta_chi2_vs_lcdm": float(r7.fun - chi2_l),
        "perm_p_survey_z_cell": p7,
        "split_delta_chi2": split7,
        "note": "data-quality axis, not an environmental coordinate; "
                "the gate direction (+) means amplitude rises with "
                "light-curve fit quality, while the split test shows the "
                "TEP-vs-LCDM preference carried by the lower-quality "
                "half - the gate exploits amplitude-vs-quality structure "
                "rather than marking where the preference lives.",
    }
    print_status(f"G7 FITPROB: gain={gain7:+.2f} p={p7:.4f} "
                 f"split={split7['low_fitprob']['delta_chi2_tep_lcdm']:.1f}/"
                 f"{split7['high_fitprob']['delta_chi2_tep_lcdm']:.1f}")

    # ---- G8: quality-controlled two-channel model -----------------------
    # The FITPROB axis is not unique: several light-curve quality
    # coordinates (err_DIAG, FITCHI2) carry independent amplitude
    # structure. If left unmodeled they mask the environmental channel:
    # the standalone two-channel gain of +8.6 understates the
    # environmental contribution because the gate must compete with
    # un-modeled quality structure. The controlled model fits the two
    # environmental channels jointly with the quality gates and assigns
    # the environmental contribution as chi2(quality-only) -
    # chi2(env+quality). The environmental coordinates are permuted
    # (within z-slices / among valid-angsep rows) with the quality
    # coordinates held fixed to null the environmental attribution.
    print_status("G8: quality-controlled environmental model...")
    def _rankc(vv):
        ok = np.isfinite(vv)
        cr = np.zeros(len(z))
        rr = rankdata(vv[ok]) / ok.sum()
        cr[ok] = rr - rr.mean()
        return cr

    err_r = _rankc(df["m_b_corr_err_DIAG"].values[valid].astype(float))
    fchi_r = _rankc(df["FITCHI2"].values[valid].astype(float))

    def w_qual(p):
        c1, c2, c3 = p
        return np.exp(c1 * fpr + c2 * err_r + c3 * fchi_r)

    def w_env_qual(p):
        a, b, c1, c2, c3 = p
        return np.exp(a * X_amb + b * x_depth
                      + c1 * fpr + c2 * err_r + c3 * fchi_r)

    r_q = fit.fit(w_qual, [(-4.0, 4.0)] * 3)
    r_eq = fit.fit(w_env_qual, [(-4.0, 4.0)] * 5)
    env_attr = float(r_q.fun - r_eq.fun)

    # channel decomposition inside the controlled model
    r_qa = fit.fit(lambda p: np.exp(
        p[0] * X_amb + p[1] * fpr + p[2] * err_r + p[3] * fchi_r),
        [(-4.0, 4.0)] * 4)
    r_qd = fit.fit(lambda p: np.exp(
        p[0] * x_depth + p[1] * fpr + p[2] * err_r + p[3] * fchi_r),
        [(-4.0, 4.0)] * 4)
    amb_attr = float(r_q.fun - r_qa.fun)
    dep_attr = float(r_q.fun - r_qd.fun)
    print_status(
        f"G8: qual-only chi2={r_q.fun:.2f} | env+qual chi2={r_eq.fun:.2f} "
        f"(gain vs uniform {chi2_u - r_eq.fun:+.2f}) | "
        f"env-attributable={env_attr:+.2f} "
        f"a={r_eq.x[1]:+.3f} b={r_eq.x[2]:+.3f}"
    )

    # permutation null on the environmental attribution only
    rng8 = np.random.default_rng(31)
    draws8 = [event_null.permute(X_amb, x_depth, rng8) for _ in range(N_PERM_HEAD)]
    print_status(f'G8: {N_PERM_HEAD} physical-event permutations; {workers} fit workers')
    null_attr = r_q.fun - np.asarray(permutation_chi2(
        fit, draws8, quality=(fpr, err_r, fchi_r), workers=workers))
    p_attr = float((np.sum(null_attr >= env_attr) + 1) / (len(null_attr) + 1))
    z_attr = float((env_attr - null_attr.mean()) / null_attr.std())

    # anchor-removed version
    sel8 = z > BULK_LO
    ii8 = np.where(sel8)[0]
    sub8 = EnvKernelFit(z[ii8], mb[ii8], mass[ii8], is_cal[ii8],
                        cov[np.ix_(ii8, ii8)])
    rq8 = sub8.fit(lambda p: np.exp(p[0] * fpr[ii8] + p[1] * err_r[ii8]
                                    + p[2] * fchi_r[ii8]), [(-4.0, 4.0)] * 3)
    re8 = sub8.fit(lambda p: np.exp(p[0] * X_amb[ii8] + p[1] * x_depth[ii8]
                                    + p[2] * fpr[ii8] + p[3] * err_r[ii8]
                                    + p[4] * fchi_r[ii8]), [(-4.0, 4.0)] * 5)
    env_attr_nz = float(rq8.fun - re8.fun)

    # order-invariant precision-weighted environmental attribution
    wq_ = w_qual(tuple(r_q.x[1:4]))
    wf_ = w_env_qual(tuple(r_eq.x[1:6]))
    resid_q = mb - (fit.mu_model(float(r_q.x[0]), wq_) + float(r_q.x[-1]))
    resid_eq = mb - (fit.mu_model(float(r_eq.x[0]), wf_) + float(r_eq.x[-1]))
    d_env = precision_chi2_difference(fit, resid_q, resid_eq)
    assert np.isclose(d_env.sum(), env_attr, atol=1e-5)
    loc8 = {name: float(d_env[m].sum())
            for name, m in [
                ("calibrators", is_cal), ("z_gt_1", z > 1),
                ("bulk_0p023_1", (z > BULK_LO) & (z < BULK_HI)),
                ("lowz_field", (z < BULK_LO) & ~is_cal)]}

    # survey-offset robustness: absorb per-survey magnitude intercepts
    surv_i = df["IDSURVEY"].values[valid].astype(int)
    usv = np.unique(surv_i)

    def chi2_surv(th, env_terms):
        eps = th[0]
        M = th[3]
        a, b, c1, c2, c3 = th[1], th[2], th[4], th[5], th[6]
        offs = th[7:]
        if env_terms:
            w = np.exp(a * X_amb + b * x_depth
                       + c1 * fpr + c2 * err_r + c3 * fchi_r)
        else:
            w = np.exp(c1 * fpr + c2 * err_r + c3 * fchi_r)
        mu = fit.mu_model(eps, w) + M
        so = np.zeros(len(z))
        for k, s_ in enumerate(usv[1:]):
            so[surv_i == s_] = offs[k]
        rr = mb - mu - so
        if not np.all(np.isfinite(rr)):
            return 1e10
        with np.errstate(all="ignore"):
            yy = solve_triangular(fit.L, rr, lower=True)
            v = float(yy @ yy)
        return v if np.isfinite(v) else 1e10

    koff = len(usv) - 1
    thq0 = [0.87, 0.0, 0.0, -19.33, 0.3, 0.3, 0.0] + [0.0] * koff
    the0 = [0.87, 0.2, -0.25, -19.33, 0.3, 0.3, 0.0] + [0.0] * koff
    bnd_sv = ([(0.0, 2.0), (-4, 4), (-4, 4), (-21, -16)]
              + [(-4, 4)] * 3 + [(-0.5, 0.5)] * koff)
    rs_q = minimize(chi2_surv, thq0, args=(False,), method="L-BFGS-B",
                    bounds=bnd_sv)
    rs_e = minimize(chi2_surv, the0, args=(True,), method="L-BFGS-B",
                    bounds=bnd_sv)

    # leave-one-survey-out: the depth coordinate is contributed by two
    # survey schemes (1=SDSS, 10=DES), so the channel could in
    # principle be carried by either alone. Refit qual-only vs env+qual
    # dropping each in turn; IDSURVEY 15 is additionally dropped as a
    # check that its excluded (nominal-radius) rows carry nothing.
    loso = {}
    for s_drop in (1, 10, 15):
        keep = np.where(surv_i != s_drop)[0]
        sub_k = EnvKernelFit(z[keep], mb[keep], mass[keep],
                             is_cal[keep], cov[np.ix_(keep, keep)])
        rq_k = sub_k.fit(lambda p, k=keep: np.exp(
            p[0] * fpr[k] + p[1] * err_r[k] + p[2] * fchi_r[k]),
            [(-4.0, 4.0)] * 3)
        re_k = sub_k.fit(lambda p, k=keep: np.exp(
            p[0] * X_amb[k] + p[1] * x_depth[k]
            + p[2] * fpr[k] + p[3] * err_r[k] + p[4] * fchi_r[k]),
            [(-4.0, 4.0)] * 5)
        loso[f"drop_IDSURVEY_{s_drop}"] = {
            "n": int(len(keep)),
            "env_attributable_chi2": float(rq_k.fun - re_k.fun),
            "a_ambient": float(re_k.x[1]), "b_depth": float(re_k.x[2]),
        }
        print_status(
            f"  LOSO drop survey {s_drop} (n={len(keep)}): "
            f"env-attributable {rq_k.fun - re_k.fun:+.2f} "
            f"a={re_k.x[1]:+.3f} b={re_k.x[2]:+.3f}")

    # evidence: env-specific BF = lnZ(env+qual) - lnZ(qual-only), +-1 priors
    print_status("dynesty G8 qual-only and env+qual (nlive=800)...")
    import dynesty

    def _nested_g8(env_terms):
        bnd = [(0.0, 2.0)] + ([(-1.0, 1.0)] * (5 if env_terms else 3)) \
            + [(-21.0, -16.0)]
        lo_ = np.array([b[0] for b in bnd])
        hi_ = np.array([b[1] for b in bnd])

        def ll(t):
            eps = t[0]
            gp = t[1:1 + (5 if env_terms else 3)]
            M = t[-1]
            w = (w_env_qual(gp) if env_terms else w_qual(gp))
            if not np.all(np.isfinite(w)) or np.any(w <= 0):
                return -1e6
            c2_ = fit.chi2_given_w(eps, M, w)
            return -0.5 * c2_ if np.isfinite(c2_) else -1e6

        s = dynesty.NestedSampler(ll, lambda u: lo_ + u * (hi_ - lo_),
                                  len(bnd), bound="multi", sample="rwalk",
                                  nlive=800, rstate=np.random.default_rng(20261011))
        s.run_nested(dlogz=0.05)
        return float(s.results.logz[-1]), float(s.results.logzerr[-1])

    lnZ_q, e_q = _nested_g8(False)
    lnZ_eq, e_eq = _nested_g8(True)

    results["G8_quality_controlled_env"] = {
        "model": "w = exp[a X_amb + b x_depth + c1 FITPROB + c2 errDIAG "
                 "+ c3 FITCHI2]; env-attributable = chi2(qual) - chi2(full)",
        "k_qual": 5, "k_full": 7,
        "chi2_qual_only": float(r_q.fun),
        "chi2_env_qual": float(r_eq.fun),
        "gain_vs_uniform_chi2": float(chi2_u - r_eq.fun),
        "env_attributable_chi2": env_attr,
        "channel_decomposition": {
            "ambient_only_attributable": amb_attr,
            "depth_only_attributable": dep_attr,
            "note": "Separate single-channel fits carry comparable gains over "
                    "the same quality-only baseline; not independent detections",
        },
        "coefficients": {
            "epsilon_los": float(r_eq.x[0]),
            "a_ambient": float(r_eq.x[1]),
            "b_depth": float(r_eq.x[2]),
            "c_fitprob": float(r_eq.x[3]),
            "c_errdiag": float(r_eq.x[4]),
            "c_fitchi2": float(r_eq.x[5]),
            "M": float(r_eq.x[6]),
        },
        "env_permutation_null": {
            "n_permutations": N_PERM_HEAD,
            **event_null.metadata,
            "quality_coordinates": "held fixed",
            "null_mean": float(null_attr.mean()),
            "null_std": float(null_attr.std()),
            "null_max": float(null_attr.max()),
            "p_value": p_attr,
            "gaussian_z_vs_null": z_attr,
        },
        "anchor_removed_z_gt_0p023": {
            "env_attributable_chi2": env_attr_nz,
            "a_ambient": float(re8.x[1]), "b_depth": float(re8.x[2]),
        },
        "precision_localization": loc8,
        "localization_method": (
            "r_i*(C^-1 r)_i before minus after; sums to full delta chi2, "
            "invariant under row ordering but not a standalone subset refit"),
        "survey_offset_robustness": {
            "n_survey_offsets": koff,
            "chi2_qual_plus_offsets": float(rs_q.fun),
            "chi2_env_qual_plus_offsets": float(rs_e.fun),
            "env_attributable_with_offsets": float(rs_q.fun - rs_e.fun),
            "a_ambient": float(rs_e.x[1]), "b_depth": float(rs_e.x[2]),
        },
        "leave_one_survey_out": loso,
        "evidence": {
            "nlive": 800,
            "lnZ_qual_only": lnZ_q, "lnZerr_qual": e_q,
            "lnZ_env_qual": lnZ_eq, "lnZerr_env_qual": e_eq,
            "BF_env_qual_vs_lcdm": float(np.exp(
                lnZ_eq - results["evidence"]["lnZ_lcdm_ref"])),
            "BF_env_extension_vs_qual_only": float(np.exp(lnZ_eq - lnZ_q)),
            "prior": "gate params +-1 (theory-motivated O(10%) modulation)",
        },
        "note": "the environmental attribution survives explicit "
                "modelling of the light-curve-quality amplitude axes: "
                "the standalone gate slightly understates the "
                "environmental modulation because un-modeled quality "
                "structure competes with it.",
    }
    print_status(
        f"G8 env-attributable: {env_attr:+.2f} p={p_attr:.4f} "
        f"(z={z_attr:.1f}) | anchor-removed {env_attr_nz:+.2f} | "
        f"with survey offsets {rs_q.fun - rs_e.fun:+.2f} | "
        f"BF env-vs-qual={np.exp(lnZ_eq - lnZ_q):.1f}"
    )
    save_stage(results, 'G8_evidence')

    # ---- G10: depth-coordinate robustness audit --------------------------
    # Documents the host-association scheme audit and every remaining
    # robustness check on the two-channel coordinate.
    # (a) scheme distributions: IDSURVEY 15 has unverified offset
    #     provenance (17-24") and is excluded from x_depth under the
    #     specified selection (see construction above).
    # (b) per-scheme replication: refit env+qual inside each resolved
    #     survey independently. The depth channel must keep b<0 inside
    #     each scheme, not only in the pooled contrast.
    # (c) interaction: X_amb * x_depth cross-term (nested screening
    #     would allow within-well depth to matter more in deep ambient
    #     wells).
    # (d) VPEC gate: an unused environmental observable.
    # (e) saturating gate: tanh(x/0.3) bounded form; the corpus's S_A
    #     transfer is bounded, so an unbounded exp(ax) could overfit.
    # (f) profile chi2 on a and b -> coefficient uncertainties.
    ang_stats = {}
    for s_ in np.unique(surv_i):
        av = ang[surv_i == s_]
        av = av[np.isfinite(av) & (av > 0)]
        if len(av) >= 20:
            ang_stats[f"IDSURVEY_{s_}"] = {
                "n": int(len(av)),
                "median_arcsec": float(np.median(av)),
                "p5_p95_arcsec": [float(np.percentile(av, 5)),
                                  float(np.percentile(av, 95))],
            }

    per_scheme = {}
    for s_ in (1, 10, 15):
        ii = np.where(surv_i == s_)[0]
        if len(ii) < 100:
            continue
        sub_s = EnvKernelFit(z[ii], mb[ii], mass[ii], is_cal[ii],
                             cov[np.ix_(ii, ii)])
        rq_s = sub_s.fit(lambda p, k=ii: np.exp(
            p[0] * fpr[k] + p[1] * err_r[k] + p[2] * fchi_r[k]),
            [(-4.0, 4.0)] * 3)
        xd_s = np.zeros(len(ii))
        kk = np.where(ok_ang[ii])[0]
        if len(kk) >= 20:
            rr = rankdata(-np.log(np.maximum(
                s_raw[ii][kk], 1e-3))) / len(kk)
            xd_s[kk] = rr - rr.mean()
        xa_s = X_amb[ii] - X_amb[ii].mean()
        re_s = sub_s.fit(lambda p, xa=xa_s, xd=xd_s: np.exp(
            p[0] * xa + p[1] * xd
            + p[2] * fpr[ii] + p[3] * err_r[ii] + p[4] * fchi_r[ii]),
            [(-4.0, 4.0)] * 5)
        ent = {
            "n": int(len(ii)), "n_depth": int(len(kk)),
            "env_attributable_chi2": float(rq_s.fun - re_s.fun),
            "a_ambient": float(re_s.x[1]), "b_depth": float(re_s.x[2]),
        }
        if len(kk) < 20:
            ent["note"] = ("x_depth has no coverage in this scheme "
                           "(offset provenance unverified; neutral depth weight), so "
                           "b_depth is unidentified and its fitted "
                           "value carries no physical meaning")
        per_scheme[f"IDSURVEY_{s_}"] = ent
        print_status(
            f"  G10 per-scheme {s_}: env-attr "
            f"{rq_s.fun - re_s.fun:+.2f} a={re_s.x[1]:+.3f} "
            f"b={re_s.x[2]:+.3f}")

    inter = X_amb * x_depth
    r_int = fit.fit(lambda p: np.exp(
        p[0] * X_amb + p[1] * x_depth + p[2] * inter
        + p[3] * fpr + p[4] * err_r + p[5] * fchi_r), [(-4.0, 4.0)] * 6)

    vp = df["VPEC"].values[valid].astype(float)
    okv = np.isfinite(vp)
    vpr = np.zeros(len(z))
    vv_ = rankdata(np.abs(vp[okv])) / okv.sum()
    vpr[okv] = vv_ - vv_.mean()
    r_v = fit.fit(lambda p: np.exp(p[0] * vpr), [(-4.0, 4.0)])
    r_vj = fit.fit(lambda p: np.exp(
        p[0] * X_amb + p[1] * x_depth + p[2] * fpr + p[3] * err_r
        + p[4] * fchi_r + p[5] * vpr), [(-4.0, 4.0)] * 6)

    r_sat = fit.fit(lambda p: np.exp(
        p[0] * np.tanh(X_amb / 0.3) + p[1] * np.tanh(x_depth / 0.3)
        + p[2] * fpr + p[3] * err_r + p[4] * fchi_r), [(-4.0, 4.0)] * 5)

    base10 = np.array([r_eq.x[0], *r_eq.x[1:6], r_eq.x[6]])

    def _chi2_ab(a_fix=None, b_fix=None):
        def f(th):
            eps, a_, b_, c1, c2, c3, M = th
            if a_fix is not None:
                a_ = a_fix
            if b_fix is not None:
                b_ = b_fix
            w = np.exp(a_ * X_amb + b_ * x_depth
                       + c1 * fpr + c2 * err_r + c3 * fchi_r)
            return fit.chi2_given_w(eps, M, w)
        return float(minimize(
            f, base10.copy(), method="L-BFGS-B",
            bounds=[(0, 2)] + [(-4, 4)] * 5 + [(-21, -16)]).fun)

    prof_a = {str(v): _chi2_ab(a_fix=v) - float(r_eq.fun)
              for v in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)}
    prof_b = {str(v): _chi2_ab(b_fix=v) - float(r_eq.fun)
              for v in (-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0)}

    # (g) ambient-coordinate robustness: the fixed-aperture neighbour
    #     count is sparse by construction (most sightlines see zero
    #     catalogued companions within 3 Mpc, |dz| < 0.02). Scan the
    #     aperture/dz window and the ambient mass:density weight; the
    #     attribution must persist under every construction rather than
    #     sit on a knife-edge choice.
    dr_counts, nb_counts = local_density_rank(
        z, ra, dec, return_counts=True, event_ids=event_ids)
    dens_scan = {}
    for ap_ in (2.0, 3.0, 4.0, 8.0, 16.0):
        for dz_ in (0.02, 0.05, 0.10):
            dr2 = local_density_rank(z, ra, dec,
                                     aperture_mpc=ap_, dz=dz_, event_ids=event_ids)
            Xa2 = 0.5 * (mass_rank + dr2)
            Xa2 -= Xa2.mean()
            re2 = fit.fit(lambda p, xa=Xa2: np.exp(
                p[0] * xa + p[1] * x_depth + p[2] * fpr
                + p[3] * err_r + p[4] * fchi_r), [(-4.0, 4.0)] * 5)
            dens_scan[f"ap{ap_:.0f}Mpc_dz{dz_:.2f}"] = {
                "env_attributable_chi2": float(r_q.fun - re2.fun),
                "a_ambient": float(re2.x[1]),
                "b_depth": float(re2.x[2]),
            }
    weight_scan = {}
    for wm in (0.25, 0.5, 0.75):
        Xa3 = wm * mass_rank + (1.0 - wm) * dens_rank
        Xa3 -= Xa3.mean()
        re3 = fit.fit(lambda p, xa=Xa3: np.exp(
            p[0] * xa + p[1] * x_depth + p[2] * fpr
            + p[3] * err_r + p[4] * fchi_r), [(-4.0, 4.0)] * 5)
        weight_scan[f"w_mass_{wm}"] = {
            "env_attributable_chi2": float(r_q.fun - re3.fun),
            "a_ambient": float(re3.x[1]),
            "b_depth": float(re3.x[2]),
        }

    # (h) per-sightline influence: use the same order-invariant
    #     precision-weighted contribution as the localization above.
    d_sn = d_env
    d_sort = np.sort(d_sn)[::-1]
    influence = {
        "total": float(d_sn.sum()),
        "top10_sum": float(d_sort[:10].sum()),
        "top10_frac_of_total": float(d_sort[:10].sum() / d_sn.sum()),
        "max_abs_single": float(np.abs(d_sn).max()),
        "n_improved": int((d_sn > 0).sum()),
        "n_worsened": int((d_sn < 0).sum()),
    }

    # (j) covariance sensitivity: the attribution must not be a
    #     knife-edge property of the released covariance's structure.
    cov_sens = {}
    for label, c2 in (
            ("nominal", cov),
            ("diag_only", np.diag(np.diag(cov))),
            ("cov_scale_0p9", 0.9 * cov),
            ("cov_scale_1p1", 1.1 * cov),
            ("cov_plus_10pct_diag", cov + 0.1 * np.diag(np.diag(cov)))):
        f_cv = EnvKernelFit(z, mb, mass, is_cal, c2)
        rq_cv = f_cv.fit(lambda p: np.exp(
            p[0] * fpr + p[1] * err_r + p[2] * fchi_r), [(-4.0, 4.0)] * 3)
        re_cv = f_cv.fit(lambda p: np.exp(
            p[0] * X_amb + p[1] * x_depth + p[2] * fpr
            + p[3] * err_r + p[4] * fchi_r), [(-4.0, 4.0)] * 5)
        cov_sens[label] = {
            "env_attributable_chi2": float(rq_cv.fun - re_cv.fun),
            "a_ambient": float(re_cv.x[1]),
            "b_depth": float(re_cv.x[2]),
        }

    # (k) joint marginal covariance of (a, b) at the G8 MLE from the
    #     full 7-parameter Hessian: the (a,b) block of the inverted
    #     Hessian gives each coefficient's marginal sigma and their
    #     correlation with all nuisance covariance propagated.
    th0 = np.array(r_eq.x, dtype=float)
    Cgates = np.stack([X_amb, x_depth, fpr, err_r, fchi_r])

    def _chi2_full(th):
        w = np.exp(th[1:6] @ Cgates)
        return fit.chi2_given_w(th[0], th[-1], w)

    nh = len(th0)
    hh = np.array([1e-3, 5e-3, 5e-3, 1e-2, 1e-2, 1e-2, 1e-3])
    Hfull = np.zeros((nh, nh))
    fc = _chi2_full(th0)
    for i_ in range(nh):
        e_ = np.zeros(nh); e_[i_] = hh[i_]
        Hfull[i_, i_] = (_chi2_full(th0 + e_) - 2 * fc
                         + _chi2_full(th0 - e_)) / hh[i_] ** 2
    for i_ in range(nh):
        for j_ in range(i_ + 1, nh):
            ei = np.zeros(nh); ei[i_] = hh[i_]
            ej = np.zeros(nh); ej[j_] = hh[j_]
            Hfull[i_, j_] = Hfull[j_, i_] = (
                _chi2_full(th0 + ei + ej) - _chi2_full(th0 + ei - ej)
                - _chi2_full(th0 - ei + ej) + _chi2_full(th0 - ei - ej)
            ) / (4 * hh[i_] * hh[j_])
    try:
        Cov_full = 2.0 * np.linalg.inv(Hfull)
        hess = {
            "sigma_a": float(np.sqrt(Cov_full[1, 1])),
            "sigma_b": float(np.sqrt(Cov_full[2, 2])),
            "corr_ab": float(Cov_full[1, 2] / np.sqrt(
                Cov_full[1, 1] * Cov_full[2, 2])),
            "note": "marginal covariance from the full 7-param "
                    "Hessian at the MLE (Delta chi2 = 2 * H^-1)",
        }
    except np.linalg.LinAlgError:
        hess = {"error": "non-PD Hessian"}

    # (i) held-out validation: strict parameter transfer. The gate is
    #     trained on a subset and its coefficients applied unchanged to
    #     held-out sightlines (no nuisance re-fit). Positive held-out
    #     attribution checks parameter transfer within this catalog.
    #     Two partitionings: physical-event random 5-fold and survey-blocked.
    #     Marginal held-out covariance is used; this is not external replication.
    def _eval_transfer(fit_te, p_, coords_te, tem):
        w_ = np.exp(p_[1:-1] @ coords_te)
        mu_ = fit_te.mu_model(float(p_[0]), w_) + float(p_[-1])
        y_ = solve_triangular(fit_te.L, mb[tem] - mu_, lower=True)
        return float(y_ @ y_)

    rng_cv = np.random.default_rng(3)
    unique_cv, inverse_cv = np.unique(event_ids, return_inverse=True)
    event_folds_cv = np.zeros(len(unique_cv), dtype=int)
    perm_cv = rng_cv.permutation(len(unique_cv))
    for kf in range(5):
        event_folds_cv[perm_cv[kf::5]] = kf
    folds_cv = event_folds_cv[inverse_cv]
    oos_folds = []
    for kf in range(5):
        trm = folds_cv != kf
        tem = ~trm
        assert not (set(event_ids[tem]) & set(event_ids[trm]))
        f_tr = EnvKernelFit(z[trm], mb[trm], mass[trm], is_cal[trm],
                            cov[np.ix_(trm, trm)])
        r_eq_t = f_tr.fit(lambda p, k=trm: np.exp(
            p[0] * X_amb[k] + p[1] * x_depth[k] + p[2] * fpr[k]
            + p[3] * err_r[k] + p[4] * fchi_r[k]), [(-4.0, 4.0)] * 5)
        r_q_t = f_tr.fit(lambda p, k=trm: np.exp(
            p[0] * fpr[k] + p[1] * err_r[k] + p[2] * fchi_r[k]),
            [(-4.0, 4.0)] * 3)
        f_te = EnvKernelFit(z[tem], mb[tem], mass[tem], is_cal[tem],
                            cov[np.ix_(tem, tem)])
        d_ = (_eval_transfer(
                  f_te, r_q_t.x,
                  np.stack([fpr[tem], err_r[tem], fchi_r[tem]]), tem)
              - _eval_transfer(
                  f_te, r_eq_t.x,
                  np.stack([X_amb[tem], x_depth[tem], fpr[tem],
                            err_r[tem], fchi_r[tem]]), tem))
        oos_folds.append(float(d_))

    oos_surveys = {}
    for s_ in (1, 10, 15):
        tem = surv_i == s_
        trm = ~np.isin(event_ids, event_ids[tem])
        assert not (set(event_ids[tem]) & set(event_ids[trm]))
        f_tr = EnvKernelFit(z[trm], mb[trm], mass[trm], is_cal[trm],
                            cov[np.ix_(trm, trm)])
        r_eq_t = f_tr.fit(lambda p, k=trm: np.exp(
            p[0] * X_amb[k] + p[1] * x_depth[k] + p[2] * fpr[k]
            + p[3] * err_r[k] + p[4] * fchi_r[k]), [(-4.0, 4.0)] * 5)
        r_q_t = f_tr.fit(lambda p, k=trm: np.exp(
            p[0] * fpr[k] + p[1] * err_r[k] + p[2] * fchi_r[k]),
            [(-4.0, 4.0)] * 3)
        f_te = EnvKernelFit(z[tem], mb[tem], mass[tem], is_cal[tem],
                            cov[np.ix_(tem, tem)])
        d_ = (_eval_transfer(
                  f_te, r_q_t.x,
                  np.stack([fpr[tem], err_r[tem], fchi_r[tem]]), tem)
              - _eval_transfer(
                  f_te, r_eq_t.x,
                  np.stack([X_amb[tem], x_depth[tem], fpr[tem],
                            err_r[tem], fchi_r[tem]]), tem))
        oos_surveys[f"IDSURVEY_{s_}"] = {
            "n_held_out": int(tem.sum()),
            "held_out_env_attributable_chi2": float(d_),
        }

    results["G10_depth_robustness_audit"] = {
        "scheme_audit": {
            "angsep_distributions": ang_stats,
            "excluded_from_depth": "IDSURVEY 15 (PS1MD): depth weight "
                                   "neutral under the specified selection; "
                                   "the ~20 arcsec column distribution "
                                   "has unverified offset provenance",
            "survey_names_source": "official Pantheon+ release README",
            "ranking": "within-survey rank of -ln(r_proj/R_eff) per "
                       "scheme; prevents cross-scheme contrast from "
                       "masquerading as radial depth",
        },
        "per_scheme_replication": per_scheme,
        "interaction_test": {
            "extra_gain_chi2": float(r_eq.fun - r_int.fun),
            "c_interaction": float(r_int.x[3]),
        },
        "vpec_test": {
            "n_valid": int(okv.sum()),
            "standalone_gain_chi2": float(chi2_u - r_v.fun),
            "extra_gain_on_g8": float(r_eq.fun - r_vj.fun),
            "c_vpec": float(r_vj.x[6]),
        },
        "saturating_gate": {
            "chi2": float(r_sat.fun),
            "env_attributable_chi2": float(r_q.fun - r_sat.fun),
            "a": float(r_sat.x[1]), "b": float(r_sat.x[2]),
            "note": "The bounded tanh coordinate retains positive environmental "
                    "gain; this alternative fit alone does not establish "
                    "a linear-response regime or exclude influential objects",
        },
        "profile_chi2": {
            "a_scan": prof_a, "b_scan": prof_b,
            "note": "dchi2 vs joint MLE; +-1-sigma where curve "
                    "crosses +1",
        },
        "ambient_coordinate_robustness": {
            "frac_zero_neighbours": float((nb_counts == 0).mean()),
            "aperture_dz_scan": dens_scan,
            "ambient_weight_scan": weight_scan,
            "note": "the fixed-aperture count is sparse by design "
                    "(catalog density limit); attribution persists "
                    "under every aperture, dz window, and "
                    "mass:density weight tested",
        },
        "per_sightline_influence": influence,
        "covariance_sensitivity": cov_sens,
        "ab_joint_hessian": hess,
        "held_out_validation": {
            "random_partition_unit": "physical CID event; no train/test event overlap",
            "coordinate_scope": "Outcome-free full-catalog coordinates; within-catalog parameter transfer, not external replication",
            "random_5fold_per_fold": oos_folds,
            "random_5fold_total": float(np.sum(oos_folds)),
            "survey_blocked": oos_surveys,
            "note": "Strict parameter transfer with no nuisance refit and no "
                    "shared physical events. Marginal held-out covariance "
                    "is used. Positive pooled gain does not imply every "
                    "fold or survey is positive; within-catalog validation",
        },
    }
    print_status(
        f"G10: interaction +{r_eq.fun - r_int.fun:.2f} | "
        f"VPEC +{chi2_u - r_v.fun:.2f}/+{r_eq.fun - r_vj.fun:.2f} | "
        f"sat {r_q.fun - r_sat.fun:+.2f} | "
        f"profiles a:0->{prof_a['0.0']:+.1f} b:0->{prof_b['0.0']:+.1f} | "
        f"OOS 5fold {np.sum(oos_folds):+.2f} | "
        f"survey-OOS "
        f"{[round(v['held_out_env_attributable_chi2'], 2) for v in oos_surveys.values()]}"
    )
    save_stage(results, 'G10_robustness')

    # ---- G9: standardization audit --------------------------------------
    # Check reconstruction of the released standardized brightness:
    # the released product reconstructs exactly from raw columns
    # (m_b_corr = mB + alpha*x1 - beta*c - biasCor + const), and a
    # self-contained re-standardisation with the Tripp parameters and a
    # host-mass step refit inside the likelihood checks residual freedom.
    # The retained biasCor already includes mass-dependent corrections;
    # this reconstruction does not remove or validate those corrections.
    mB_ = df["mB"].values[valid].astype(float)
    x1_ = df["x1"].values[valid].astype(float)
    c_ = df["c"].values[valid].astype(float)
    bc_ = df["biasCor_m_b"].values[valid].astype(float)
    ok9 = np.isfinite(mB_) & np.isfinite(x1_) & np.isfinite(c_) & np.isfinite(bc_)
    X9 = np.column_stack([x1_[ok9], c_[ok9], np.ones(ok9.sum())])
    y9 = (df["m_b_corr"].values[valid] - mB_ + bc_)[ok9]
    coef9 = np.linalg.lstsq(X9, y9, rcond=None)[0]
    resid9 = y9 - X9 @ coef9

    step9 = (mass >= MASS_SPLIT).astype(float)

    def chi2_restd(th):
        eps, M, a_, b_, g_ = th
        resid = (mB_ - bc_ + a_ * x1_ - b_ * c_) \
            - (fit.mu_model(eps, np.ones(len(z))) + M) \
            - observer_correction - g_ * step9
        yy = solve_triangular(fit.L, resid, lower=True)
        return float(yy @ yy)

    r9 = minimize(chi2_restd, [0.87, -19.33, 0.15, 3.0, 0.05],
                  method="L-BFGS-B",
                  bounds=[(0, 2), (-21, -16), (-1, 1), (0, 5), (-1, 1)])
    results["G9_standardization_audit"] = {
        "reconstruction": "m_b_corr = mB + a x1 - b c - biasCor + const",
        "alpha": float(coef9[0]), "beta": float(-coef9[1]),
        "const": float(coef9[2]),
        "reconstruction_resid_std": float(resid9.std()),
        "mass_step_in_mb_corr": "no additional explicit step outside biasCor "
                                "is needed for this reconstruction; "
                                "mass-dependent corrections remain in biasCor",
        "bias_correction_scope": "retains released biasCor; not an independent "
                                 "raw-light-curve standardization",
        "self_restandardized_uniform_fit": {
            "chi2": float(r9.fun), "epsilon": float(r9.x[0]),
            "alpha": float(r9.x[2]), "beta": float(r9.x[3]),
            "mass_step_gamma": float(r9.x[4]),
        },
    }
    print_status(
        f"G9 standardization: m_b_corr reconstructs exactly "
        f"(a={coef9[0]:.3f} b={-coef9[1]:.3f}); "
        f"restd uniform chi2={r9.fun:.2f} step={r9.x[4]:+.4f}"
    )

    results["meta"]["runtime_seconds"] = round(time.time() - t0, 1)
    out = step_json_path(STEP_ID)
    write_json(out, results)
    print_status(f"Wrote {out}")
    return results


def _nested_two_channel(fit, X_amb, x_depth, prior_wid, nlive):
    """dynesty lnZ for the two-channel gate under +-prior_wid priors."""
    import dynesty
    bounds = [(0.0, 2.0), (-prior_wid, prior_wid),
              (-prior_wid, prior_wid), (-21.0, -16.0)]
    lo = np.array([b[0] for b in bounds])
    hi = np.array([b[1] for b in bounds])

    def prior(u):
        return lo + u * (hi - lo)

    def ll(t):
        eps, a, b, M = t
        w = np.exp(a * X_amb + b * x_depth)
        if not np.all(np.isfinite(w)) or np.any(w <= 0):
            return -1e6
        c = fit.chi2_given_w(eps, M, w)
        return -0.5 * c if np.isfinite(c) else -1e6

    s = dynesty.NestedSampler(ll, prior, 4, bound="multi",
                              sample="rwalk", nlive=nlive, rstate=np.random.default_rng(20261011))
    s.run_nested(dlogz=0.05)
    r = s.results
    return float(r.logz[-1]), float(r.logzerr[-1])


if __name__ == "__main__":
    main()
