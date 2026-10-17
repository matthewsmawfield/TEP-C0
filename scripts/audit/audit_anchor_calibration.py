#!/usr/bin/env python3
"""Audit calibration: is the anchor-removal BF collapse a null test?

The corrected evidence run reports BF_M1_zT5_vs_M0a falling from ~4.3
(full) to ~1.17 after removing the 77 IS_CALIBRATOR anchors and to ~0.66
on the z>0.023 subset (330 objects removed). Whether that collapse is
informative depends on two uncalibrated questions:

A. Does removing an equal number of *random* non-anchor SNe collapse the
   BF the same way?  If yes -> the loss is a calibration/degeneracy
   effect of the subset size, not a calibrator-specific signal.
   If no  -> the evidence localizes in the anchor block.

B. Under TEP-true synthetic data (TEP distances + real covariance
   noise), what BF does the *same* anchor-removed subset return?  If
   TEP-true mocks recover BF>1 on the subset, the observed 0.66 is
   adverse; if they collapse similarly, the subset test is powerless
   and 0.66 carries no information either way.

   GR-true injection on the full sample provides the other side of the
   calibration: whether the observed full-sample BF sits off the
   GR-true null distribution.

Method: identical machinery to verify_corrected_pipeline.py (zHD,
n_T=2, dynesty nlive/dlogz via AUDIT_NLIVE/AUDIT_DLOGZ env).  Writes
results/audit_anchor_calibration.json.
"""

import os
import sys
import json
import time
import numpy as np
from pathlib import Path
from multiprocessing import Pool

os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("TEP_DYNESTY_PROGRESS", "0")

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "steps"))

import step_03_01_three_model_comparison as s301
from verify_corrected_pipeline import subset_data, calibrator_mask, ModelTEPnT
from scipy.linalg import cholesky

NLIVE = int(os.getenv("AUDIT_NLIVE", "500"))
DLOGZ = float(os.getenv("AUDIT_DLOGZ", "0.1"))
N_INJ = int(os.getenv("N_INJ", "12"))
N_RAND = int(os.getenv("N_RAND", "8"))
SEED = int(os.getenv("INJ_SEED", "20261011"))
OUT = ROOT / "results" / "audit_anchor_calibration.json"

# TEP-true generating parameters: corrected full-sample MLE
TEP_TRUE = {"epsilon_shear_los": 0.8735737701811535, "M": -19.331609913021715}
GR_TRUE = {"Om0": 0.3612049791988382, "M": -19.33935146362162}

_DATA = None


def _init_worker():
    global _DATA
    d = s301.PantheonData()
    d.load()
    _DATA = d


def _ev(task):
    """task = (model_tag, z, zhel, mb, dmb, cov2d) -> (tag, logZ, logZ_err)"""
    model_tag, z, zhel, mb, dmb, cov = task
    d = s301.PantheonData()
    d.z = z
    d.zhel = zhel
    d.mb = mb
    d.dmb = dmb
    d.cov = cov
    d.cov_cholesky = cholesky(cov, lower=True, check_finite=False)
    d.cov_logdet = float(2.0 * np.sum(np.log(np.diag(d.cov_cholesky))))
    if model_tag == "M0a_LCDM":
        model = s301.ModelLCDM()
    elif model_tag == "M1_zT5":
        model = ModelTEPnT(z_T=5.0, n_T=2.0)
    else:
        raise ValueError(model_tag)
    ev = s301.run_nested_evidence(model, d, nlive=NLIVE, dlogz=DLOGZ)
    return model_tag, ev["log_evidence"], ev["log_evidence_error"], ev.get("nested_ncall")


def _bf(ev_lcdm, ev_tep):
    return float(np.exp(ev_tep[1] - ev_lcdm[1]))


def main():
    t_start = time.time()
    base = s301.PantheonData()
    base.load()
    rng = np.random.default_rng(SEED)

    cal_mask = calibrator_mask("zHD")
    anchor_window = base.z <= 0.023
    n_cal = int(cal_mask.sum())
    n_low = int(anchor_window.sum())
    print(f"[setup] N={len(base.z)}, calibrators={n_cal}, z<=0.023 removed={n_low}")

    # ---------- truth vectors ----------
    tep = ModelTEPnT(z_T=5.0, n_T=2.0)
    mu_tep = tep.predict(base.z, TEP_TRUE) + base.observer_correction
    lcdm = s301.ModelLCDM()
    mu_gr = lcdm.predict(base.z, GR_TRUE) + base.observer_correction

    L = base.cov_cholesky

    jobs = []   # list of (label, task)
    labels = []

    # ---------- Test A: random-matched removal ----------
    pool_idx = np.where(~cal_mask)[0]           # non-calibrator pool
    flow_idx = np.where(~anchor_window)[0]      # non-anchor-window pool
    for r in range(N_RAND):
        sel = np.sort(rng.choice(pool_idx, n_cal, replace=False))
        keep = np.ones(len(base.z), bool)
        keep[sel] = False
        labels.append((f"rand77_{r}", keep))
    for r in range(N_RAND):
        sel = np.sort(rng.choice(flow_idx, n_low, replace=False))
        keep = np.ones(len(base.z), bool)
        keep[sel] = False
        labels.append((f"rand330_{r}", keep))
    # also: same-count removal drawn ONLY from low-z flow (non-calibrator, z<=0.023)
    lowflow_idx = np.where(anchor_window & ~cal_mask)[0]
    for r in range(N_RAND):
        if len(lowflow_idx) > n_cal:
            sel = np.sort(rng.choice(lowflow_idx, n_cal, replace=False))
        else:
            sel = lowflow_idx
        keep = np.ones(len(base.z), bool)
        keep[sel] = False
        labels.append((f"randlow77_{r}", keep))

    # ---------- Test B: injections ----------
    for r in range(N_INJ):
        g = rng.standard_normal(len(base.z))
        mb_sim = mu_tep + L @ g
        labels.append((f"injTEP_full_{r}", np.ones(len(base.z), bool), mb_sim))
        labels.append((f"injTEP_zcut_{r}", ~anchor_window, mb_sim))
    for r in range(N_INJ):
        g = rng.standard_normal(len(base.z))
        mb_sim = mu_gr + L @ g
        labels.append((f"injGR_full_{r}", np.ones(len(base.z), bool), mb_sim))
        labels.append((f"injGR_zcut_{r}", ~anchor_window, mb_sim))

    # build task list
    for item in labels:
        name, keep = item[0], item[1]
        mb_src = item[2] if len(item) > 2 else base.mb
        for mtag in ("M0a_LCDM", "M1_zT5"):
            jobs.append((name + "|" + mtag,
                         (mtag, base.z[keep], None if base.zhel is None else base.zhel[keep],
                          mb_src[keep], base.dmb[keep], base.cov[np.ix_(keep, keep)])))

    print(f"[run] {len(jobs)} nested-sampling evals, nlive={NLIVE}")
    with Pool(4, initializer=_init_worker) as pool:
        results = pool.map(_ev, [j[1] for j in jobs])

    # collect
    bf = {}
    detail = {}
    for (jname, _), res in zip(jobs, results):
        label, mtag = jname.split("|")
        detail.setdefault(label, {})[mtag] = {
            "logZ": float(res[1]), "logZ_err": float(res[2]), "ncall": res[3]}
    for label, m in detail.items():
        if "M0a_LCDM" in m and "M1_zT5" in m:
            bf[label] = float(np.exp(m["M1_zT5"]["logZ"] - m["M0a_LCDM"]["logZ"]))

    def summ(prefix):
        vals = np.array([v for k, v in bf.items() if k.startswith(prefix)])
        return {
            "n": int(len(vals)),
            "median": float(np.median(vals)),
            "mean": float(vals.mean()),
            "frac_gt1": float(np.mean(vals > 1.0)),
            "q16": float(np.quantile(vals, 0.16)),
            "q84": float(np.quantile(vals, 0.84)),
        }

    out = {
        "description": "Calibration of the anchor-removal BF collapse",
        "nlive": NLIVE, "dlogz": DLOGZ, "seed": SEED,
        "n_calibrators": n_cal, "n_anchor_window": n_low,
        "reference": {
            "full_BF_M1_zT5_vs_M0a": 4.280228490246799,
            "no_calibrators_BF": 1.1683556568180087,
            "z_gt_0.023_BF": 0.657195214000326,
        },
        "summaries": {
            "rand77_remove_random_noncal": summ("rand77_"),
            "rand330_remove_random_nonanchor": summ("rand330_"),
            "randlow77_remove_random_lowz_noncal": summ("randlow77_"),
            "injTEP_full": summ("injTEP_full_"),
            "injTEP_zcut": summ("injTEP_zcut_"),
            "injGR_full": summ("injGR_full_"),
            "injGR_zcut": summ("injGR_zcut_"),
        },
        "per_realization_BF": bf,
        "runtime_s": time.time() - t_start,
    }
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1))
    print(f"[done] {time.time()-t_start:.0f}s -> {OUT}")
    for k, v in out["summaries"].items():
        print(f"  {k:38s} BF med={v['median']:.2f}  [{v['q16']:.2f},{v['q84']:.2f}]  frac>1={v['frac_gt1']:.2f}")


if __name__ == "__main__":
    main()
