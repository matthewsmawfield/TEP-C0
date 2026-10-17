#!/usr/bin/env python3
"""Audit verification: corrected-pipeline Bayes factors.

Runs the pipeline's own model classes (step_03_01 ModelLCDM / ModelTEP)
under the corrected conventions (zHD, n_T=2) at research-grade nested
sampling settings (nlive=500, dlogz=0.1), plus isolation runs that vary
one convention at a time so the contribution of each corrected bug is
visible.

Grid:
  corrected  : zHD + n_T=2   -> full, z>0.023, no_calibrators
  zcmb_iso   : zCMB + n_T=2  -> full          (isolates redshift bug)
  nt1_iso    : zHD + n_T=1   -> full          (isolates n_T bug)
  legacy     : zCMB + n_T=1  -> full          (reproduces pre-fix run)

Writes results/audit_corrected_evidence.json.
"""

import os
import sys
import json
import time
import numpy as np
from pathlib import Path

os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("TEP_DYNESTY_PROGRESS", "0")

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "steps"))

import step_03_01_three_model_comparison as s301
from scipy.linalg import cholesky

NLIVE = int(os.getenv("AUDIT_NLIVE", "500"))
DLOGZ = float(os.getenv("AUDIT_DLOGZ", "0.1"))
OUT = ROOT / "results" / "audit_corrected_evidence.json"


class ModelTEPnT(s301.ModelTEP):
    """ModelTEP with explicit n_T control (everything else identical)."""

    def __init__(self, z_T=5.0, n_T=2.0, free_z_T=False):
        super().__init__(pure_shear=False, free_z_T=free_z_T)
        self.n_T = n_T
        if not free_z_T:
            self.z_T = z_T

    def predict(self, z, params, data=None):
        M = params.get("M", 0.0)
        eps = params["epsilon_shear_los"]
        z_T = params.get("z_T", self.z_T) if self.free_z_T else self.z_T
        tep = s301.TEPCosmology(
            H0=self.H0_ref, Omega_m=1.0, epsilon_T=eps, z_T=z_T, n_T=self.n_T
        )
        return tep.distance_modulus(z) + M


def subset_data(data, mask):
    """Return a PantheonData restricted to mask (covariance sub-matrix)."""
    d = s301.PantheonData()
    d.z = data.z[mask]
    d.zhel = None if data.zhel is None else data.zhel[mask]
    d.mb = data.mb[mask]
    d.dmb = data.dmb[mask]
    d.cov = data.cov[np.ix_(mask, mask)]
    d.cov_cholesky = cholesky(d.cov, lower=True, check_finite=False)
    d.cov_logdet = float(2.0 * np.sum(np.log(np.diag(d.cov_cholesky))))
    return d


def calibrator_mask(zcol):
    """IS_CALIBRATOR flags aligned to the filtered PantheonData arrays."""
    import pandas as pd

    df = pd.read_csv(ROOT / "data/raw/pantheon_plus_shoes.dat", sep=r"\s+", comment="#")
    z_raw = df[zcol].values
    mb_raw = df["m_b_corr"].values
    valid = np.isfinite(z_raw) & np.isfinite(mb_raw) & (z_raw > 0.001) & (z_raw < 3.0)
    return df["IS_CALIBRATOR"].values[valid].astype(bool)


def run_case(tag, model, data):
    t0 = time.time()
    mle, logl, _cov = s301.fit_mle(model, data)
    resid = data.mb - model.predict(data.z, mle) - data.observer_correction
    chi2 = data.chi2(resid)
    ev = s301.run_nested_evidence(model, data, nlive=NLIVE, dlogz=DLOGZ)
    dt = time.time() - t0
    out = {
        "model": tag,
        "n_sne": int(len(data.z)),
        "mle": {k: float(v) for k, v in mle.items()},
        "chi2_mle": float(chi2),
        "logZ": ev["log_evidence"],
        "logZ_err": ev["log_evidence_error"],
        "ncall": ev.get("nested_ncall"),
        "runtime_s": round(dt, 1),
    }
    print(f"[done] {tag}: chi2={chi2:.2f} logZ={ev['log_evidence']:.3f}+/-{ev['log_evidence_error']:.3f} ({dt:.0f}s)", flush=True)
    return out


def main():
    record = {
        "description": "Audit rerun of step_03_01 evidence under corrected conventions",
        "nlive": NLIVE,
        "dlogz": DLOGZ,
        "configs": {},
    }

    for zcol, nt, cfg in [
        ("zHD", 2.0, "corrected"),
        ("zCMB", 2.0, "zcmb_iso"),
        ("zHD", 1.0, "nt1_iso"),
        ("zCMB", 1.0, "legacy"),
    ]:
        os.environ["PANTHEON_REDSHIFT_COL"] = zcol
        data = s301.PantheonData()
        data.load()
        assert data.redshift_column_used == zcol

        subsets = {"full": np.ones(len(data.z), dtype=bool)}
        if cfg == "corrected":
            subsets["z_gt_0.023"] = data.z > 0.023
            cal = calibrator_mask("zHD")
            subsets["no_calibrators"] = ~cal

        cfg_res = {"redshift_column": zcol, "n_T": nt, "subsets": {}}
        for sub_name, mask in subsets.items():
            d_sub = subset_data(data, mask)
            sub_res = {}
            m0 = s301.ModelLCDM()
            sub_res["M0a_LCDM"] = run_case(f"{cfg}/{sub_name}/M0a", m0, d_sub)
            m1 = ModelTEPnT(z_T=5.0, n_T=nt)
            sub_res["M1_zT5"] = run_case(f"{cfg}/{sub_name}/M1_zT5", m1, d_sub)
            if cfg == "corrected" and sub_name == "full":
                m1f = ModelTEPnT(n_T=nt, free_z_T=True)
                sub_res["M1_free_zT"] = run_case(f"{cfg}/{sub_name}/M1_free_zT", m1f, d_sub)

            dlnz = sub_res["M1_zT5"]["logZ"] - sub_res["M0a_LCDM"]["logZ"]
            sub_res["BF_M1_zT5_vs_M0a"] = float(np.exp(dlnz))
            sub_res["dlnZ"] = float(dlnz)
            sub_res["delta_chi2_mle"] = sub_res["M1_zT5"]["chi2_mle"] - sub_res["M0a_LCDM"]["chi2_mle"]
            if "M1_free_zT" in sub_res:
                dlnz_f = sub_res["M1_free_zT"]["logZ"] - sub_res["M0a_LCDM"]["logZ"]
                sub_res["BF_M1_free_zT_vs_M0a"] = float(np.exp(dlnz_f))
            cfg_res["subsets"][sub_name] = sub_res

        record["configs"][cfg] = cfg_res
        OUT.write_text(json.dumps(record, indent=2))
        print(f"[checkpoint] wrote {OUT}", flush=True)

    print("all done")


if __name__ == "__main__":
    main()
