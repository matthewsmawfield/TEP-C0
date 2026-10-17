#!/usr/bin/env python3
"""Decisive sampler-configuration test for the shared-intercept evidence.

The main pipeline (step_03_01, bound='multi', sample='rwalk') reports
ln BF ~ 1.87 for M1_zT5 vs M0a, while audit_anchored_likelihood.py
(dynesty defaults: bound='single', sample='auto') reports ~1.22 for the
identical likelihood, data, and priors. This script runs both sampler
configurations at nlive=1000, dlogz=0.05 to determine whether the
difference is a systematic bounding bias or realization noise.
"""

import os, sys, json, time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "steps"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import importlib
m = importlib.import_module("step_03_01_three_model_comparison")

from scipy.linalg import cholesky, solve_triangular

NLIVE = int(os.getenv("TEP_VERIFY_NLIVE", "1000"))
DLOGZ = float(os.getenv("TEP_VERIFY_DLOGZ", "0.05"))

data = m.PantheonData()
assert data.load()
z, mb, cov = data.z, data.mb, data.cov

L = cholesky(cov, lower=True, check_finite=False)
logdet = float(2.0 * np.sum(np.log(np.diag(L))))
n = len(mb)

def loglike_resid(r):
    y = solve_triangular(L, r, lower=True, check_finite=False)
    return float(-0.5 * (y @ y + logdet + n * np.log(2.0 * np.pi)))

models = {
    "M0a_LCDM": (m.ModelLCDM(), ["Om0"], [(0.05, 0.9), (-21.0, -16.0)]),
    "M1_zT5":   (m.ModelTEP(pure_shear=False), ["epsilon_shear_los"], [(0.0, 2.0), (-21.0, -16.0)]),
}

def resid_shared(model, shape_names, x):
    ns = len(shape_names)
    mu = model.predict(z, {**dict(zip(shape_names, x[:ns])), "M": 0.0}, data) + x[ns]
    return mb - mu

import dynesty

def run(model, snames, bounds, bound, sample, tag):
    b = np.array(bounds, dtype=float)
    lo, w = b[:, 0], b[:, 1] - b[:, 0]
    def pt(u):
        return lo + u * w
    def ll(u):
        return loglike_resid(resid_shared(model, snames, u))
    t0 = time.time()
    s = dynesty.NestedSampler(ll, pt, len(bounds), nlive=NLIVE, bound=bound, sample=sample)
    s.run_nested(dlogz=DLOGZ, print_progress=False)
    dt = time.time() - t0
    res = {
        "logZ": float(s.results.logz[-1]),
        "logZ_err": float(s.results.logzerr[-1]),
        "ncall": int(np.sum(s.results.ncall)),
        "eff": float(s.results.eff),
        "minutes": round(dt / 60, 1),
    }
    print(f"{tag}: lnZ={res['logZ']:.3f}+-{res['logZ_err']:.3f} ncall={res['ncall']} ({res['minutes']}m)", flush=True)
    return res

out = {"nlive": NLIVE, "dlogz": DLOGZ, "runs": {}}
for cfg_name, bound, sample in [("multi_rwalk", "multi", "rwalk"), ("single_auto", "single", "auto")]:
    for mname, (model, snames, bounds) in models.items():
        tag = f"{mname}__{cfg_name}"
        out["runs"][tag] = run(model, snames, bounds, bound, sample, tag)

for cfg in ["multi_rwalk", "single_auto"]:
    dz = out["runs"][f"M1_zT5__{cfg}"]["logZ"] - out["runs"][f"M0a_LCDM__{cfg}"]["logZ"]
    out[f"dlnZ_{cfg}"] = dz
    out[f"BF_{cfg}"] = float(np.exp(dz))
    print(f"{cfg}: dlnZ={dz:.3f}  BF={np.exp(dz):.2f}", flush=True)

Path("results/verify_evidence_samplers.json").write_text(json.dumps(out, indent=2))
print("Wrote results/verify_evidence_samplers.json", flush=True)
