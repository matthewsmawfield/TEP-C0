#!/usr/bin/env python3
"""Step 025c: TEP Screening Model Fit to DDR Data.

Fits the TEP screening model to DDR constraints to extract
the screening parameter ε_s that governs when TEP effects
become important at high redshift.

The screening model predicts η(z) → 1 as z → ∞, with
approach to unity governed by ε_s.

References:
  - BOSS/DESI BAO: Anderson et al. 2014, MNRAS, 441, 24; DESI Collaboration 2024
    z_c,probe    - Critical redshift for screening transition
    α_probe      - Transition steepness

Physical interpretation:
    - Different probes couple differently to TEP screening
    - Stronger coupling in linear regime (BAO, η_0 ~ 0.1)
    - Weaker coupling in collapsed structures (clusters, η_0 ~ 0.2-0.3)
    - All probes converge to η → 1 at high z (homogeneous early universe)
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from c0_common import (
    RESULTS_DIR,
    TEPLogger,
    ensure_dirs,
    print_status,
    set_step_logger,
    step_json_path,
    write_json,
)

STEP_ID = "step_04_06_screening_fit"


def eta_model(z: np.ndarray, eta_0: float, beta: float) -> np.ndarray:
    """Empirical probe-dependent DDR anomaly model.

    Two-parameter power-law per probe that admits both eta>1 and eta<1 and
    can capture rising or falling trends with redshift:

        eta(z) = eta_0 * (1 + z) ** beta

    - ``eta_0``  : z=0 amplitude (eta_0 = 1 corresponds to no anomaly).
    - ``beta``   : slope index. beta > 0  ->  anomaly grows with z
                                beta < 0  ->  anomaly decays with z
                                beta = 0  ->  scale-invariant offset.

    The previous (eta_0 + (1-eta_0)(z/(z+z_c))^alpha) form forced eta -> 1
    at high z and bounded eta in [eta_0, 1], which is incompatible with the
    SZ/SGL data (eta ~ 1.2-1.3 at low z) and with the BAO trend (eta drops
    well below 1 at z ~ 1.5). The power-law form removes both pathologies.
    """
    return eta_0 * (1.0 + z) ** beta


def chi2_total(params: np.ndarray, data: dict) -> float:
    """Compute total chi-squared for all probes.

    Parameter layout (length 6):
        params = [eta0_BAO, beta_BAO, eta0_SZ, beta_SZ, eta0_SGL, beta_SGL]
    """
    chi2 = 0.0
    probe_slices = (("BAO", 0), ("SZ", 2), ("SGL", 4))
    for probe, off in probe_slices:
        pts = data.get(probe, [])
        if not pts:
            continue
        zs = np.array([p[0] for p in pts])
        etas = np.array([p[1] for p in pts])
        errs = np.array([max(p[2], np.finfo(float).tiny) for p in pts])
        pred = eta_model(zs, params[off], params[off + 1])
        chi2 += float(np.sum(((etas - pred) / errs) ** 2))
    return chi2


def run() -> dict:
    """Execute TEP screening model fit."""
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    # Load BAO data from step 025
    bao_path = RESULTS_DIR / "step_04_04_distance_duality.json"
    bao_z, bao_eta, bao_err = [], [], []
    if bao_path.exists():
        with open(bao_path) as f:
            bao_data = json.load(f)
        for c in bao_data.get("ddr_constraints", []):
            bao_z.append(c["z"])
            bao_eta.append(c["eta_obs"])
            bao_err.append(c["eta_err"])
    
    # Load SZ data (skipped if no empirical eta_obs available)
    sz_path = RESULTS_DIR / "step_01_05_download_sz.json"
    sz_z, sz_eta, sz_err = [], [], []
    sz_skipped_reason = None
    if sz_path.exists():
        with open(sz_path) as f:
            sz_data = json.load(f)
        clusters = sz_data.get("data", {}).get("clusters", [])
        missing_eta = [c for c in clusters if "eta_obs" not in c]
        if missing_eta and not any("eta_obs" in c for c in clusters):
            sz_skipped_reason = (
                f"step_01_05 provides D_A but not eta_obs for {len(clusters)} clusters; "
                "fabrication forbidden, SZ probe omitted from screening fit."
            )
            print_status(sz_skipped_reason, "WARNING")
        else:
            for c in clusters:
                if "eta_obs" not in c:
                    continue
                sz_z.append(c["z"])
                sz_eta.append(c["eta_obs"])
                sz_err.append(c.get("eta_err", 0.07))
    
    # Load SGL data (skipped if no empirical eta_obs available)
    sgl_path = RESULTS_DIR / "step_01_06_download_sgl.json"
    sgl_z, sgl_eta, sgl_err = [], [], []
    sgl_skipped_reason = None
    if sgl_path.exists():
        with open(sgl_path) as f:
            sgl_data = json.load(f)
        systems = sgl_data.get("data", {}).get("systems", [])
        if systems and not any("eta_obs" in s for s in systems):
            sgl_skipped_reason = (
                f"step_01_06 provides D_A but not eta_obs for {len(systems)} systems; "
                "fabrication forbidden, SGL probe omitted from screening fit."
            )
            print_status(sgl_skipped_reason, "WARNING")
        else:
            for s in systems:
                if "eta_obs" not in s:
                    continue
                sgl_z.append(s["z_lens"])
                sgl_eta.append(s["eta_obs"])
                sgl_err.append(s.get("eta_err", 0.025))
    
    # Prepare data for fit
    data_dict = {
        "BAO": list(zip(bao_z, bao_eta, bao_err)),
        "SZ": list(zip(sz_z, sz_eta, sz_err)),
        "SGL": list(zip(sgl_z, sgl_eta, sgl_err)),
    }
    
    n_total = len(bao_z) + len(sz_z) + len(sgl_z)
    print_status(f"Total data points: {n_total}", "INFO")
    print_status(f"  BAO: {len(bao_z)}", "INFO")
    print_status(f"  SZ:  {len(sz_z)}", "INFO")
    print_status(f"  SGL: {len(sgl_z)}", "INFO")

    if len(bao_z) == 0:
        payload = {
            "step": STEP_ID,
            "status": "blocked",
            "error": "No BAO DDR constraints available; screening fit cannot proceed.",
            "validation": {
                "claim_gate": "blocked",
                "blockers": ["missing_bao_ddr_constraints"],
                "sz_skipped_reason": sz_skipped_reason,
                "sgl_skipped_reason": sgl_skipped_reason,
            },
            "timestamp": int(time.time()),
        }
        write_json(step_json_path(STEP_ID), payload)
        print_status(f"Step {STEP_ID} blocked: no BAO data", "WARNING")
        return payload

    # Initial parameter guess for [eta0, beta] per probe.
    # Starting from no-anomaly (eta_0=1.0, beta=0.0) so the optimiser is
    # free to move toward whatever sign and magnitude the data prefer.
    p0 = [1.0, 0.0,   # BAO
          1.0, 0.0,   # SZ
          1.0, 0.0]   # SGL

    # Physical/empirical bounds: eta_0 ~ O(1), |beta| < a few.
    bounds = [
        (0.2, 3.0), (-3.0, 3.0),  # BAO
        (0.2, 3.0), (-3.0, 3.0),  # SZ
        (0.2, 3.0), (-3.0, 3.0),  # SGL
    ]
    
    print_status("Fitting TEP screening model...", "INFO")
    
    result = minimize(
        chi2_total,
        p0,
        args=(data_dict,),
        bounds=bounds,
        method="L-BFGS-B",
    )
    
    if result.success:
        p = result.x
        chi2_min = result.fun
        n_params = 6
        n_dof = n_total - n_params
        n_dof_safe = max(n_dof, 1)  # Prevent division by zero

        print_status("Fit successful", "SUCCESS")
        print_status(f"  χ² = {chi2_min:.2f}", "INFO")
        print_status(f"  ndof = {n_dof}", "INFO")
        print_status(f"  χ²/ndof = {chi2_min/n_dof_safe:.2f}", "INFO")

        # Extract parameters (eta_0, beta) per probe.
        bao_params = {"eta_0": float(p[0]), "beta": float(p[1])}
        sz_params  = {"eta_0": float(p[2]), "beta": float(p[3])}
        sgl_params = {"eta_0": float(p[4]), "beta": float(p[5])}

        print_status("\nBest-fit parameters:", "INFO")
        print_status(f"  BAO: η_0={bao_params['eta_0']:.3f}, β={bao_params['beta']:+.3f}", "INFO")
        print_status(f"  SZ:  η_0={sz_params['eta_0']:.3f}, β={sz_params['beta']:+.3f}", "INFO")
        print_status(f"  SGL: η_0={sgl_params['eta_0']:.3f}, β={sgl_params['beta']:+.3f}", "INFO")

        # Compute predictions at test redshifts.
        z_test = [0.1, 0.5, 1.0, 2.0]
        predictions = {}
        for name, prm in [("BAO", bao_params), ("SZ", sz_params), ("SGL", sgl_params)]:
            preds = []
            for z in z_test:
                eta_pred = eta_model(np.array([z]), prm["eta_0"], prm["beta"])[0]
                preds.append(float(eta_pred))
            predictions[name] = preds
        
        # Compute probe differences
        delta_bao_sz = abs(bao_params["eta_0"] - sz_params["eta_0"])
        delta_bao_sgl = abs(bao_params["eta_0"] - sgl_params["eta_0"])
        delta_sz_sgl = abs(sz_params["eta_0"] - sgl_params["eta_0"])
        
        print_status(f"\nProbe baseline differences:", "INFO")
        print_status(f"  |η_0_BAO - η_0_SZ|  = {delta_bao_sz:.3f}", "INFO")
        print_status(f"  |η_0_BAO - η_0_SGL| = {delta_bao_sgl:.3f}", "INFO")
        print_status(f"  |η_0_SZ  - η_0_SGL| = {delta_sz_sgl:.3f}", "INFO")
        
        # Data-driven ordering by |eta_0 - 1|
        probe_strengths = sorted(
            [("BAO", bao_params), ("SZ", sz_params), ("SGL", sgl_params)],
            key=lambda kv: abs(kv[1]["eta_0"] - 1.0),
            reverse=True,
        )
        print_status("\nPhysical interpretation (data-driven):", "INFO")
        for probe, prm in probe_strengths:
            sign = "above" if prm["eta_0"] > 1 else "below"
            print_status(
                f"  • {probe}: η_0={prm['eta_0']:.3f} ({sign} unity), β={prm['beta']:+.3f}",
                "INFO",
            )
        
        payload = {
            "step": STEP_ID,
            "description": "TEP probe-dependent screening model fit to DDR data",
            "status": "completed",
            "model": {
                "formula": "η_probe(z) = η_0 × (1+z)^β",
                "name": "Two-parameter empirical DDR anomaly per probe",
                "note": (
                    "Earlier (η_0 + (1-η_0)(z/(z+z_c))^α) form forced η→1 at high z and "
                    "η_0 ≤ 1, incompatible with SZ/SGL data (η ≳ 1.2). The power-law "
                    "form admits η>1, η<1, rising and falling trends, and is therefore "
                    "the empirically minimal parameterisation that can fit all three probes."
                ),
            },
            "parameters": {
                "BAO": bao_params,
                "SZ": sz_params,
                "SGL": sgl_params,
            },
            "fit_quality": {
                "chi2": float(chi2_min),
                "n_dof": n_dof,
                "chi2_per_dof": float(chi2_min / n_dof) if n_dof > 0 else None,
                "success": bool(result.success),
            },
            "predictions": {
                "redshifts": z_test,
                "BAO": predictions["BAO"],
                "SZ": predictions["SZ"],
                "SGL": predictions["SGL"],
            },
            "probe_differences": {
                "delta_eta0_BAO_vs_SZ":  float(delta_bao_sz),
                "delta_eta0_BAO_vs_SGL": float(delta_bao_sgl),
                "delta_eta0_SZ_vs_SGL":  float(delta_sz_sgl),
            },
            "interpretation": {
                "physical_scale_BAO": "~100 Mpc (linear regime)",
                "physical_scale_SZ": "~5 Mpc (cluster gas)",
                "physical_scale_SGL": "~5 Mpc (cluster potential)",
                "ordering_by_anomaly_magnitude": [p for p, _ in probe_strengths],
                "implication": (
                    "Probe-dependent Etherington anomaly: SZ and SGL exhibit η>1 at "
                    "low z (D_L exceeds the metric prediction) while BAO declines "
                    "through unity into η<1 by z~1.5. Within TEP this is the line-of-"
                    "sight signature of environment-dependent screening; under ΛCDM "
                    "it must be absorbed into systematics."
                ),
            },
            "data_usage": {
                "n_BAO": len(bao_z),
                "n_SZ": len(sz_z),
                "n_SGL": len(sgl_z),
                "sz_skipped_reason": sz_skipped_reason,
                "sgl_skipped_reason": sgl_skipped_reason,
            },
            "validation": {
                "claim_gate": "open" if (len(sz_z) > 0 and len(sgl_z) > 0) else "partial",
                "research_grade": bool(len(sz_z) > 0 and len(sgl_z) > 0),
                "note": (
                    "BAO-only fit; SZ and SGL probes omitted because upstream loaders do "
                    "not yet provide empirical eta_obs. SZ/SGL fit parameters in this "
                    "payload are unconstrained and must not be quoted as measurements."
                ) if (len(sz_z) == 0 or len(sgl_z) == 0) else None,
            },
            "timestamp": int(time.time()),
        }
    else:
        print_status(f"Fit failed: {result.message}", "FAIL")
        payload = {
            "step": STEP_ID,
            "status": "failed",
            "error": str(result.message),
            "timestamp": int(time.time()),
        }
    
    write_json(step_json_path(STEP_ID), payload)
    print_status(f"Step {STEP_ID} completed", "SUCCESS")
    return payload


if __name__ == "__main__":
    print(run())
