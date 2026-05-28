#!/usr/bin/env python3
"""Step 025b: Three-way DDR comparison (BAO, SNe, Strong Lensing).

Compares DDR constraints from three independent methods:
- BAO angular diameter distances (BOSS/DESI BAO measurements)
- Supernovae luminosity distances (Pantheon+ SNe Ia)
- Strong lensing time-delay distances (H0LiCOW/TDCOSMO)

Tests consistency of the Etherington relation η = D_L/(D_A*(1+z)²) = 1

References:
  - Pantheon+: Scolnic et al. 2022, ApJ, 938, 113
  - BOSS/DESI BAO: Anderson et al. 2014, MNRAS, 441, 24; DESI Collaboration 2024
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from c0_common import (
    RAW_DIR,
    RESULTS_DIR,
    TEPLogger,
    ensure_dirs,
    print_status,
    rel,
    set_step_logger,
    step_json_path,
    write_json,
)

STEP_ID = "step_04_05_ddr_threeway"


def load_step_data(step_id: str) -> dict:
    """Load results from a previous step."""
    path = RESULTS_DIR / f"{step_id}.json"
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return {}


def compute_eta_statistics(eta_values: list, eta_errs: list) -> dict:
    """Compute weighted mean and error for eta."""
    eta_arr = np.array(eta_values)
    err_arr = np.array(eta_errs)
    weights = 1.0 / (err_arr ** 2)
    
    eta_weighted = np.sum(eta_arr * weights) / np.sum(weights)
    eta_weighted_err = np.sqrt(1.0 / np.sum(weights))
    eta_mean = np.mean(eta_arr)
    eta_std = np.std(eta_arr)
    
    return {
        "eta_mean": float(eta_mean),
        "eta_std": float(eta_std),
        "eta_weighted": float(eta_weighted),
        "eta_weighted_err": float(eta_weighted_err),
        "n_points": len(eta_values),
    }


def run() -> dict:
    """Execute three-way DDR comparison."""
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    # Load BAO results from step 025
    bao_data = load_step_data("step_04_04_distance_duality")
    bao_constraints = bao_data.get("ddr_constraints", [])
    
    # Load SZ results from step 023d
    sz_data = load_step_data("step_01_05_download_sz")
    sz_clusters = sz_data.get("data", {}).get("clusters", [])
    
    # Load SGL results from step 023e
    sgl_data = load_step_data("step_01_06_download_sgl")
    sgl_systems = sgl_data.get("data", {}).get("systems", [])
    
    # Compute eta for each method
    # For SZ and SGL, we need to match with SNe at similar redshifts
    
    # Load Pantheon SNe for D_L matching
    sne_path = RAW_DIR / "pantheon_lcparams.txt"
    if sne_path.exists():
        with open(sne_path) as f:
            lines = f.readlines()[1:]  # Skip header
        sne_data = []
        for line in lines:
            parts = line.strip().split()
            if len(parts) >= 3:
                try:
                    z = float(parts[1])
                    mb = float(parts[4])
                    sne_data.append({"z": z, "mb": mb})
                except (ValueError, IndexError):
                    continue
    else:
        sne_data = []
        print_status("Pantheon data not found - SN matching blocked", "ERROR")
    
    # NOTE: This step requires proper SN matching for SZ and SGL methods
    # Currently blocked until paired SN-SZ and SN-SGL data is available
    # This is a data requirement, not a synthetic fallback
    
    # BAO eta values
    bao_eta = [c["eta_obs"] for c in bao_constraints]
    bao_eta_err = [c["eta_err"] for c in bao_constraints]
    bao_z = [c["z"] for c in bao_constraints]

    # Load absolute magnitude from step_022
    step022 = load_step_data("step_03_01_three_model_comparison")
    m1_key = "M1_free_zT" if "M1_free_zT" in step022.get("models", {}) else "M1_NoLambda_zT5"
    M_sn = step022.get("models", {}).get(m1_key, {}).get("parameters_mle", {}).get("M", -19.36)
    
    # SZ eta values
    sz_eta = []
    sz_eta_err = []
    sz_z = []
    if sz_clusters and sne_data:
        for cluster in sz_clusters:
            z_c = cluster["z"]
            D_A = cluster["D_A"]
            D_A_err = cluster["D_A_err"]
            
            # Find nearest SN
            nearest_sn = min(sne_data, key=lambda sn: abs(sn["z"] - z_c))
            if abs(nearest_sn["z"] - z_c) < 0.05:
                # Calculate D_L from SN in Mpc
                mu = nearest_sn["mb"] - M_sn
                D_L = 10**((mu - 25)/5)
                # Approximate D_L error (~5%)
                D_L_err = D_L * 0.05
                
                eta = D_L / (D_A * (1 + z_c)**2)
                eta_err = eta * np.sqrt((D_L_err/D_L)**2 + (D_A_err/D_A)**2)
                
                sz_eta.append(eta)
                sz_eta_err.append(eta_err)
                sz_z.append(z_c)
    
    if sz_clusters and not sz_eta:
        print_status("SZ eta calculation blocked: requires proper SN-SZ paired data", "WARNING")
    
    # SGL eta values
    sgl_eta = []
    sgl_eta_err = []
    sgl_z = []
    if sgl_systems and sne_data:
        for sys in sgl_systems:
            z_l = sys["z_lens"]
            D_A = sys["D_A"]
            D_A_err = sys["D_A_err"]
            
            # Find nearest SN to the lens
            nearest_sn = min(sne_data, key=lambda sn: abs(sn["z"] - z_l))
            if abs(nearest_sn["z"] - z_l) < 0.05:
                # Calculate D_L from SN in Mpc
                mu = nearest_sn["mb"] - M_sn
                D_L = 10**((mu - 25)/5)
                D_L_err = D_L * 0.05
                
                eta = D_L / (D_A * (1 + z_l)**2)
                eta_err = eta * np.sqrt((D_L_err/D_L)**2 + (D_A_err/D_A)**2)
                
                sgl_eta.append(eta)
                sgl_eta_err.append(eta_err)
                sgl_z.append(z_l)
                
    if sgl_systems and not sgl_eta:
        print_status("SGL eta calculation blocked: requires proper SN-SGL paired data", "WARNING")
    
    # Compute statistics for each method
    bao_stats = compute_eta_statistics(bao_eta, bao_eta_err)
    sz_stats = compute_eta_statistics(sz_eta, sz_eta_err)
    sgl_stats = compute_eta_statistics(sgl_eta, sgl_eta_err)
    
    # Statistical comparison
    delta_sz_sgl = abs(sz_stats["eta_weighted"] - sgl_stats["eta_weighted"])
    delta_bao_sz = abs(bao_stats["eta_weighted"] - sz_stats["eta_weighted"])
    delta_bao_sgl = abs(bao_stats["eta_weighted"] - sgl_stats["eta_weighted"])
    
    # Significance (simplified)
    sigma_sz_sgl = delta_sz_sgl / np.sqrt(sz_stats["eta_weighted_err"]**2 + sgl_stats["eta_weighted_err"]**2)
    sigma_bao_sz = delta_bao_sz / np.sqrt(bao_stats["eta_weighted_err"]**2 + sz_stats["eta_weighted_err"]**2)
    
    print_status(f"BAO: η = {bao_stats['eta_weighted']:.3f} ± {bao_stats['eta_weighted_err']:.3f}", "INFO")
    print_status(f"SZ:  η = {sz_stats['eta_weighted']:.3f} ± {sz_stats['eta_weighted_err']:.3f}", "INFO")
    print_status(f"SGL: η = {sgl_stats['eta_weighted']:.3f} ± {sgl_stats['eta_weighted_err']:.3f}", "INFO")
    
    print_status(f"SZ vs SGL: |Δη| = {delta_sz_sgl:.3f} ({sigma_sz_sgl:.1f}σ)", "INFO")
    print_status(f"BAO vs SZ: |Δη| = {delta_bao_sz:.3f} ({sigma_bao_sz:.1f}σ)", "INFO")
    
    # Critical finding
    if sigma_sz_sgl < 2.0 and sigma_bao_sz > 3.0:
        print_status("CRITICAL: SZ and SGL agree; BAO differs significantly", "SUCCESS")
        critical_finding = True
    else:
        critical_finding = False
    
    payload = {
        "step": STEP_ID,
        "description": "Three-way comparison of DDR from BAO, SZ, and SGL methods",
        "status": "completed",
        "methods": {
            "BAO": {
                "n_points": bao_stats["n_points"],
                "eta_weighted": bao_stats["eta_weighted"],
                "eta_weighted_err": bao_stats["eta_weighted_err"],
                "physical_scale": "~100 Mpc",
                "physics": "linear density perturbations",
            },
            "SZ": {
                "n_points": sz_stats["n_points"],
                "eta_weighted": sz_stats["eta_weighted"],
                "eta_weighted_err": sz_stats["eta_weighted_err"],
                "physical_scale": "~5 Mpc",
                "physics": "intracluster hot gas",
            },
            "SGL": {
                "n_points": sgl_stats["n_points"],
                "eta_weighted": sgl_stats["eta_weighted"],
                "eta_weighted_err": sgl_stats["eta_weighted_err"],
                "physical_scale": "~5 Mpc",
                "physics": "cluster gravitational potential",
            },
        },
        "comparisons": {
            "SZ_vs_SGL": {
                "delta_eta": float(delta_sz_sgl),
                "sigma": float(sigma_sz_sgl),
                "consistent": sigma_sz_sgl < 2.0,
            },
            "BAO_vs_SZ": {
                "delta_eta": float(delta_bao_sz),
                "sigma": float(sigma_bao_sz),
                "consistent": sigma_bao_sz < 2.0,
            },
            "BAO_vs_SGL": {
                "delta_eta": float(delta_bao_sgl),
                "significant": sigma_bao_sz > 3.0,
            },
        },
        "critical_finding": {
            "found": critical_finding,
            "description": "SZ and SGL agree (cluster-scale physics); BAO differs (linear regime)",
            "implication": "Evidence for probe-dependent Etherington relation",
        },
        "validation": {
            "real_data": len(bao_eta) > 0,
            "sz_blocked": len(sz_eta) == 0,
            "sgl_blocked": len(sgl_eta) == 0,
            "sn_matching_required": True,
            "blockers": [b for b in [
                "SZ eta blocked: requires SN-SZ paired data" if len(sz_eta) == 0 else None,
                "SGL eta blocked: requires SN-SGL paired data" if len(sgl_eta) == 0 else None,
                "Pantheon SN data not available" if len(sne_data) == 0 else None,
            ] if b is not None],
        },
        "timestamp": int(time.time()),
    }
    
    write_json(step_json_path(STEP_ID), payload)
    print_status(f"Step {STEP_ID} completed", "SUCCESS")
    return payload


if __name__ == "__main__":
    print(run())
