#!/usr/bin/env python3
"""Step 013: Explanatory Evidence Matrix - Strictly Empirical Summary.

Note: This step does not perform FLRW calculations directly. It compiles results
from other steps that use CosmologyFLRW (which includes radiation component Or0)
and TEPModulatedCosmology (which inherits radiation from CosmologyFLRW).
All cosmological calculations use explicit curvature parameters (Ok0=0.0 for flat universe).
The FLRW framework is used in the underlying cosmology modules with proper radiation handling.
"""

from __future__ import annotations

import pandas as pd
from pathlib import Path
from c0_common import (
    TEPLogger, ensure_dirs, print_status, read_json, rel, 
    rounded, set_step_logger, step_csv_path, step_json_path, 
    write_csv, write_json
)

STEP_ID = "step_08_04_evidence_matrix"

def run() -> dict:
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    # Load real results
    comp = read_json(step_json_path("step_03_01_three_model_comparison"))
    m1_key = "M1_NoLambda_zT1" if "M1_NoLambda_zT1" in comp.get("models", {}) else "M1_NoLambda_zT5"
    m1 = comp["models"][m1_key]
    m0 = comp["models"]["M0a_LCDM"]
    comparison_open = comp.get("validation", {}).get("model_comparison_gate") == "open"
    
    transport = read_json(step_json_path("step_02_01_transport_kernel"))
    bbn = read_json(step_json_path("step_05_06_bbn_registry"))
    
    # Handle optional dependencies gracefully
    cmb_peak_path = step_json_path("step_05_08_cmb_acoustic")
    cmb_peak = read_json(cmb_peak_path) if cmb_peak_path.exists() else {}
    
    falsifier_path = step_json_path("step_08_01_expansion_falsifier")
    falsifier = read_json(falsifier_path) if falsifier_path.exists() else {}

    matrix = [
        {
            "observable": "Cosmic Redshift (z)",
            "tep_explanation": "Projection of Temporal Shear gradients along the photon path.",
            "numerical_anchor": f"FLRW error {transport['metrics']['max_reconstruction_error_z']:.2e}",
            "status": "Verified as a transport identity",
        },
        {
            "observable": "SN 1a Distance-Redshift",
            "tep_explanation": "TEP matter model incorporating matter-frame effects.",
            "numerical_anchor": f"Dimensionless distance comparison (H0_ref = 70.0 km/s/Mpc); M_TEP = {m1['parameters_mle'].get('M', 'N/A')} vs M_LCDM = {m0['parameters_mle'].get('M', 'N/A')}",
            "status": "diagnostic" if not comparison_open else "validated",
        },
        {
            "observable": "Hubble Constant agreement",
            "tep_explanation": "Resolution via temporal shear effects.",
            "numerical_anchor": f"Dimensionless distances fix H0_ref = 70.0 km/s/Mpc; no H0 parameter in fit.",
            "status": "diagnostic",
        },
        {
            "observable": "Distance Duality (Xi_T)",
            "tep_explanation": "Measurable residuals in high-z luminosity/angular diameter relations.",
            "numerical_anchor": f"Tolman Δχ²: {str(falsifier.get('tolman_delta_chi2', 'N/A'))[:30]}...",
            "status": falsifier.get('confidence', 'Preliminary'),
        },
        {
            "observable": "CMB Acoustic Scales",
            "tep_explanation": "Preserved matter-frame acoustic horizon properly projected through shear connection.",
            "numerical_anchor": f"Consistency: {cmb_peak.get('consistency', 'N/A')}",
            "status": "blocked pending full CMB likelihood",
        }
    ]

    write_csv(step_csv_path(STEP_ID), matrix)
    
    # Save payload
    payload = {
        "step": STEP_ID,
        "matrix": matrix,
        "artifacts": {
            "csv": rel(step_csv_path(STEP_ID))
        },
        "validation": {
            "claim_gate": "open" if comparison_open else "blocked",
            "blockers": comp.get("validation", {}).get("blockers", []),
        },
    }
    write_json(step_json_path(STEP_ID), payload)
    print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
    return payload

if __name__ == "__main__":
    run()
