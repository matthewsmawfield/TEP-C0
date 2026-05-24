#!/usr/bin/env python3
"""Step 009: evidence gate summary and global synthesis.

References: consumes cited Pantheon+, FIRAS/Planck, BAO, and BBN step artifacts.
"""

from __future__ import annotations

import pandas as pd
from pathlib import Path
from c0_common import (
    TEPLogger, ensure_dirs, print_status, read_json, rel, 
    rounded, set_step_logger, step_csv_path, step_json_path, 
    write_csv, write_json, RESULTS_DIR
)

STEP_ID = "step_009_evidence_gate_summary"

def run() -> dict:
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    print_status("Compiling summary from all pipeline steps", "PROCESS")
    
    # Load core artifacts using project-standard utility
    try:
        transport = read_json(step_json_path("step_001_transport_kernel"))
        global_like = read_json(step_json_path("step_016_global_likelihood_synthesis"))
        audit = read_json(step_json_path("step_008_claim_consistency_audit"))
        ingestion = read_json(step_json_path("step_000_data_ingestion"))
        mcmc = read_json(step_json_path("step_018_tep_mcmc_inference"))
        comparison = read_json(step_json_path("step_022_three_model_comparison"))
    except FileNotFoundError as e:
        print_status(f"Missing dependency: {e}", "ERROR")
        raise

    # Evidence Gates Consolidation
    comparison_validation = comparison.get("validation", {})
    research_grade_open = comparison_validation.get("model_comparison_gate") == "open"
    m1_key = "M1_NoLambda_zT1" if "M1_NoLambda_zT1" in comparison.get("models", {}) else "M1_NoLambda_zT5"
    m1_payload = comparison["models"][m1_key]
    epsilon_t = next(
        (row["median"] for row in mcmc["results"] if row["parameter"] in {"epsilon_T", "ft", "f_T"}),
        None,
    )
    gates = [
        {
            "gate": "SN_Pantheon_Model_Comparison",
            "status": "pass" if research_grade_open else "blocked",
            "metric": m1_payload.get("log_likelihood_mle", m1_payload.get("log_likelihood")),
        },
        {
            "gate": "Hubble_Tension_Diagnostic",
            "status": "diagnostic",
            "metric": global_like["metrics"]["tension_resolution_factor"],
        },
        {
            "gate": "Temporal_Shear_Posterior_Diagnostic",
            "status": "diagnostic",
            "metric": epsilon_t,
        },
    ]
    write_csv(RESULTS_DIR / "tep_c0_evidence_gates.csv", gates)
    print_status(f"Consolidated {len(gates)} evidence gates", "SUCCESS")

    # Global Synthesis Payload
    payload = {
        "step": STEP_ID,
        "summary": {
            "h0_mcmc": global_like["metrics"]["mcmc_h0"],
            "f_t_mcmc": global_like["metrics"]["mcmc_ft"],
            "trf": global_like["metrics"]["tension_resolution_factor"],
            "audit_status": audit["metrics"]["audit_status"],
            "pass_rate": audit["metrics"]["pass_rate"],
            "research_grade_claim_gate": "open" if research_grade_open else "blocked",
            "claim_blockers": comparison_validation.get("blockers", []),
        },
        "reconstruction_metrics": {
            "transport_recovery_error": transport["metrics"]["max_reconstruction_error_z"],
            "transport_low_z_ht": transport["metrics"]["low_z_recovered_HT_km_s_Mpc"],
        },
        "interpretation": "TEP-C0 records the current empirical diagnostics and exposes claim gates explicitly. Expansion-replacement, CMB-replacement, and BBN-replacement claims remain blocked unless the research-grade likelihood, CMB, BBN, and robustness gates are open."
    }
    write_json(step_json_path(STEP_ID), payload)
    
    print_status(f"Final TEP-C0 Summary generated", "SUCCESS")
    return payload

if __name__ == "__main__":
    run()
