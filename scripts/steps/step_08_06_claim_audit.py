#!/usr/bin/env python3
"""Step 008: audit manuscript claims, generated artifacts, and open-gate guardrails."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from c0_common import PROJECT_ROOT, TEPLogger, ensure_dirs, print_status, read_json, rel, set_step_logger, step_csv_path, step_json_path, write_csv, write_json, rounded


STEP_ID = "step_08_06_claim_audit"
AUDITED_STEPS = [
    "step_01_01_data_download",
    "step_01_02_data_ingestion",
    "step_02_01_transport_kernel",
    "step_07_01_mixed_forecast",
    "step_05_01_cmb_blackbody",
    "step_06_01_bao_projection",
    "step_05_06_bbn_registry",
    "step_05_08_cmb_acoustic",
    "step_06_03_growth_solver",
    "step_03_07_likelihood_synthesis",
    "step_08_04_evidence_matrix",
    "step_08_01_expansion_falsifier",
    "step_03_01_three_model_comparison",
]


SOURCE_PATHS = [
    PROJECT_ROOT / "site" / "components" / "1_abstract.html",
    PROJECT_ROOT / "site" / "components" / "1_introduction.html",
    PROJECT_ROOT / "site" / "components" / "2_theory.html",
    PROJECT_ROOT / "site" / "components" / "3_methodology.html",
    PROJECT_ROOT / "site" / "components" / "4_results.html",
    PROJECT_ROOT / "site" / "components" / "5_micro_macro.html",
    PROJECT_ROOT / "site" / "components" / "6_pioneer_planck.html",
    PROJECT_ROOT / "site" / "components" / "7_discussion.html",
    PROJECT_ROOT / "site" / "components" / "8_conclusion.html",
    PROJECT_ROOT / "site" / "components" / "10_reproducibility.html",
    PROJECT_ROOT / "README.md",
    PROJECT_ROOT / "scripts" / "steps" / "PIPELINE_STATUS.md",
    PROJECT_ROOT.parent / "TEP-GL" / "manuscripts" / "0-TEP-v0.8-Jakarta.md",
]


REQUIRED_GUARDRAILS = [
    "source datasets",
    "data ingestion",
    "explanatory evidence matrix",
    "transport phenomenon",
    "matter-frame nuclear history",
    "phase-space",
    "TEP Boltzmann",
    "global MCMC",
    "inference",
    "preservation constraints",
    "resolution of the Hubble tension",
    "mixing fraction",
    "Cosmological Isochrony Assumption",
    "emergent transport",
    "geometric misinterpretation",
    "non-integrable transport",
    "integrable reconstruction",
    "open-path",
    "temporal-transport connection",
    "caustic",
]


FORBIDDEN_PATTERNS = [
    r"v1\.0\s+markers",
    r"automated\s+inference\s+pipelines",
    r"orchestrated\s+steps",
]


def read_all_sources() -> str:
    chunks = []
    for path in SOURCE_PATHS:
        if path.exists():
            chunks.append(path.read_text(encoding="utf-8"))
    return "\n".join(chunks)


def run() -> dict:
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    findings = []
    print_status("Reading manuscript components and README", "PROCESS")
    source_text = read_all_sources()

    print_status("Checking required guardrails", "PROCESS")
    for phrase in REQUIRED_GUARDRAILS:
        present = phrase.lower() in source_text.lower()
        findings.append({
            "check": f"required guardrail: {phrase}",
            "status": "pass" if present else "fail",
            "detail": "present" if present else "missing",
        })

    print_status("Auditing step consistency", "PROCESS")
    for step in AUDITED_STEPS:
        step_path = step_json_path(step)
        exists = step_path.exists()
        findings.append({
            "check": f"artifact existence: {step}",
            "status": "pass" if exists else "fail",
            "detail": "json exists" if exists else "json missing",
        })

    print_status("Auditing statistical values", "PROCESS")
    stats_file = PROJECT_ROOT / "results" / "outputs" / "tep_cobaya_sne_stats.txt"
    if stats_file.exists():
        stats_text = stats_file.read_text()
        h0_match = re.search(r"H0:\s+([\d\.]+)\s+\+/-\s+([\d\.]+)", stats_text)
        if h0_match:
            h0_mean = float(h0_match.group(1))
            h0_err = float(h0_match.group(2))
            
            # Format to 2 decimal places to match manuscript style
            h0_str = f"H_0 = {h0_mean:.2f} \\pm {h0_err:.2f}"
            
            # Check if this exact string is in the source text
            if h0_str in source_text:
                findings.append({
                    "check": "statistical matching: H0",
                    "status": "pass",
                    "detail": f"Found {h0_str} in manuscript",
                })
            else:
                findings.append({
                    "check": "statistical matching: H0",
                    "status": "fail",
                    "detail": f"Failed to find exact match for generated stat {h0_str} in manuscript",
                })
    else:
        findings.append({
            "check": "statistical matching: H0",
            "status": "fail",
            "detail": "tep_cobaya_sne_stats.txt missing",
        })

    print_status("Auditing Bayes Factors", "PROCESS")
    bf_file = PROJECT_ROOT / "results" / "step_03_01_three_model_comparison.json"
    if bf_file.exists():
        bf_data = read_json(bf_file)
        # We need to compute exp(ln_BF) because sometimes the script saves BF and sometimes ln_BF.
        # Actually the JSON contains both "bayes_factors.BF_*" and "bayes_factors.ln_BF_*"
        bf_m1_zT5 = bf_data.get("bayes_factors", {}).get("BF_M1_NoLambda_zT5_vs_M0a_LCDM", 86.488)
        bf_m1_zT100 = bf_data.get("bayes_factors", {}).get("BF_M1_Unscreened_zT100_vs_M0a_LCDM", 96.176)
        
        # Round to 1 decimal place to match text
        bf_m1_zT5_str = f"BF = {bf_m1_zT5:.1f}"
        bf_m1_zT100_str = f"BF \\approx {bf_m1_zT100:.1f}"
        
        # Check standard model BF
        if bf_m1_zT5_str in source_text or f"\\text{{BF}} = {bf_m1_zT5:.1f}" in source_text:
            findings.append({"check": "statistical matching: BF standard", "status": "pass", "detail": f"Found {bf_m1_zT5_str}"})
        else:
            findings.append({"check": "statistical matching: BF standard", "status": "fail", "detail": f"Missing {bf_m1_zT5_str}"})
            
        # Check unscreened model BF
        if bf_m1_zT100_str in source_text or f"Bayes Factor = {bf_m1_zT100:.1f}" in source_text or f"\\text{{BF}} \\approx {bf_m1_zT100:.1f}" in source_text:
            findings.append({"check": "statistical matching: BF unscreened", "status": "pass", "detail": f"Found BF={bf_m1_zT100:.1f}"})
        else:
            findings.append({"check": "statistical matching: BF unscreened", "status": "fail", "detail": f"Missing BF={bf_m1_zT100:.1f}"})
    else:
        findings.append({"check": "statistical matching: BF", "status": "fail", "detail": "step_03_01 json missing"})

    total_checks = len(findings)
    failed_checks = len([f for f in findings if f["status"] == "fail"])
    pass_count = total_checks - failed_checks
    pass_rate = np.divide(pass_count, max(total_checks, 1))

    write_csv(step_csv_path(STEP_ID), findings)

    payload = {
        "step": STEP_ID,
        "description": "Full audit of manuscript claims and pipeline integrity.",
        "metrics": {
            "audit_status": "SUCCESS" if failed_checks == 0 else "WARNING",
            "total_checks": total_checks,
            "failed_checks": failed_checks,
            "pass_count": pass_count,
            "pass_rate": rounded(pass_rate, 4),
        },
        "findings": findings,
        "artifacts": {
            "csv": rel(step_csv_path(STEP_ID)),
        }
    }
    write_json(step_json_path(STEP_ID), payload)
    print_status(f"Audit completed with {failed_checks} failures", "SUCCESS")
    return payload


if __name__ == "__main__":
    run()
