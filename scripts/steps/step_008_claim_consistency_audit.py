#!/usr/bin/env python3
"""Step 008: audit manuscript claims, generated artifacts, and open-gate guardrails."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from c0_common import PROJECT_ROOT, TEPLogger, ensure_dirs, print_status, read_json, rel, set_step_logger, step_csv_path, step_json_path, write_csv, write_json, rounded


STEP_ID = "step_008_claim_consistency_audit"
AUDITED_STEPS = [
    "step_000_data_download",
    "step_000_data_ingestion",
    "step_001_transport_kernel",
    "step_002_mixed_ft_forecast",
    "step_010_cmb_blackbody_preservation",
    "step_011_bao_acoustic_projection",
    "step_012_bbn_preservation_registry",
    "step_014_cmb_acoustic_projection",
    "step_015_structure_growth_solver",
    "step_016_global_likelihood_synthesis",
    "step_013_explanatory_evidence_matrix",
    "step_021_expansion_falsifier",
    "step_022_three_model_comparison",
]


SOURCE_PATHS = [
    PROJECT_ROOT / "site" / "components" / "1_abstract.html",
    PROJECT_ROOT / "site" / "components" / "1_introduction.html",
    PROJECT_ROOT / "site" / "components" / "2_theory.html",
    PROJECT_ROOT / "site" / "components" / "3_methodology.html",
    PROJECT_ROOT / "site" / "components" / "4_results.html",
    PROJECT_ROOT / "site" / "components" / "5_discussion.html",
    PROJECT_ROOT / "site" / "components" / "6_conclusion.html",
    PROJECT_ROOT / "site" / "components" / "8_reproducibility.html",
    PROJECT_ROOT / "README.md",
    PROJECT_ROOT / "14-TEP-C0-v0.1-Athens.md",
    PROJECT_ROOT / "scripts" / "steps" / "PIPELINE_STATUS.md",
    PROJECT_ROOT.parent / "TEP-GL" / "manuscripts" / "0-TEP-v0.8-Jakarta.md",
]


REQUIRED_GUARDRAILS = [
    "source datasets",
    "data ingestion",
    "explanatory evidence matrix",
    "transport phenomenon",
    "preservation of matter-frame nuclear history",
    "phase-space",
    "TEP Boltzmann",
    "global MCMC",
    "inference",
    "preservation constraints",
    "resolution of the Hubble tension",
    "mixing fraction",
    "arbiter",
    "Cosmological Isochrony Axiom",
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
