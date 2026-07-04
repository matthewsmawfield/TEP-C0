#!/usr/bin/env python3
"""Step 028: CMB full-spectra resolver gate.

Delegates to step 017 but writes its own validation payload so the pipeline can
separate "CLASS diagnostic ran" from "research-grade CMB replacement passed".
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from c0_common import ensure_dirs, step_json_path, write_json
from step_05_03_cmb_boltzmann import run as run_step017

STEP_ID = "step_05_04_cmb_spectra"


def run() -> dict:
    ensure_dirs()
    source = run_step017()
    validation = source.get("validation", {})
    payload = {
        "step": STEP_ID,
        "status": "completed" if validation.get("research_grade_cmb") else "blocked",
        "source_step": "step_05_03_cmb_boltzmann",
        "cmb_results": source.get("cmb_results", {}),
        "validation": {
            "class_available": validation.get("class_available", False),
            "tep_class_available": validation.get("tep_class_available", False),
            "local_class_source_tep_patch_present": validation.get("local_class_source_tep_patch_present", False),
            "tep_zero_limit_ok": validation.get("tep_zero_limit_ok", False),
            "research_grade_cmb": validation.get("research_grade_cmb", False),
            "claim_gate": "open" if validation.get("research_grade_cmb") else "blocked",
            "blockers": validation.get("blockers", [
                "Step 017 did not report CMB validation status.",
            ]),
        },
    }
    write_json(step_json_path(STEP_ID), payload)
    return payload


if __name__ == "__main__":
    result = run()
    print({"step": result["step"], "validation": result["validation"]})
