#!/usr/bin/env python3
"""Step 036: JWST Early Universe Test (z > 6).

Plan for testing TEP at z > 6 using JWST observations:
- High-z galaxy luminosity functions
- UV surface brightness
- Angular size measurements

TEP Predictions at z > 6:
-------------------------
- η(z) → 1 (unscreened regime)
- S(z) → 1 (no temporal coupling)
- Standard distance duality relation restored
- Growth of structure converges to ΛCDM

This is a planning step - actual data to be acquired from JWST surveys.

Tier: PLANNING / RESEARCH

References:
  - JWST CEERS survey (Finkelstein et al. 2022)
  - JWST JADES survey (Robertson et al. 2023)
  - ASTRODEEP-JWST catalogue (Merlin et al. 2024)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

sys.path.insert(0, str(Path(__file__).parent))
from c0_common import TEPLogger, set_step_logger, ensure_dirs, print_status, step_json_path, write_json

STEP_ID = "step_07_03_jwst_test"


def plan_jwst_tests() -> Dict[str, Any]:
    """Create plan for JWST TEP tests."""
    return {
        "test_1_angular_sizes": {
            "description": "Galaxy angular sizes at z > 6",
            "tep_prediction": "Standard angular diameter distance",
            "lcdm_prediction": "Standard angular diameter distance",
            "discriminating": "No - both predict same at high z",
            "data_source": "JWST CEERS, JADES surveys",
            "priority": "Medium",
        },
        "test_2_surface_brightness": {
            "description": "UV surface brightness evolution",
            "tep_prediction": "Convergence to standard (1+z)^-4 dimming",
            "lcdm_prediction": "(1+z)^-4 Tolman dimming",
            "discriminating": "Yes - TEP predicts convergence",
            "data_source": "ASTRODEEP-JWST catalogue",
            "priority": "High",
        },
        "test_3_luminosity_functions": {
            "description": "Galaxy luminosity functions",
            "tep_prediction": "Standard evolution at high z",
            "lcdm_prediction": "Standard evolution",
            "discriminating": "No - both agree at high z",
            "data_source": "JWST deep fields",
            "priority": "Low",
        },
    }


def run() -> dict:
    """Run JWST planning step."""
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()
    
    tests = plan_jwst_tests()
    
    print_status("Planned JWST tests:", "INFO")
    for test_name, test_info in tests.items():
        print_status(f"  {test_name}: {test_info['priority']} priority", "INFO")
        print_status(f"    {test_info['description']}", "INFO")
    
    payload = {
        "step": STEP_ID,
        "status": "planning",
        "planned_tests": tests,
        "data_acquisition": {
            "astodeep_jwst": "Available (6,860 galaxies)",
            "ceers": "To be downloaded",
            "jades": "To be downloaded",
        },
    }
    
    write_json(step_json_path(STEP_ID), payload)
    print_status(f"Step {STEP_ID} completed (planning)", "SUCCESS")
    return payload


if __name__ == "__main__":
    run()
