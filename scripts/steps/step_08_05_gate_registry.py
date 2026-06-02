#!/usr/bin/env python3
"""Step 007: Evidence Gate Registry - aggregate all test results.

Collects and summarizes all pipeline test results into unified evidence matrix.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
from c0_common import TEPLogger, ensure_dirs, print_status, read_json, rounded, set_step_logger, step_json_path, write_json

STEP_ID = "step_08_05_gate_registry"


def load_step_result(step_name: str) -> dict:
    """Load result from a previous step."""
    result_file = Path(f"results/{step_name}.json")
    if result_file.exists():
        return read_json(result_file)
    return None


def run():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    # Collect results from all major steps
    steps_to_check = [
        'step_03_01_three_model_comparison',
        'step_04_01_sn_time_dilation',
        'step_04_03_tolman_sb',
        'step_04_04_distance_duality',
        'step_07_02_redshift_drift',
        'step_04_08_environment_residuals'
    ]
    
    evidence_summary = {}
    all_passed = True
    
    for step in steps_to_check:
        result = load_step_result(step)
        
        # Check if the step genuinely passed its validation checks
        is_completed = False
        if result:
            val = result.get('validation', {})
            # Look for strong validation signals
            if val.get('research_grade', False) or val.get('all_validated', False) or result.get('test_passed', False):
                is_completed = True
            elif result.get('status') in ['completed', 'validated', 'passed']:
                is_completed = True
                
        if is_completed:
            evidence_summary[step] = {
                'completed': True,
                'key_finding': result.get('key_finding', 'N/A'),
                'interpretation': result.get('interpretation', 'N/A')[:100] + '...'
            }
        else:
            evidence_summary[step] = {'completed': False, 'reason': 'Failed validation or missing'}
            all_passed = False
    
    results = {
        'step': STEP_ID,
        'description': 'Evidence gate registry aggregating all test results',
        'evidence_summary': evidence_summary,
        'all_critical_tests_completed': all_passed,
        'n_tests_completed': sum(1 for s in evidence_summary.values() if s.get('completed', False)),
        'n_tests_total': len(steps_to_check)
    }
    
    write_json(step_json_path(STEP_ID), results)
    print_status(f"Step {STEP_ID} completed: {results['n_tests_completed']}/{results['n_tests_total']} tests", "SUCCESS")
    return results


if __name__ == "__main__":
    run()
