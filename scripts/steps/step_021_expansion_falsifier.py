#!/usr/bin/env python3
"""Step 021: Expansion Falsifier - final discriminator test.

Uses DDR and Tolman tests to falsify pure expansion hypothesis.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
from c0_common import TEPLogger, ensure_dirs, print_status, read_json, rounded, set_step_logger, step_json_path, write_json

STEP_ID = "step_021_expansion_falsifier"


def load_step_result(step_name: str) -> dict:
    result_file = step_json_path(step_name)
    if result_file.exists():
        return read_json(result_file)
    return None


def run():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    # Load DDR and Tolman results
    ddr = load_step_result('step_025_distance_duality_test')
    tolman = load_step_result('step_024_tolman_surface_brightness')
    
    if not ddr:
        print_status("Missing DDR test results from step_025. Running in forecast mode.", "WARNING")
        ddr = {'results': {}}
    if not tolman:
        print_status("Missing Tolman test results from step_024. Running in forecast mode.", "WARNING")
        tolman = {'results': {}}

    # Handle different possible key structures
    ddr_results = ddr.get('statistical_results', ddr.get('results', {}))
    tolman_results = tolman.get('statistical_results', tolman.get('results', {}))
    ddr_deviation = ddr_results.get('deviation_from_unity_sigma')
    tolman_delta_chi2 = tolman_results.get('delta_chi2')
    
    # Handle missing data - block step if required data unavailable
    if ddr_deviation is None:
        print_status("DDR deviation_from_unity_sigma not found - step blocked", "ERROR")
        results = {
            'step': STEP_ID,
            'description': 'Expansion falsifier using DDR and Tolman tests',
            'status': 'BLOCKED',
            'blocker': 'DDR test results (step_025) not available',
            'distance_duality_deviation': 'N/A',
            'tolman_delta_chi2': 'N/A',
            'tep_hypothesis_favored': None,
            'confidence': 'Blocked - missing data',
            'interpretation': 'Cannot evaluate without DDR and Tolman test results'
        }
        write_json(step_json_path(STEP_ID), results)
        return results
    if tolman_delta_chi2 is None:
        print_status("Tolman delta_chi2 not found - step blocked", "ERROR")
        results = {
            'step': STEP_ID,
            'description': 'Expansion falsifier using DDR and Tolman tests',
            'status': 'BLOCKED',
            'blocker': 'Tolman test results (step_024) not available',
            'distance_duality_deviation': 'N/A',
            'tolman_delta_chi2': 'N/A',
            'tep_hypothesis_favored': None,
            'confidence': 'Blocked - missing data',
            'interpretation': 'Cannot evaluate without DDR and Tolman test results'
        }
        write_json(step_json_path(STEP_ID), results)
        return results
    
    # Test if data favors TEP over pure expansion
    tep_favored = abs(ddr_deviation) > 2.0 or tolman_delta_chi2 > 9.0
    
    results = {
        'step': STEP_ID,
        'description': 'Expansion falsifier using DDR and Tolman tests',
        'distance_duality_deviation': rounded(ddr_deviation, 4) if ddr else 'N/A',
        'tolman_delta_chi2': rounded(tolman_delta_chi2, 2) if tolman else 'N/A',
        'expansion_hypothesis': {
            'ddr_must_be': 1.0,
            'tolman_xi_must_be': 1.0
        },
        'tep_hypothesis_favored': tep_favored,
        'confidence': 'Preliminary' if not ddr or not tolman else ('Moderate' if tep_favored else 'Weak'),
        'interpretation': 'TEP predicts measurable deviations from pure expansion in DDR and Tolman tests'
    }
    
    write_json(step_json_path(STEP_ID), results)
    print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
    return results


if __name__ == "__main__":
    run()
