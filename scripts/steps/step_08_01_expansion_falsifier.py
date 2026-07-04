#!/usr/bin/env python3
"""Step 021: Expansion Falsifier - final discriminator test.

Uses DDR and Tolman tests to falsify pure expansion hypothesis.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
from c0_common import TEPLogger, ensure_dirs, print_status, read_json, rounded, set_step_logger, step_json_path, write_json

STEP_ID = "step_08_01_expansion_falsifier"


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
    ddr = load_step_result('step_04_04_distance_duality')
    tolman = load_step_result('step_04_03_tolman_sb')
    
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
    
    # Neither DDR nor Tolman currently provides a clean TEP discriminator:
    # - DDR: both LCDM and TEP predict eta=1 by construction; the 6.6σ deviation
    #   reflects systematic tension in the compiled D_L/D_A sample (Planck D_L vs BAO D_A)
    # - Tolman: TEP predicts n_TEP ≈ 4.8, LCDM predicts n=4.0, data shows n=3.375;
    #   TEP is further from the data than LCDM, and evolution/K-correction systematics dominate
    tep_favored = False
    
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
        'confidence': 'Blocked',
        'interpretation': 'DDR and Tolman sectors are currently blocked as clean discriminators. DDR compilation mixes inconsistent D_L/D_A sources. Tolman is dominated by galaxy-evolution systematics. A self-consistent TEP-derived compilation and evolution-corrected Tolman analysis are required before these tests can discriminate.'
    }
    
    write_json(step_json_path(STEP_ID), results)
    print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
    return results


if __name__ == "__main__":
    run()
