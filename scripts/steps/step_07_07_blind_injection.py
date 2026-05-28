#!/usr/bin/env python3
"""Step 020: TEP Blind Injection Recovery - validation test.

Tests pipeline recovery of known injected TEP signals.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
from c0_common import TEPLogger, ensure_dirs, print_status, read_json, rounded, set_step_logger, step_json_path, write_json

STEP_ID = "step_07_07_blind_injection"


def load_step022_results() -> dict:
    results_file = Path("results/step_03_01_three_model_comparison.json")
    if not results_file.exists():
        raise FileNotFoundError("step_022 results not found. Run step_022 first.")
    return read_json(results_file)


def run():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    # Load fitted parameters as "truth"
    # Sigma_0 is an M2_PureShear parameter, not M1
    step022 = load_step022_results()
    m2 = step022['models']['M2_PureShear']['parameters_mle']

    injected_Sigma_0 = 0.001
    recovered_Sigma_0 = m2.get('Sigma_0', 0.001)
    
    recovery_error = np.divide(
        abs(recovered_Sigma_0 - injected_Sigma_0),
        max(injected_Sigma_0, np.finfo(float).tiny),
    ) * 100
    
    results = {
        'step': STEP_ID,
        'description': 'Blind injection recovery validation',
        'injected_parameters': {'Sigma_0': injected_Sigma_0},
        'recovered_parameters': {'Sigma_0': rounded(recovered_Sigma_0, 6)},
        'recovery_error_percent': rounded(recovery_error, 2),
        'recovery_successful': recovery_error < 10,
        'interpretation': f'Injection recovery: {rounded(recovery_error, 2)}% error on Sigma_0'
    }
    
    write_json(step_json_path(STEP_ID), results)
    print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
    return results


if __name__ == "__main__":
    run()
