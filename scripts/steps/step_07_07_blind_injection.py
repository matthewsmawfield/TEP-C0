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
    results_file = Path("results/observer_injection_recovery.json")
    if not results_file.exists():
        raise FileNotFoundError("Current signed recovery result not found. Run scripts/steps/refresh_observer_injection_recovery.py first.")
    return read_json(results_file)


def run():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    # Report the real injection-recovery test performed in step_03_01:
    # a Pantheon-realistic mock generated at epsilon_shear_los = +0.23 (-0.23)
    # was refit blind and recovered to ~1e-7 relative precision.
    step022 = load_step022_results()
    pos = step022.get('positive_injection_test', {})
    neg = step022.get('negative_injection_test', {})

    injected = float(pos.get('injected_epsilon_shear_los', 0.23))
    recovered = float(pos.get('recovered_epsilon_shear_los', np.nan))
    injected_neg = float(neg.get('injected_epsilon_shear_los', -0.23))
    recovered_neg = float(neg.get('recovered_epsilon_shear_los', np.nan))

    def rel_err(rec, inj):
        return float(np.divide(
            abs(rec - inj), max(abs(inj), np.finfo(float).tiny))) * 100

    recovery_error_pos = rel_err(recovered, injected)
    recovery_error_neg = rel_err(recovered_neg, injected_neg)
    recovery_error = recovery_error_pos

    results = {
        'step': STEP_ID,
        'description': 'Blind injection recovery validation (M1 line-of-sight shear)',
        'injected_parameters': {
            'epsilon_shear_los_positive': injected,
            'epsilon_shear_los_negative': neg.get('injected_epsilon_shear_los'),
        },
        'recovered_parameters': {
            'epsilon_shear_los_positive': rounded(recovered, 6),
            'epsilon_shear_los_negative': rounded(neg.get('recovered_epsilon_shear_los', np.nan), 6),
        },
        'recovery_error_percent': recovery_error,
        'recovery_error_percent_positive': recovery_error_pos,
        'recovery_error_percent_negative': recovery_error_neg,
        'recovery_successful': bool(pos.get('passed') and neg.get('passed')),
        'interpretation': (
            f'Injection recovery: {recovery_error_pos:.3e}% (positive) '
            f'and {recovery_error_neg:.3e}% (negative) error on '
            'epsilon_shear_los (M1 mock refits on the current observer likelihood)'
        ),
    }
    
    write_json(step_json_path(STEP_ID), results)
    print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
    return results


if __name__ == "__main__":
    run()
