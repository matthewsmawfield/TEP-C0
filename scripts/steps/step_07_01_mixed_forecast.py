#!/usr/bin/env python3
"""Step 002: Mixed f_T Forecast - compute evolution of TEP fraction from fitted model.

Uses fitted M1 parameters to predict how TEP fraction evolves with redshift.
Computes actual predictions, not templates.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
from c0_common import TEPLogger, ensure_dirs, print_status, read_json, rounded, set_step_logger, step_csv_path, step_json_path, write_csv, write_json

STEP_ID = "step_07_01_mixed_forecast"


def load_step022_results() -> dict:
    results_file = Path("results/step_03_01_three_model_comparison.json")
    if not results_file.exists():
        raise FileNotFoundError("step_022 results not found. Run step_022 first.")
    return read_json(results_file)


def compute_ft_evolution(z: np.ndarray, ft: float, z_transition: float = 0.5) -> np.ndarray:
    """Compute TEP fraction evolution with redshift.
    
    f_T(z) = ft * (1 + tanh((z - z_transition) / 0.3)) / 2
    
    This models TEP becoming more dominant at higher redshift.
    """
    return ft * (1 + np.tanh((z - z_transition) / 0.3)) / 2


def run():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    step022 = load_step022_results()
    m1_key = "M1_free_zT" if "M1_free_zT" in step022.get("models", {}) else "M1_NoLambda_zT5"
    m1_params = step022['models'][m1_key]['parameters_mle']
    ft = m1_params.get('ft', m1_params.get('epsilon_T', 0.1))  # epsilon_T is the TEP coupling parameter
    
    z_grid = np.linspace(0.0, 3.0, 100)
    ft_evolution = compute_ft_evolution(z_grid, ft)
    
    results = {
        'step': STEP_ID,
        'description': 'TEP fraction evolution forecast from fitted M1 model',
        'model_parameters': {'ft_at_z0': rounded(ft, 3)},
        'predictions': {
            'z_grid': z_grid.tolist(),
            'ft_evolution': ft_evolution.tolist(),
            'ft_at_z1': rounded(compute_ft_evolution(np.array([1.0]), ft)[0], 3),
            'ft_at_z2': rounded(compute_ft_evolution(np.array([2.0]), ft)[0], 3)
        },
        'interpretation': f'TEP fraction f_T evolves from {rounded(ft_evolution[0], 3)} at z=0 to {rounded(ft_evolution[-1], 3)} at z=3'
    }
    
    write_json(step_json_path(STEP_ID), results)
    
    csv_data = [{"z": rounded(z_grid[i], 3), "f_T": rounded(ft_evolution[i], 4)} for i in range(len(z_grid))]
    write_csv(step_csv_path(STEP_ID), csv_data)
    
    print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
    return results


if __name__ == "__main__":
    run()
