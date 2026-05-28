#!/usr/bin/env python3
"""Step 011: BAO Acoustic Projection - Strictly Empirical Mode."""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from c0_common import (
    TEPLogger, ensure_dirs, print_status, read_json, 
    rounded, set_step_logger, step_json_path, step_csv_path, 
    write_csv, write_json, PROCESSED_DIR
)

STEP_ID = "step_06_01_bao_projection"

def run():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    step022 = read_json(step_json_path("step_03_01_three_model_comparison"))
    m1_key = "M1_NoLambda_zT1" if "M1_NoLambda_zT1" in step022.get("models", {}) else "M1_NoLambda_zT5"
    m1 = step022['models'][m1_key]['parameters_mle']
    H0 = 70.0  # Dimensionless-distance models fix H0_ref; no H0 parameter in fit

    # Load real BAO data
    bao_path = PROCESSED_DIR / "tep_c0_bao_uncorrelated_compilation.csv"
    if not bao_path.exists():
        raise FileNotFoundError(f"BAO compilation missing at {bao_path}")
    bao_df = pd.read_csv(bao_path)

    # TEP prediction for BAO angular scale diagnostic.
    # theta_BAO = r_s,drag / D_M(z)
    z_drag = 1059.0
    r_s_drag = 147.1 # Mpc
    z_grid = np.linspace(0.1, 2.5, 100)
    # Comoving distance via proper FLRW integration
    # Note: CosmologyFLRW includes radiation component (Or0 computed from CMB temperature)
    # This uses the FLRW framework with proper radiation and curvature handling
    from core.cosmology import CosmologyFLRW
    cosmo = CosmologyFLRW(H0=H0, Om0=m1.get('Om0', 0.3), Ok0=0.0)  # Explicitly flat universe
    d_m = cosmo.luminosity_distance(z_grid) / (1 + z_grid)  # D_M = D_L/(1+z)
    theta_tep = np.divide(r_s_drag, np.maximum(d_m, np.finfo(float).tiny))

    results = {
        'step': STEP_ID,
        'description': 'BAO angular-scale diagnostic; not a covariance likelihood',
        'metrics': {
            'z_drag_reference': z_drag,
            'r_s_drag_mpc': r_s_drag,
            'n_bao_points': len(bao_df),
            'max_theta_tep_deg': rounded(np.max(theta_tep) * 180 / np.pi, 4)
        },
        'validation': {
            'uses_real_bao_compilation': True,
            'research_grade_bao': False,
            'claim_gate': 'blocked',
            'blockers': [
                'Current BAO step is an angular-scale diagnostic, not a full BAO likelihood with covariance.',
                'Research-grade BAO requires D_M/r_d, D_H/r_d, correlations, and survey covariance propagation.',
            ],
        },
    }
    
    write_json(step_json_path(STEP_ID), results)
    print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
    return results

if __name__ == "__main__":
    run()
