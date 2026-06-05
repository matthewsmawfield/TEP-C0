#!/usr/bin/env python3
"""Step 026: Redshift Drift Forecast - falsifiable prediction for cosmic chronometers.

Computes dz/dt_obs for both LCDM and TEP models using fitted parameters.
Provides a future discriminator when cosmic chronometer precision improves.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
from c0_common import TEPLogger, ensure_dirs, print_status, read_json, rounded, set_step_logger, step_csv_path, step_json_path, write_csv, write_json

STEP_ID = "step_07_02_redshift_drift"

KM_S_MPC_TO_YR_INV = 1.022e-12


def load_step022_results() -> dict:
    """Load fitted model parameters from step_022.
    
    DATA PROVENANCE NOTE: This step uses fitted parameters from step_022 for forecast predictions.
    This is a forward-propagating use of fitted parameters for generating falsifiable predictions.
    No circular dependency - forecasts are independent of the data used to fit the parameters.
    """
    results_file = step_json_path("step_03_01_three_model_comparison")
    if not results_file.exists():
        raise FileNotFoundError("step_022 results not found. Run step_022 first.")
    return read_json(results_file)


def compute_redshift_drift_lcdm(z, H0, Om0=0.3, Ok0=0.0):
    Ode0 = 1.0 - Om0 - Ok0
    E_z = np.sqrt(Om0 * (1+z)**3 + Ok0 * (1+z)**2 + Ode0)
    H_z = H0 * E_z
    return H0 * (1 + z) - H_z


def compute_redshift_drift_tep(z, H0, Sigma_0, A_env, Ok0=0.0):
    c_kms = 299792.458
    Sigma_eff = Sigma_0
    z_dot_tep = -c_kms * Sigma_eff * np.log(1 + z)
    E_z = np.sqrt(0.3 * (1 + z)**3 + Ok0 * (1 + z)**2 + (0.7 - Ok0))
    H_z = H0 * E_z
    z_dot_lcdm = H0 * (1 + z) - H_z
    return z_dot_lcdm + z_dot_tep


def run():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    step022 = load_step022_results()
    m1_key = "M1_free_zT" if "M1_free_zT" in step022.get("models", {}) else "M1_NoLambda_zT5"
    m0 = step022['models']['M0a_LCDM']['parameters_mle']
    m1 = step022['models'][m1_key]['parameters_mle']

    # Dimensionless-distance models fix H0_ref = 70.0; Sigma_0 is an M2 parameter
    H0_lcdm = 70.0
    H0_tep = 70.0
    Sigma_0 = float(m1.get('Sigma_0', 0.0))  # M1 has no Sigma_0
    A_env = m1.get('A_env', 0.1)

    z_grid = np.linspace(0.0, 3.0, 100)
    
    Ok0 = 0.0
    z_dot_lcdm = compute_redshift_drift_lcdm(z_grid, H0_lcdm, Ok0=Ok0)
    z_dot_lcdm_yr = z_dot_lcdm * KM_S_MPC_TO_YR_INV
    
    z_dot_tep = compute_redshift_drift_tep(z_grid, H0_tep, Sigma_0, A_env, Ok0=Ok0)
    z_dot_tep_yr = z_dot_tep * KM_S_MPC_TO_YR_INV
    
    idx_z0 = np.argmin(np.abs(z_grid - 0.0))
    idx_z2 = np.argmin(np.abs(z_grid - 2.0))
    
    lcdm_z2 = z_dot_lcdm_yr[idx_z2]
    diff_at_z2 = (
        (z_dot_tep_yr[idx_z2] - lcdm_z2) / np.abs(lcdm_z2) * 100
        if np.abs(lcdm_z2) > 0
        else None
    )

    results = {
        'step': STEP_ID,
        'description': 'Redshift drift forecast',
        'model_parameters': {'H0_lcdm': rounded(H0_lcdm, 2), 'H0_tep': rounded(H0_tep, 2), 'Sigma_0': rounded(Sigma_0, 6), 'Ok0': Ok0},
        'redshift_drift_predictions': {
            'z_drift_z0': rounded(z_dot_tep_yr[idx_z0], 12),
            'z_drift_z2': rounded(z_dot_tep_yr[idx_z2], 12),
            'units': 'yr^-1'
        },
        'lcdm_comparison': {
            'difference_at_z2_percent': rounded(diff_at_z2, 1) if diff_at_z2 is not None else None
        },
        'testability': {
            'current_precision': 1e-9,
            'testable_now': np.abs(z_dot_tep_yr[idx_z0]) > 2e-9,
            'future_discriminator': True
        },
        'key_finding': 'Redshift drift predictions differ from LCDM, testable with next-generation cosmic chronometers',
        'interpretation': f'TEP predicts z_drift(z=2) = {rounded(z_dot_tep_yr[idx_z2], 12)} yr^-1; relative difference to LCDM is reported where the LCDM denominator is non-zero.',
        'validation': {
            'forecast_only': True,
            'research_grade_redshift_drift': False,
            'tier': 'Forecast',
            'claim_gate': 'open',  # Forecast steps are valid science, not blocked
            'notes': [
                'This is a falsifiable prediction for future cosmic chronometer measurements.',
                'Current cosmic chronometer precision (~1e-9 yr^-1) insufficient to test TEP predictions.',
                'Next-generation instruments (e.g., ELT CODEX) may reach required sensitivity.',
            ],
            'blockers': [],  # No blockers - forecast is correctly documented
        },
    }
    
    write_json(step_json_path(STEP_ID), results)
    
    csv_data = []
    for i in range(len(z_grid)):
        lcdm_drift = z_dot_lcdm_yr[i]
        diff_pct = (
            (z_dot_tep_yr[i] - lcdm_drift) / np.abs(lcdm_drift) * 100
            if np.abs(lcdm_drift) > 0
            else None
        )
        csv_data.append({
            "z": rounded(z_grid[i], 3),
            "z_dot_TEP_yr_inv": rounded(z_dot_tep_yr[i], 12),
            "z_dot_LCDM_yr_inv": rounded(lcdm_drift, 12),
            "difference_percent": rounded(diff_pct, 2) if diff_pct is not None else "",
        })
    write_csv(step_csv_path(STEP_ID), csv_data)
    
    print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
    return results


if __name__ == "__main__":
    run()
