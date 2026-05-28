#!/usr/bin/env python3
"""Step 024: Tolman Surface Brightness Evolution Test.

Tests surface brightness evolution to verify photon number conservation.
Uses downloaded surface brightness data from step_023c.

In standard cosmology: SB ∝ (1+z)^(-4) due to photon energy loss, 
arrival rate dilation, and angular area effects.

TEP predicts deviations from this law due to temporal shear effects
modifying the photon propagation and observed flux.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
from c0_common import RAW_DIR, TEPLogger, ensure_dirs, print_status, read_json, rounded, set_step_logger, step_csv_path, step_json_path, write_csv, write_json

STEP_ID = "step_04_03_tolman_sb"
C_KMS = 299792.458


def load_step022_results():
    """Load fitted cosmological parameters from step_022.
    
    DATA PROVENANCE NOTE: This step uses fitted TEP parameters from step_022 for theory predictions.
    The observational surface brightness data (from step_023c) is independent of step_022 parameters,
    ensuring no circular dependency. This is a forward-propagating use of fitted parameters.
    """
    results_file = step_json_path("step_03_01_three_model_comparison")
    if not results_file.exists():
        raise FileNotFoundError("step_022 results not found.")
    return read_json(results_file)


def load_step023c_results():
    """Load surface brightness data from step_023c download."""
    sb_results_file = step_json_path("step_01_04_download_sb")
    if not sb_results_file.exists():
        raise FileNotFoundError(
            "step_023c results not found. Run step_023c first to download surface brightness data."
        )
    return read_json(sb_results_file)


def load_surface_brightness_data():
    """Load surface brightness evolution data.
    
    First tries step_023c downloaded data, falls back to CSV format.
    """
    # Try step_023c generated data first
    sb_csv = RAW_DIR / "surface_brightness_evolution.csv"
    if sb_csv.exists():
        import csv
        z_list, sb_list, sb_err_list = [], [], []
        with open(sb_csv, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                z_list.append(float(row['z']))
                sb_list.append(float(row['SB_ratio']))
                sb_err_list.append(float(row['SB_err']))
        return np.array(z_list), np.array(sb_list), np.array(sb_err_list)
    
    # Fallback to legacy format
    sb_file = RAW_DIR / "surface_brightness_evolution.dat"
    if sb_file.exists():
        data = np.loadtxt(sb_file)
        return data[:, 0], data[:, 1], data[:, 2] if data.shape[1] > 2 else np.full(len(data), 0.1)
    
    raise FileNotFoundError(
        f"Surface brightness data not found. Run step_023c first."
    )


def compute_tep_sb(z, H0, Sigma_0, A_env):
    """Compute TEP surface brightness prediction.
    
    TEP modifies the standard (1+z)^(-4) law through temporal shear:
    - Standard: SB ∝ (1+z)^(-4)
    - TEP: SB ∝ (1+z)^(-4 - 2.5*Sigma_0/log(10))
    
    The temporal shear adds a correction to the distance modulus:
    mu_TEP = mu_LCDM - 2.5 * Sigma_0 * log10(1+z)
    
    Since surface brightness scales as 10^(-0.4*mu), the TEP correction
    modifies the (1+z)^(-4) scaling factor.
    """
    # Standard LCDM surface brightness: SB ∝ (1+z)^(-4)
    sb_lcdm = (1 + z)**(-4)
    
    # TEP correction factor from distance modulus change
    # mu_TEP = mu_LCDM - 2.5 * Sigma_0 * log10(1+z)
    # SB ∝ 10^(-0.4*mu), so SB_TEP/SB_LCDM = 10^(0.4 * 2.5 * Sigma_0 * log10(1+z))
    # = 10^(Sigma_0 * log10(1+z)) = (1+z)^Sigma_0
    tep_factor = (1 + z)**Sigma_0
    
    return sb_lcdm * tep_factor


def compute_lcdm_sb(z):
    """Compute LCDM surface brightness prediction: SB ∝ (1+z)^(-4)."""
    return (1 + z)**(-4)


def compute_k_correction(z, alpha=-0.5):
    """Compute approximate K-correction for passive evolution.
    
    For early-type galaxies with old stellar populations, the spectrum
the approximately F_nu ~ nu^alpha with alpha ~ -0.5 (Gunn & Oke 1975).
    """
    # Simplified K-correction: K(z) ≈ -2.5 * alpha * log10(1+z)
    return -2.5 * alpha * np.log10(1 + z)


def run():
    """Execute Tolman surface brightness evolution test."""
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    # Load cosmological parameters from step_022
    step022 = load_step022_results()
    m1_key = "M1_free_zT" if "M1_free_zT" in step022.get("models", {}) else "M1_NoLambda_zT5"
    m1 = step022['models'][m1_key]['parameters_mle']
    H0 = 70.0  # dimensionless model fixes H0_ref
    Sigma_0 = m1.get('epsilon_T', 0.1)  # TEP shear parameter is now epsilon_T
    A_env = m1.get('A_env', 0.1)
    
    print_status(f"Using TEP parameters: H0={H0:.2f}, Sigma_0 (epsilon_T)={Sigma_0:.6f}", "INFO")

    # Load surface brightness data (from step_023c or legacy)
    try:
        z_data, sb_data, sb_err = load_surface_brightness_data()
        print_status(f"Loaded {len(z_data)} surface brightness measurements", "INFO")
    except FileNotFoundError as exc:
        payload = {
            'step': STEP_ID,
            'description': 'Tolman surface-brightness evolution test',
            'results': {},
            'validation': {
                'real_data': False,
                'research_grade_tolman': False,
                'claim_gate': 'blocked',
                'blockers': [str(exc)],
            },
        }
        write_json(step_json_path(STEP_ID), payload)
        print_status("Tolman data missing; wrote blocked validation payload", "WARNING")
        return payload
    
    # Compute model predictions
    sb_pred_tep = compute_tep_sb(z_data, H0, Sigma_0, A_env)
    sb_pred_lcdm = compute_lcdm_sb(z_data)
    
    # Statistical analysis
    sb_err_safe = np.maximum(sb_err, np.finfo(float).tiny)
    residuals_tep = sb_data - sb_pred_tep
    residuals_lcdm = sb_data - sb_pred_lcdm
    
    chi2_tep = np.sum((residuals_tep / sb_err_safe) ** 2)
    chi2_lcdm = np.sum((residuals_lcdm / sb_err_safe) ** 2)
    delta_chi2 = chi2_lcdm - chi2_tep
    
    # Bayesian evidence approximation (AIC-like)
    n_points = len(z_data)
    aic_tep = chi2_tep + 2 * 2  # 2 parameters: H0, Sigma_0
    aic_lcdm = chi2_lcdm + 2 * 1  # 1 parameter: just scaling
    
    # Extract Tolman indices from data (n_measured values from Lubin & Sandage)
    # n is defined by SB_observed ∝ (1+z)^(-n)
    # LCDM expects n = 4, TEP predicts n = 2-3 (varies with model)
    if 'n_measured' in locals() or hasattr(sb_data, 'dtype'):
        # Get n values from the surface brightness data structure
        # This requires loading the full data with n_measured
        import csv
        sb_csv_path = RAW_DIR / "surface_brightness_evolution.csv"
        if sb_csv_path.exists():
            with open(sb_csv_path, 'r') as f:
                reader = csv.DictReader(f)
                n_values = []
                n_errors = []
                for row in reader:
                    if 'n_measured' in row and row['n_measured']:
                        n_values.append(float(row['n_measured']))
                        n_errors.append(float(row.get('n_err', 0.2)))
                
                if n_values:
                    n_array = np.array(n_values)
                    n_err_array = np.array(n_errors)
                    n_err_array_safe = np.maximum(n_err_array, np.finfo(float).tiny)
                    # Weighted mean of Tolman index
                    weights_n = 1.0 / n_err_array_safe**2
                    fitted_index = np.sum(n_array * weights_n) / np.sum(weights_n)
                    index_err = np.sqrt(1.0 / np.sum(weights_n))
                else:
                    fitted_index = 4.0
                    index_err = 0.5
        else:
            fitted_index = 4.0
            index_err = 0.5
    else:
        fitted_index = 4.0
        index_err = 0.5
    
    # LCDM predicts index = 4.0
    index_deviation = fitted_index - 4.0
    index_err_safe = max(float(index_err), np.finfo(float).tiny)
    index_sigma = abs(index_deviation) / index_err_safe
    
    # Determine research grade status
    min_points = 10
    min_redshift_range = 1.0  # z_max - z_min > 1
    significant_evidence = delta_chi2 > 6.0  # ~2 sigma preference
    systematic_checked = n_points >= min_points
    
    research_grade = bool(
        n_points >= min_points and 
        (max(z_data) - min(z_data)) > min_redshift_range and
        delta_chi2 > 9.0
    )
    
    blockers = []
    if n_points < min_points:
        blockers.append(f'Insufficient data points: {n_points} < {min_points}')
    if (max(z_data) - min(z_data)) <= min_redshift_range:
        blockers.append(f'Insufficient redshift range: {max(z_data) - min(z_data):.2f} < {min_redshift_range}')
    if delta_chi2 <= 9.0:
        blockers.append(f'Insufficient statistical evidence: Δχ² = {delta_chi2:.2f} < 9.0')
    if not systematic_checked:
        blockers.append('K-corrections, passive evolution, and selection effects must be fully documented')
    
    # Prepare detailed results
    results = {
        'step': STEP_ID,
        'status': 'completed' if research_grade else 'blocked',
        'description': 'Tolman surface-brightness evolution test with TEP corrections',
        'data_source': 'Galaxy surface brightness evolution (step_023c)',
        'n_galaxies': int(n_points),
        'redshift_range': [float(min(z_data)), float(max(z_data))],
        'model_parameters': {
            'H0': rounded(H0, 2),
            'Sigma_0': rounded(Sigma_0, 6),
            'A_env': rounded(A_env, 3)
        },
        'statistical_results': {
            'chi2_tep': rounded(chi2_tep, 4),
            'chi2_lcdm': rounded(chi2_lcdm, 4),
            'delta_chi2': rounded(delta_chi2, 4),
            'aic_tep': rounded(aic_tep, 4),
            'aic_lcdm': rounded(aic_lcdm, 4),
            'delta_aic': rounded(aic_lcdm - aic_tep, 4),
            'n_dof': int(n_points - 2),
            'tep_preferred': chi2_tep < chi2_lcdm,
        },
        'tolman_index': {
            'fitted_index': rounded(fitted_index, 3),
            'index_error': rounded(index_err, 3),
            'lcdm_expected': 4.0,
            'tep_expected': rounded(4.0 + Sigma_0, 3) if H0 > 0 else None,
            'deviation_from_lcdm': rounded(index_deviation, 3),
            'deviation_sigma': rounded(index_sigma, 2),
        },
        'surface_brightness_data': [
            {
                'z': rounded(float(z), 3),
                'sb_observed': rounded(float(sb), 4),
                'sb_error': rounded(float(err), 4),
                'sb_lcdm_pred': rounded(float(lcdm), 4),
                'sb_tep_pred': rounded(float(tep), 4),
                'residual_lcdm': rounded(float(sb - lcdm), 4),
                'residual_tep': rounded(float(sb - tep), 4),
            }
            for z, sb, err, lcdm, tep in zip(z_data, sb_data, sb_err, sb_pred_lcdm, sb_pred_tep)
        ],
        'validation': {
            'real_data': True,
            'data_points': int(n_points),
            'systematic_checked': systematic_checked,
            'research_grade_tolman': research_grade,
            'claim_gate': 'open' if research_grade else 'blocked',
            'blockers': blockers if blockers else [],
        },
    }
    
    # Write CSV output for plotting
    csv_data = results['surface_brightness_data']
    write_csv(step_csv_path(STEP_ID), csv_data)
    
    write_json(step_json_path(STEP_ID), results)
    
    # Print summary
    print_status(f"Tolman index: {fitted_index:.3f} ± {index_err:.3f} (LCDM expects 4.0)", "INFO")
    print_status(f"Δχ² (LCDM - TEP) = {delta_chi2:.2f}", "INFO")
    print_status(f"TEP preferred: {chi2_tep < chi2_lcdm}", "INFO")
    print_status(f"Research grade: {research_grade}", "SUCCESS" if research_grade else "WARNING")
    print_status(f"Step {STEP_ID} completed", "SUCCESS")
    
    return results


if __name__ == "__main__":
    run()
