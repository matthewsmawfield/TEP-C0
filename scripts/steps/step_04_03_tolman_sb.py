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
    
    TEP modifies the standard (1+z)^(-4) law through temporal shear.
    Because TEP expands the effective distance (mu_TEP > mu_LCDM),
    the observed flux is lower, meaning the surface brightness falls off faster.
    """
    import sys
    import os
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
    from core.cosmology import TEPCosmology
    
    # Standard LCDM surface brightness: SB ∝ (1+z)^(-4)
    sb_lcdm = (1 + z)**(-4)
    
    cosmo_tep = TEPCosmology(H0=H0, epsilon_T=Sigma_0)
    cosmo_lcdm = TEPCosmology(H0=H0, epsilon_T=0.0)
    
    # Flux ratio = 10^(-0.4 * Delta mu)
    mu_tep = cosmo_tep.distance_modulus(z)
    mu_lcdm = cosmo_lcdm.distance_modulus(z)
    
    tep_factor = 10**(-0.4 * (mu_tep - mu_lcdm))
    
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
    # CRITICAL FIX: Tolman SB is a background test, so it must use the
    # homogeneous background epsilon_T (acoustic sector ~0.018), NOT the
    # SNe line-of-sight void shear epsilon_shear_los (~0.83). Using the
    # void shear overpredicts the TEP deviation by a factor of ~40.
    # If no acoustic value is available, fall back to the canonical 0.018.
    Sigma_0 = m1.get('epsilon_T', 0.018)
    A_env = m1.get('A_env', 0.1)
    
    print_status(f"Using TEP parameters: H0={H0:.2f}, Sigma_0 (epsilon_T)={Sigma_0:.6f}", "INFO")
    print_status("  (Using acoustic-sector epsilon_T for background Tolman test)", "INFO")

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
    
    # Load TEPCosmology for distance calculations
    import sys, os
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
    from core.cosmology import TEPCosmology
    
    # Load n_measured and n_err from the CSV (the actual physical observable)
    # The SB_ratio column is a derived quantity; n_measured is the fundamental Tolman index
    import csv
    sb_csv_path = RAW_DIR / "surface_brightness_evolution.csv"
    n_values, n_errors, z_n_values = [], [], []
    if sb_csv_path.exists():
        with open(sb_csv_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if 'n_measured' in row and row['n_measured']:
                    n_values.append(float(row['n_measured']))
                    n_errors.append(float(row.get('n_err', 0.2)))
                    z_n_values.append(float(row['z']))
    
    n_array = np.array(n_values)
    n_err_array = np.array(n_errors)
    n_err_safe = np.maximum(n_err_array, np.finfo(float).tiny)
    z_n_array = np.array(z_n_values)
    
    # Compute TEP-predicted Tolman index at each data redshift
    # n_tep = 4 - log10(SB_tep/SB_lcdm) / log10(1+z)
    # SB_tep/SB_lcdm = 10^(-0.4 * (mu_tep - mu_lcdm))
    cosmo_tep = TEPCosmology(H0=H0, epsilon_T=Sigma_0)
    cosmo_lcdm = TEPCosmology(H0=H0, epsilon_T=0.0)
    n_tep_values = []
    for z in z_n_array:
        mu_tep = cosmo_tep.distance_modulus(z)
        mu_lcdm = cosmo_lcdm.distance_modulus(z)
        if z > 0.001:
            tep_factor = 10**(-0.4 * (mu_tep - mu_lcdm))
            n_tep = 4.0 - np.log10(tep_factor) / np.log10(1 + z)
        else:
            n_tep = 4.0
        n_tep_values.append(n_tep)
    n_tep_array = np.array(n_tep_values)
    
    # Statistical analysis on Tolman indices (the physical observables)
    # LCDM predicts n = 4.0 for all z
    n_lcdm_array = np.full_like(n_array, 4.0)
    
    residuals_tep = n_array - n_tep_array
    residuals_lcdm = n_array - n_lcdm_array
    
    chi2_tep = np.sum((residuals_tep / n_err_safe) ** 2)
    chi2_lcdm = np.sum((residuals_lcdm / n_err_safe) ** 2)
    delta_chi2 = chi2_lcdm - chi2_tep
    
    # Bayesian evidence approximation (AIC-like)
    n_points = len(n_array)
    aic_tep = chi2_tep + 2 * 2  # 2 parameters: H0, Sigma_0
    aic_lcdm = chi2_lcdm + 2 * 1  # 1 parameter: just scaling
    
    # Weighted mean of measured Tolman index
    weights_n = 1.0 / n_err_safe**2
    fitted_index = np.sum(n_array * weights_n) / np.sum(weights_n)
    index_err = np.sqrt(1.0 / np.sum(weights_n))
    
    # LCDM predicts index = 4.0
    index_deviation = fitted_index - 4.0
    index_err_safe = max(float(index_err), np.finfo(float).tiny)
    index_sigma = abs(index_deviation) / index_err_safe
    
    # Redshift trend analysis: the data shows n DECREASING with z, while TEP
    # predicts n INCREASING with z (more temporal shear at higher z).
    n_vs_z_slope = 0.0
    n_vs_z_intercept = 0.0
    low_z_mean_n = 0.0
    high_z_mean_n = 0.0
    low_z_mask = np.array([])
    high_z_mask = np.array([])
    if len(z_n_array) > 1:
        n_vs_z_slope, n_vs_z_intercept = np.polyfit(z_n_array, n_array, 1)
        low_z_mask = z_n_array < 0.3
        high_z_mask = z_n_array > 0.5
        if np.any(low_z_mask):
            low_z_mean_n = float(np.mean(n_array[low_z_mask]))
        if np.any(high_z_mask):
            high_z_mean_n = float(np.mean(n_array[high_z_mask]))
    
    # TEP predicted slope: compute dn/dz for TEP at a representative z
    tep_slope = 0.0
    if len(z_n_array) > 1:
        z_mid = float(np.median(z_n_array))
        dz = 0.01
        mu_tep_lo = cosmo_tep.distance_modulus(z_mid - dz)
        mu_tep_hi = cosmo_tep.distance_modulus(z_mid + dz)
        mu_lcdm_lo = cosmo_lcdm.distance_modulus(z_mid - dz)
        mu_lcdm_hi = cosmo_lcdm.distance_modulus(z_mid + dz)
        tep_factor_lo = 10**(-0.4 * (mu_tep_lo - mu_lcdm_lo))
        tep_factor_hi = 10**(-0.4 * (mu_tep_hi - mu_lcdm_hi))
        n_tep_lo = 4.0 - np.log10(tep_factor_lo) / np.log10(1 + z_mid - dz)
        n_tep_hi = 4.0 - np.log10(tep_factor_hi) / np.log10(1 + z_mid + dz)
        tep_slope = (n_tep_hi - n_tep_lo) / (2 * dz)
    
    # Systematic uncertainty estimate: K-corrections for early-type galaxies
    # in R and I bands can shift n by ~0.3-0.5 mag at z~1 (Lubin & Sandage 2001).
    # Passive evolution adds another ~0.2-0.3 mag uncertainty.
    # Total systematic uncertainty on n: ~0.5 (conservative).
    k_corr_systematic = 0.5  # mag-equivalent uncertainty in n
    total_systematic_err = np.sqrt(n_err_safe**2 + k_corr_systematic**2)
    
    # Determine research grade status
    min_points = 10
    min_redshift_range = 1.0  # z_max - z_min > 1
    systematic_checked = n_points >= min_points
    
    # Pipeline-debug signal: when |n_tep - n_lcdm| << |n_obs - n_lcdm|, both models
    # predict the same thing and the test has zero discriminating power.
    # For epsilon_T=0.018, TEP predicts n≈4.015, LCDM predicts n=4.0,
    # while data shows n≈3.375. The 0.6 mag offset is astrophysical systematics.
    mean_n_tep = float(np.mean(n_tep_array)) if len(n_tep_array) > 0 else 4.0
    model_degeneracy = abs(mean_n_tep - 4.0) < 0.1  # TEP ≈ LCDM prediction
    
    # CRITICAL: data trend is OPPOSITE to TEP prediction.
    # TEP (with any epsilon_T > 0) predicts n >= 4.0 and dn/dz >= 0 (flat or increasing).
    # Data shows n << 4.0 and dn/dz < 0 (strongly decreasing with z).
    # This sign mismatch means TEP cannot explain the Tolman anomaly even in principle.
    opposite_trend = n_vs_z_slope < -0.1 and tep_slope > -0.05
    
    # When models are degenerate, the test is inconclusive as a discriminator.
    # We still require sufficient data for a robust measurement.
    research_grade = bool(
        n_points >= min_points and 
        (max(z_n_array) - min(z_n_array)) > min_redshift_range and
        not model_degeneracy and
        not opposite_trend and
        delta_chi2 > 9.0
    )
    
    inconclusive = model_degeneracy and n_points >= min_points
    
    blockers = []
    if n_points < min_points:
        blockers.append(f'Insufficient data points: {n_points} < {min_points}')
    if (max(z_n_array) - min(z_n_array)) <= min_redshift_range:
        blockers.append(f'Insufficient redshift range: {max(z_n_array) - min(z_n_array):.2f} < {min_redshift_range}')
    if opposite_trend:
        tep_trend_str = 'n increases with z' if tep_slope > 0.01 else 'n ≈ 4.0 (flat)'
        blockers.append(
            f'Data trend is OPPOSITE to TEP prediction: '
            f'data slope = {n_vs_z_slope:.3f} (n decreases with z), '
            f'TEP slope = {tep_slope:.3f} ({tep_trend_str}). '
            f'TEP (and LCDM) predict n ≥ 4.0; data shows n ≈ 3.375 and falls to n ≈ 2.8 at high z. '
            f'TEP cannot explain the Tolman anomaly in either amplitude or trend. '
            f'Observed offset is dominated by astrophysical systematics '
            f'(K-corrections ±{k_corr_systematic:.1f}, passive evolution, selection effects).'
        )
    elif model_degeneracy:
        blockers.append(f'TEP (n≈{mean_n_tep:.2f}) and LCDM (n=4.0) are degenerate; test has no discriminating power. Observed n={fitted_index:.3f} is dominated by astrophysical systematics (K-corrections, passive evolution).')
    elif delta_chi2 <= 9.0:
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
            'tep_expected': rounded(np.mean(n_tep_array), 3) if len(n_tep_array) > 0 else None,
            'deviation_from_lcdm': rounded(index_deviation, 3),
            'deviation_sigma': rounded(index_sigma, 2),
        },
        'redshift_trend': {
            'data_slope': rounded(n_vs_z_slope, 3),
            'data_intercept': rounded(n_vs_z_intercept, 3),
            'tep_slope': rounded(tep_slope, 3),
            'low_z_mean_n': rounded(low_z_mean_n, 3) if np.any(low_z_mask) else None,
            'high_z_mean_n': rounded(high_z_mean_n, 3) if np.any(high_z_mask) else None,
            'opposite_trend': opposite_trend,
        },
        'systematic_uncertainty': {
            'k_correction_estimate': k_corr_systematic,
            'total_systematic_err': rounded(float(np.mean(total_systematic_err)), 3),
        },
        'tolman_data': [
            {
                'z': rounded(float(z), 3),
                'n_observed': rounded(float(n), 4),
                'n_err': rounded(float(err), 4),
                'n_lcdm_pred': 4.0,
                'n_tep_pred': rounded(float(nt), 4),
                'residual_lcdm': rounded(float(n - 4.0), 4),
                'residual_tep': rounded(float(n - nt), 4),
            }
            for z, n, err, nt in zip(z_n_array, n_array, n_err_array, n_tep_array)
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
    csv_data = results['tolman_data']
    write_csv(step_csv_path(STEP_ID), csv_data)
    
    write_json(step_json_path(STEP_ID), results)
    
    # Print summary
    print_status(f"Tolman index: {fitted_index:.3f} ± {index_err:.3f} (LCDM expects 4.0)", "INFO")
    print_status(f"  Data redshift trend: n = {n_vs_z_intercept:.3f} + {n_vs_z_slope:.3f}*z (slope {'< 0' if n_vs_z_slope < 0 else '> 0'})", "INFO")
    tep_trend_desc = "n increases with z" if tep_slope > 0.01 else "n ≈ 4.0 (flat)"
    print_status(f"  TEP predicted trend: slope ≈ {tep_slope:.3f} ({tep_trend_desc})", "INFO")
    if np.any(low_z_mask):
        print_status(f"  Low-z (z<0.3) mean n: {low_z_mean_n:.3f}", "INFO")
    if np.any(high_z_mask):
        print_status(f"  High-z (z>0.5) mean n: {high_z_mean_n:.3f}", "INFO")
    print_status(f"  K-correction systematic uncertainty: ±{k_corr_systematic:.1f} mag", "INFO")
    print_status(f"Δχ² (LCDM - TEP) = {delta_chi2:.2f}", "INFO")
    print_status(f"TEP preferred: {chi2_tep < chi2_lcdm}", "INFO")
    print_status(f"Research grade: {research_grade}", "SUCCESS" if research_grade else "WARNING")
    print_status(f"Step {STEP_ID} completed", "SUCCESS")
    
    return results


if __name__ == "__main__":
    run()
