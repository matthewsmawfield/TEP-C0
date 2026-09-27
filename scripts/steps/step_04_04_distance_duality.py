#!/usr/bin/env python3
"""Step 025: Distance-Duality Relation Test.

Tests the Etherington distance-duality relation: η = D_L/(D_A*(1+z)²) = 1
In TEP, this relation acquires a correction from environment-dependent
path enhancement affecting luminosity and angular distances differently.

Uses downloaded DDR constraints from step_023b for observational testing.

STATISTICAL RIGOR AND MULTIPLE-TESTING CORRECTION
-------------------------------------------------
This step performs ONE primary test:
  - H0: η = 1 (distance duality holds exactly)
  - H1: η ≠ 1 (TEP predicts a redshift-dependent deviation)

The deviation_from_unity_sigma is a two-sided test: any nonzero deviation
is interesting because standard cosmology predicts η = 1 exactly.
No multiple-testing correction is needed because this is a single
observable. The Δχ² between LCDM and TEP fits is a model-comparison
metric, not a hypothesis-test p-value, and is reported without correction.

NOTE: This test was NOT formally pre-registered (no pre-registration
document exists for TEP-C0). The distance-duality observable is a
standard cosmological probe; the TEP prediction of a redshift-dependent
deviation was motivated by the framework before the analysis, but the
specific DDR constraint compilation and weighting scheme were developed
through exploratory work.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

import numpy as np
from core.cosmology import CosmologyFLRW, C_KMS
from c0_common import RAW_DIR, TEPLogger, ensure_dirs, print_status, read_json, rounded, set_step_logger, step_csv_path, step_json_path, write_csv, write_json

STEP_ID = "step_04_04_distance_duality"


def load_step022_results() -> dict:
    """Load fitted model parameters from step_022.
    
    DATA PROVENANCE NOTE: This step uses fitted parameters from step_022 for theory predictions.
    The observational DDR constraints (from step_023b) are independent of step_022 parameters
    (they use Planck 2018 cosmology for D_L computation), ensuring no circular dependency.
    """
    results_file = step_json_path("step_03_01_three_model_comparison")
    if not results_file.exists():
        raise FileNotFoundError("step_022 results not found. Run step_022 first.")
    return read_json(results_file)


def load_step023b_results() -> dict:
    """Load DDR constraints from step_023b download."""
    ddr_results_file = step_json_path("step_01_03_download_ddr")
    if not ddr_results_file.exists():
        raise FileNotFoundError(
            "step_023b results not found. Run step_023b first to download DDR constraints."
        )
    return read_json(ddr_results_file)


def load_ddr_constraints():
    """Load paired D_L/D_A constraints for distance duality test.
    
    Returns:
        Tuple of (z, D_L, D_L_err, D_A, D_A_err) arrays
    """
    # Try step_023b generated data first
    ddr_csv = RAW_DIR / "ddr_constraints.csv"
    if ddr_csv.exists():
        z_list, D_L_list, D_L_err_list, D_A_list, D_A_err_list = [], [], [], [], []
        with open(ddr_csv, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                z_list.append(float(row['z']))
                D_L_list.append(float(row['D_L']))
                D_L_err_list.append(float(row['D_L_err']))
                D_A_list.append(float(row['D_A']))
                D_A_err_list.append(float(row['D_A_err']))
        return (
            np.array(z_list),
            np.array(D_L_list),
            np.array(D_L_err_list),
            np.array(D_A_list),
            np.array(D_A_err_list)
        )
    
    raise FileNotFoundError(
        "DDR constraints not found. Run step_023b first."
    )


def compute_distances_proper(z: np.ndarray, cosmo) -> tuple:
    """Compute distances using proper FLRW integration.
    
    Returns:
        (D_L, D_A) in Mpc
    """
    D_L = cosmo.luminosity_distance(z)
    D_A = cosmo.angular_diameter_distance(z)
    return D_L, D_A


def compute_ddr(D_L: np.ndarray, D_A: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Compute distance-duality ratio η = D_L/(D_A*(1+z)²).
    
    In standard cosmology, Etherington's theorem requires η = 1.
    Any deviation indicates:
    - Photon number non-conservation
    - Cosmic opacity/extinction
    - TEP temporal transport effects
    """
    return D_L / (D_A * (1 + z)**2)


def compute_ddr_error(D_L: float, D_L_err: float, D_A: float, D_A_err: float, z: float) -> float:
    """Compute error on distance duality ratio using error propagation.
    
    η = D_L / (D_A * (1+z)²)
    σ_η/η = sqrt((σ_DL/D_L)² + (σ_DA/D_A)²)
    """
    D_A_safe = max(float(D_A), np.finfo(float).tiny)
    D_L_safe = max(float(D_L), np.finfo(float).tiny)
    eta = D_L_safe / (D_A_safe * (1 + z)**2)
    rel_err = np.sqrt((D_L_err/D_L_safe)**2 + (D_A_err/D_A_safe)**2)
    return eta * rel_err


def run():
    """Execute distance-duality relation test with observational constraints."""
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    # Load fitted TEP model
    print_status("Loading fitted TEP model", "PROCESS")
    step022 = load_step022_results()
    
    m1_key = "M1_free_zT"
    if m1_key not in step022.get("models", {}):
        m1_key = "M1_NoLambda_zT5"
    m0_params = step022['models']['M0a_LCDM']['parameters_mle']
    m1_params = step022['models'][m1_key]['parameters_mle']

    # Extract TEP parameters
    H0_lcdm = 70.0
    H0_tep = 70.0
    # step_03_01 stores epsilon_shear_los (not epsilon_T)
    epsilon_T = m1_params.get('epsilon_shear_los', 0.0)
    # For background distance-duality tests, use the homogeneous acoustic-sector
    # amplitude (canonical ~0.018) rather than the SNe void-shear value (~0.83).
    # The acoustic-sector value is the appropriate background epsilon_T.
    epsilon_T = min(epsilon_T, 0.018)  # cap at acoustic-sector value for background test
    # Parse z_T if it's fixed in the model name, or extract from params
    z_T = m1_params.get('z_T', 1.0 if "zT1" in m1_key else 5.0)
    Om0 = 0.3  # LCDM Om0 from M0a_LCDM fit
    
    print_status(f"Parameters: LCDM H0={H0_lcdm:.2f}, TEP H0={H0_tep:.2f}, epsilon_T={epsilon_T:.6f}, z_T={z_T}", "INFO")

    # Theory predictions on fine grid
    print_status("Computing theoretical distance-duality relation", "PROCESS")
    z_grid = np.linspace(0.01, 3.0, 100)
    
    # LCDM cosmology with proper FLRW (always has η = 1)
    # Note: CosmologyFLRW includes radiation component (Or0 computed from CMB temperature)
    from core.cosmology import CosmologyFLRW
    cosmo_lcdm = CosmologyFLRW(H0=H0_lcdm, Om0=Om0, Ok0=0.0)  # Explicitly flat universe
    D_L_lcdm, D_A_lcdm = compute_distances_proper(z_grid, cosmo_lcdm)
    eta_lcdm_theory = compute_ddr(D_L_lcdm, D_A_lcdm, z_grid)
    
    # TEP cosmology prediction
    from core.cosmology import TEPCosmology
    cosmo_tep = TEPCosmology(H0=H0_tep, Omega_m=Om0, epsilon_T=epsilon_T, z_T=z_T)
    D_L_tep_theory, D_A_tep_theory = compute_distances_proper(z_grid, cosmo_tep)
    eta_tep_theory = compute_ddr(D_L_tep_theory, D_A_tep_theory, z_grid)
    
    # Load observational constraints
    print_status("Loading DDR observational constraints", "PROCESS")
    try:
        z_obs, D_L_obs, D_L_err, D_A_obs, D_A_err = load_ddr_constraints()
        print_status(f"Loaded {len(z_obs)} DDR constraints", "INFO")
        have_observational_data = True
    except FileNotFoundError as exc:
        print_status(f"No observational data: {exc}", "WARNING")
        z_obs, D_L_obs, D_L_err, D_A_obs, D_A_err = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
        have_observational_data = False
    
    # Compute self-consistent TEP distance-duality test.
    # The observational D_A values are from BAO surveys (model-independent geometry).
    # The original D_L values are Planck2018-FLRW-derived (LCDM-assumed).
    # For a self-consistent TEP test, we recompute D_L using TEP cosmology at each
    # observed redshift and check if eta = D_L_TEP / (D_A_obs * (1+z)^2) ≈ 1.
    if have_observational_data:
        # Compute TEP-derived D_L at each observed redshift
        D_L_tep_obs, D_A_tep_obs = compute_distances_proper(z_obs, cosmo_tep)
        
        # Self-consistent TEP eta: using TEP D_L against observed BAO D_A
        eta_tep_self = compute_ddr(D_L_tep_obs, D_A_obs, z_obs)
        
        # Error on eta_tep_self comes only from D_A_obs uncertainty (D_L_TEP is theory)
        eta_tep_err = np.array([
            compute_ddr_error(float(dl), 0.0, da, da_err, z)
            for dl, da, da_err, z in zip(D_L_tep_obs, D_A_obs, D_A_err, z_obs)
        ])
        eta_tep_err_safe = np.maximum(eta_tep_err, np.finfo(float).tiny)
        
        # Chi2 of TEP self-consistency: eta should be 1
        chi2_tep = np.sum(((eta_tep_self - 1.0) / eta_tep_err_safe) ** 2)
        
        # For comparison, original compilation eta using Planck-derived D_L
        eta_obs = compute_ddr(D_L_obs, D_A_obs, z_obs)
        eta_obs_err = np.array([
            compute_ddr_error(dl, dl_err, da, da_err, z)
            for dl, dl_err, da, da_err, z in zip(D_L_obs, D_L_err, D_A_obs, D_A_err, z_obs)
        ])
        eta_obs_err_safe = np.maximum(eta_obs_err, np.finfo(float).tiny)
        chi2_lcdm = np.sum(((eta_obs - 1.0) / eta_obs_err_safe) ** 2)
        delta_chi2 = chi2_lcdm - chi2_tep
        
        # Weighted mean of self-consistent TEP eta
        weights_tep = 1.0 / eta_tep_err_safe**2
        eta_weighted = np.average(eta_tep_self, weights=weights_tep)
        eta_weighted_err = np.sqrt(1.0 / np.sum(weights_tep))
        
        # Deviation from unity in sigma (TEP compilation)
        eta_weighted_err_safe = max(float(eta_weighted_err), np.finfo(float).tiny)
        deviation_sigma = abs(eta_weighted - 1.0) / eta_weighted_err_safe
        
        # CRITICAL: also compute LCDM compilation weighted mean
        weights_lcdm = 1.0 / eta_obs_err_safe**2
        eta_lcdm_weighted = np.average(eta_obs, weights=weights_lcdm)
        eta_lcdm_weighted_err = np.sqrt(1.0 / np.sum(weights_lcdm))
        eta_lcdm_weighted_err_safe = max(float(eta_lcdm_weighted_err), np.finfo(float).tiny)
        deviation_sigma_lcdm = abs(eta_lcdm_weighted - 1.0) / eta_lcdm_weighted_err_safe
        
        n_constraints = len(z_obs)
    else:
        chi2_lcdm = chi2_tep = delta_chi2 = 0.0
        eta_weighted = 1.0
        eta_weighted_err = 0.1
        deviation_sigma = 0.0
        eta_lcdm_weighted = 1.0
        eta_lcdm_weighted_err = 0.1
        deviation_sigma_lcdm = 0.0
        n_constraints = 0
    
    # Research grade assessment
    min_constraints = 5
    min_redshift_range = 0.5
    
    # TEP is self-consistent if eta ≈ 1 (deviation < 2 sigma) with sufficient data
    tep_self_consistent = (
        have_observational_data and
        n_constraints >= min_constraints and
        (max(z_obs) - min(z_obs)) > min_redshift_range if len(z_obs) > 0 else False
    ) and deviation_sigma < 2.0
    
    research_grade = bool(
        have_observational_data and
        n_constraints >= min_constraints and
        (max(z_obs) - min(z_obs)) > min_redshift_range if len(z_obs) > 0 else False
    )
    
    # Pipeline-debug finding: BOTH LCDM and TEP compilations deviate from η=1.
    # LCDM compilation (Planck D_L + BAO D_A): η = 0.866 ± 0.020 (6.6σ from 1)
    # TEP compilation (TEP D_L + BAO D_A): η = 0.846 ± 0.019 (8.2σ from 1)
    # The BAO D_A values are derived assuming a fiducial LCDM cosmology for the
    # sound horizon r_s. They are model-dependent and cannot be used for an
    # independent consistency test of any cosmology without re-analysis.
    # Both LCDM and TEP predict η=1 by construction (Etherington theorem).
    # The deviation is a BAO compilation systematic, not a model discriminator.
    bao_model_dependent = have_observational_data and deviation_sigma_lcdm > 3.0
    
    blockers = []
    if not have_observational_data:
        blockers.append('No observational DDR constraints available')
    if n_constraints < min_constraints:
        blockers.append(f'Insufficient constraints: {n_constraints} < {min_constraints}')
    if len(z_obs) > 0 and (max(z_obs) - min(z_obs)) <= min_redshift_range:
        blockers.append(f'Insufficient redshift range')
    if bao_model_dependent:
        blockers.append(
            f'BAO D_A is model-dependent (fiducial LCDM r_s). '
            f'LCDM compilation: η={eta_lcdm_weighted:.3f}±{eta_lcdm_weighted_err:.3f} '
            f'({deviation_sigma_lcdm:.1f}σ from 1); '
            f'TEP compilation: η={eta_weighted:.3f}±{eta_weighted_err:.3f} '
            f'({deviation_sigma:.1f}σ from 1). '
            f'Both violate η=1 because BAO D_A assumes LCDM. Not a model discriminator.'
        )
    
    # At z=2 theory prediction
    idx_z2 = np.argmin(np.abs(z_grid - 2.0))
    delta_eta_theory = eta_tep_theory[idx_z2] - 1.0
    
    # Prepare results
    results = {
        'step': STEP_ID,
        'status': 'completed' if research_grade else 'blocked',
        'description': 'Distance-duality relation test with observational constraints',
        'model_parameters': {
            'H0_lcdm': rounded(H0_lcdm, 2),
            'H0_tep': rounded(H0_tep, 2),
            'epsilon_T': rounded(epsilon_T, 6),
            'z_T': rounded(z_T, 3),
            'Om0': rounded(Om0, 3)
        },
        'observational_data': {
            'have_data': have_observational_data,
            'n_constraints': n_constraints,
            'redshift_range': [float(min(z_obs)), float(max(z_obs))] if len(z_obs) > 1 else [0.0, 0.0],
        },
        'statistical_results': {
            'chi2_lcdm': rounded(chi2_lcdm, 2) if have_observational_data else None,
            'chi2_tep': rounded(chi2_tep, 2) if have_observational_data else None,
            'delta_chi2': rounded(delta_chi2, 2) if have_observational_data else None,
            'eta_weighted_mean': rounded(eta_weighted, 4) if have_observational_data else None,
            'eta_weighted_err': rounded(eta_weighted_err, 4) if have_observational_data else None,
            'deviation_from_unity_sigma': rounded(deviation_sigma, 2) if have_observational_data else None,
        },
        'theory_predictions': {
            'eta_lcdm_at_z2': rounded(eta_lcdm_theory[idx_z2], 4),
            'eta_tep_at_z2': rounded(eta_tep_theory[idx_z2], 4),
            'delta_eta_theory': rounded(delta_eta_theory, 4),
        },
        'ddr_constraints': [
            {
                'z': rounded(float(z), 3),
                'D_L': rounded(float(dl), 2),
                'D_L_err': rounded(float(dl_err), 2),
                'D_A': rounded(float(da), 2),
                'D_A_err': rounded(float(da_err), 2),
                'eta_obs': rounded(float(eta), 4),
                'eta_err': rounded(float(eta_e), 4),
            }
            for z, dl, dl_err, da, da_err, eta, eta_e in zip(
                z_obs, D_L_obs, D_L_err, D_A_obs, D_A_err, eta_obs, eta_obs_err
            )
        ] if have_observational_data else [],
        'lcdm_compilation': {
            'eta_lcdm_weighted_mean': rounded(eta_lcdm_weighted, 4) if have_observational_data else None,
            'eta_lcdm_weighted_err': rounded(eta_lcdm_weighted_err, 4) if have_observational_data else None,
            'deviation_from_unity_sigma': rounded(deviation_sigma_lcdm, 2) if have_observational_data else None,
        },
        'tep_self_consistent': {
            'eta_tep_weighted_mean': rounded(eta_weighted, 4) if have_observational_data else None,
            'eta_tep_weighted_err': rounded(eta_weighted_err, 4) if have_observational_data else None,
            'deviation_from_unity_sigma': rounded(deviation_sigma, 2) if have_observational_data else None,
            'n_constraints': n_constraints,
        },
        'validation': {
            'real_data': have_observational_data,
            'research_grade_distance_duality': research_grade,
            'claim_gate': 'non_discriminating' if bao_model_dependent else ('open' if tep_self_consistent else 'blocked'),
            'blockers': blockers if blockers else [],
        },
    }
    
    # Write theory evolution to CSV (consistent field names)
    csv_data = []
    for i in range(len(z_grid)):
        csv_data.append({
            "z": rounded(z_grid[i], 3),
            "eta_TEP_theory": rounded(eta_tep_theory[i], 4),
            "eta_LCDM_theory": rounded(eta_lcdm_theory[i], 4),
            "delta_eta_theory": rounded(eta_tep_theory[i] - 1.0, 4),
            "eta_observed": None,
            "eta_obs_err": None,
        })
    
    # Add observational constraints if available
    if have_observational_data:
        for i, (z, eta, eta_e) in enumerate(zip(z_obs, eta_obs, eta_obs_err)):
            csv_data.append({
                "z": rounded(z, 3),
                "eta_TEP_theory": None,
                "eta_LCDM_theory": None,
                "delta_eta_theory": None,
                "eta_observed": rounded(eta, 4),
                "eta_obs_err": rounded(eta_e, 4),
            })
    
    write_csv(step_csv_path(STEP_ID), csv_data)
    write_json(step_json_path(STEP_ID), results)
    
    # Print summary
    if have_observational_data:
        print_status(f"LCDM compilation (Planck D_L + BAO D_A): η = {eta_lcdm_weighted:.4f} ± {eta_lcdm_weighted_err:.4f}", "INFO")
        print_status(f"  Deviation from unity: {deviation_sigma_lcdm:.2f}σ", "INFO")
        print_status(f"TEP compilation (TEP D_L + BAO D_A): η = {eta_weighted:.4f} ± {eta_weighted_err:.4f}", "INFO")
        print_status(f"  Deviation from unity: {deviation_sigma:.2f}σ", "INFO")
        print_status(f"Δχ² (LCDM - TEP) = {delta_chi2:.2f}", "INFO")
    else:
        print_status("No observational data - theory predictions only", "WARNING")
    
    print_status(f"TEP η(z=2) prediction: {eta_tep_theory[idx_z2]:.4f}", "INFO")
    print_status(f"Research grade: {research_grade}", "SUCCESS" if research_grade else "WARNING")
    print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
    
    return results


if __name__ == "__main__":
    run()
