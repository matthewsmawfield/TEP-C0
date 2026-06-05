#!/usr/bin/env python3
"""Step 025d: High-redshift Distance Duality Relation Test.

Tests the TEP prediction that η(z) → 1 at high redshift (z > 2).
The TEP screening model predicts convergence to unity as we probe
the unscreened high-z regime.

Uses Lyman-α BAO at z ~ 2.3-2.5 paired with high-z SNe Ia.

TEP Prediction: η(z=2.5) ≈ 0.9-0.95 (close to unity)
ΛCDM Prediction: η(z) ≡ 1 (exact)

Tier: RESEARCH GRADE

Requirements:
  - DESI/eBOSS Lyman-α BAO data (step_023f)
  - High-z SNe Ia from Pantheon+ (z > 2)
  - TEP screening model parameters

References:
  - TEP screening model: η_probe(z) = η₀ + (1-η₀) × (z/(z+z_c))^α
  - Step 025c TEP screening fit results
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from c0_common import TEPLogger, set_step_logger, ensure_dirs, print_status, step_json_path, write_json

STEP_ID = "step_04_07_highz_ddr"


def load_highz_constraints() -> List[Dict[str, Any]]:
    """Load high-z BAO constraints from step_023f."""
    constraints = []
    
    # Load from step_023f output
    json_path = Path("results/step_01_07_download_desi.json")
    if json_path.exists():
        with open(json_path) as f:
            data = json.load(f)
            constraints = data.get("ddr_constraints", [])
            
    return constraints


def compute_eta(
    D_L: float,
    D_A: float,
    z: float,
    D_L_err: float = 0.0,
    D_A_err: float = 0.0,
) -> Tuple[float, float]:
    """Compute η = D_L / (D_A × (1+z)²).
    
    Args:
        D_L: Luminosity distance
        D_A: Angular diameter distance
        z: Redshift
        D_L_err: Error on D_L
        D_A_err: Error on D_A
        
    Returns:
        (η, η_err)
    """
    D_A_safe = max(D_A, np.finfo(float).tiny)
    D_L_safe = max(D_L, np.finfo(float).tiny)
    eta = D_L_safe / (D_A_safe * (1 + z)**2)
    
    # Error propagation
    if D_L_err > 0 and D_A_err > 0:
        rel_err_L = D_L_err / D_L_safe
        rel_err_A = D_A_err / D_A_safe
        eta_err = eta * np.sqrt(rel_err_L**2 + rel_err_A**2)
    else:
        eta_err = 1e10  # Large error if not specified
        
    return eta, eta_err


def tep_screening_prediction(
    z: float,
    probe: str = "BAO",
) -> float:
    """Compute TEP screening model prediction for η(z).
    
    Uses best-fit parameters from step_025c.
    
    Formula: η(z) = η₀ + (1-η₀) × (z/(z+z_c))^α
    
    Args:
        z: Redshift
        probe: Probe type (BAO, SZ, SGL)
        
    Returns:
        Predicted η value
    """
    # Best-fit parameters from step_025c
    params = {
        "BAO": {"eta_0": 0.05, "z_c": 0.038, "alpha": 5.0},
        "SGL": {"eta_0": 0.284, "z_c": 0.168, "alpha": 3.0},
        "SZ": {"eta_0": 0.359, "z_c": 3.0, "alpha": 3.0},
    }
    
    p = params.get(probe, params["BAO"])
    eta_0 = p["eta_0"]
    z_c = p["z_c"]
    alpha = p["alpha"]
    
    # TEP screening formula
    eta = eta_0 + (1 - eta_0) * (z / (z + z_c))**alpha
    
    return eta


def test_eta_convergence(
    z_values: List[float],
    eta_values: List[float],
    eta_errors: List[float],
) -> Dict[str, Any]:
    """Test if η → 1 at high redshift.
    
    Args:
        z_values: Redshifts
        eta_values: Measured η values
        eta_errors: Errors on η
        
    Returns:
        Dictionary with test results
    """
    results = {
        "tep_prediction": {},
        "lcdm_prediction": {},
        "test_results": {},
    }
    
    for z, eta, eta_err in zip(z_values, eta_values, eta_errors):
        # TEP prediction
        eta_tep = tep_screening_prediction(z, "BAO")
        
        # ΛCDM prediction
        eta_lcdm = 1.0
        
        # Deviation from unity
        deviation = eta - 1.0
        sigma = abs(deviation) / eta_err if eta_err > 0 else 0
        
        # Deviation from TEP
        deviation_tep = eta - eta_tep
        sigma_tep = abs(deviation_tep) / eta_err if eta_err > 0 else 0
        
        results["tep_prediction"][f"z_{z:.2f}"] = {
            "z": z,
            "eta_tep": eta_tep,
            "eta_observed": eta,
            "eta_err": eta_err,
            "deviation": deviation_tep,
            "sigma": sigma_tep,
        }
        
        results["lcdm_prediction"][f"z_{z:.2f}"] = {
            "z": z,
            "eta_lcdm": eta_lcdm,
            "eta_observed": eta,
            "eta_err": eta_err,
            "deviation": deviation,
            "sigma": sigma,
        }
        
    # Overall test
    if len(z_values) > 0:
        # Check if η increases with z (convergence to unity)
        if len(z_values) >= 2:
            z_sorted = np.argsort(z_values)
            eta_sorted = [eta_values[i] for i in z_sorted]
            increasing = all(eta_sorted[i] <= eta_sorted[i+1] for i in range(len(eta_sorted)-1))
        else:
            increasing = None
            
        # Check if η is consistent with unity at highest z
        max_z_idx = np.argmax(z_values)
        eta_at_highz = eta_values[max_z_idx]
        eta_err_at_highz = eta_errors[max_z_idx]
        
        # Within 2σ of unity?
        consistent_with_unity = abs(eta_at_highz - 1.0) < 2 * eta_err_at_highz
        
        results["test_results"] = {
            "n_constraints": len(z_values),
            "z_range": [min(z_values), max(z_values)] if z_values else [0, 0],
            "eta_increasing_with_z": increasing,
            "consistent_with_unity_at_highz": consistent_with_unity,
            "eta_at_max_z": float(eta_at_highz),
            "eta_err_at_max_z": float(eta_err_at_highz),
            "deviation_from_unity": float(eta_at_highz - 1.0),
            "sigma_from_unity": float(abs(eta_at_highz - 1.0) / eta_err_at_highz) if eta_err_at_highz > 0 else 0,
        }
        
    return results


def compute_model_dl(z: float, h0: float = 70.0, om0: float = 0.3) -> float:
    """Compute model D_L using proper FLRW integration."""
    # Note: CosmologyFLRW includes radiation component (Or0 computed from CMB temperature)
    from core.cosmology import CosmologyFLRW
    cosmo_lcdm = CosmologyFLRW(H0=h0, Om0=om0, Ok0=0.0)  # Explicitly flat universe
    z_arr = np.array([z])
    return float(cosmo_lcdm.luminosity_distance(z_arr)[0])


def run() -> dict:
    """Run high-z DDR test."""
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()
    
    # Load high-z constraints
    constraints = load_highz_constraints()
    print_status(f"Loaded {len(constraints)} high-z constraints", "INFO")
    
    # Load TEP parameters from step_022
    try:
        step022 = json.loads(Path("results/step_03_01_three_model_comparison.json").read_text())
        m1_key = "M1_free_zT"
        if m1_key not in step022.get("models", {}):
            m1_key = "M1_NoLambda_zT5"
        m1_params = step022['models'][m1_key]['parameters_mle']
        h0_tep = 70.0  # dimensionless model fixes H0_ref
        om0_tep = 1.0  # M1_NoLambda is matter-only (no Lambda)
    except (KeyError, FileNotFoundError, json.JSONDecodeError):
        h0_tep = 70.0
        om0_tep = 0.3
    
    # Compute η for each constraint
    eta_results = []
    z_values = []
    eta_values = []
    eta_errors = []
    
    for constraint in constraints:
        z = constraint.get("z", 0)
        D_A = constraint.get("D_A", 0)
        D_A_err = constraint.get("D_A_err", 0)
        D_L_obs = constraint.get("D_L")
        
        # If no D_L from SNe, use model prediction
        if D_L_obs is None:
            # Use TEP model to predict D_L from D_A
            # For flat universe: D_L = D_A * (1+z)^2 / η
            # But we don't know η yet... use iterative approach
            # Start with ΛCDM assumption η=1, then apply TEP correction
            D_L_lcdm = D_A * (1+z)**2  # ΛCDM prediction
            # TEP predicts η < 1 at low z, so D_L_tep = η * D_A * (1+z)^2 < D_L_lcdm
            eta_tep = tep_screening_prediction(z, "BAO")
            D_L_model = D_L_lcdm * eta_tep  # TEP-modified D_L
            D_L_obs = D_L_model
            D_L_err = D_L_model * 0.1  # 10% systematic uncertainty from model
            model_dependent = True
            print_status(f"Using model D_L for z={z:.2f}: D_L={D_L_model:.0f} Mpc (η_TEP={eta_tep:.3f})", "INFO")
        else:
            D_L_err = constraint.get("D_L_err", D_L_obs * 0.05)
            model_dependent = False
        
        # Compute η
        eta, eta_err = compute_eta(D_L_obs, D_A, z, D_L_err, D_A_err)
        
        # TEP prediction
        eta_tep = tep_screening_prediction(z, "BAO")
        
        result = {
            "z": z,
            "D_A": D_A,
            "D_A_err": D_A_err,
            "D_L": D_L_obs,
            "D_L_err": D_L_err,
            "D_L_model_dependent": model_dependent,
            "eta": float(eta),
            "eta_err": float(eta_err),
            "eta_tep_pred": float(eta_tep),
            "source": constraint.get("source", ""),
        }
        
        eta_results.append(result)
        z_values.append(z)
        eta_values.append(eta)
        eta_errors.append(eta_err)
        
        print_status(f"z={z:.3f}: η = {eta:.3f} ± {eta_err:.3f} (TEP pred: {eta_tep:.3f})", "INFO")
    
    # Test η → 1 convergence
    test_results = test_eta_convergence(z_values, eta_values, eta_errors)
    
    # Critical test: Does η approach unity at high z?
    tep_prediction_z2 = tep_screening_prediction(2.5, "BAO")
    
    print_status(f"TEP prediction at z=2.5: η = {tep_prediction_z2:.3f}", "INFO")
    print_status(f"ΛCDM prediction at z=2.5: η = 1.000", "INFO")
    
    if test_results.get("test_results"):
        tr = test_results["test_results"]
        print_status(f"η at max z ({tr['z_range'][1]:.2f}): {tr['eta_at_max_z']:.3f} ± {tr['eta_err_at_max_z']:.3f}", "INFO")
        print_status(f"Deviation from unity: {tr['deviation_from_unity']:.3f}σ", "INFO")
        print_status(f"Consistent with unity: {tr['consistent_with_unity_at_highz']}", "INFO")
    
    # Determine status - if model-dependent D_L used, mark appropriately
    has_model_dependent = any(r.get("D_L_model_dependent", False) for r in eta_results)
    status = "completed" if len(eta_results) > 0 else "insufficient_data"
    if has_model_dependent:
        status = "model_dependent"
    
    # Prepare output
    payload = {
        "step": STEP_ID,
        "description": "High-z DDR test of TEP η → 1 prediction",
        "status": status,
        "n_constraints": len(eta_results),
        "n_model_dependent": sum(1 for r in eta_results if r.get("D_L_model_dependent", False)),
        "eta_results": eta_results,
        "test_results": test_results,
        "predictions": {
            "tep_eta_at_z2.5": float(tep_prediction_z2),
            "lcdm_eta_at_z2.5": 1.0,
            "tep_convergence_z": "As z → ∞, η → 1 for all probes",
        },
        "validation": {
            "sufficient_constraints": len(eta_results) >= 2,
            "z_max": max(z_values) if z_values else 0,
            "has_model_dependent_dl": has_model_dependent,
            "test_passed": test_results.get("test_results", {}).get("consistent_with_unity_at_highz", False) if test_results.get("test_results") else False,
        },
    }
    
    write_json(step_json_path(STEP_ID), payload)
    print_status(f"Step {STEP_ID} completed", "SUCCESS")
    
    return payload


if __name__ == "__main__":
    run()
