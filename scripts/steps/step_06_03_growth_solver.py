#!/usr/bin/env python3
"""Step 015: Structure Growth Solver - CLASS-based TEP implementation.

Uses TEP-CLASS v2.0 to compute proper matter power spectrum and growth functions.

Note: CLASS parameters include explicit curvature (Omega_k=0.0 for flat universe)
and radiation component (N_ur=2.0328 for neutrinos, plus photons from CMB temperature).
"""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

# Add TEP-CLASS to path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
TEP_CLASS_PATH = PROJECT_ROOT / "external" / "class"
import glob
_build_dir = PROJECT_ROOT / "external" / "class" / "build"
_lib_dirs = glob.glob(str(_build_dir / "lib.*"))
if _lib_dirs:
    CLASS_BUILD_PATH_STR = str(_lib_dirs[0])
else:
    CLASS_BUILD_PATH_STR = str(_build_dir / "lib.macosx-11.1-arm64-cpython-313")
TEP_CLASS_BUILD = Path(CLASS_BUILD_PATH_STR)
if str(TEP_CLASS_BUILD) not in sys.path:
    sys.path.insert(0, str(TEP_CLASS_BUILD))

from c0_common import (
    TEPLogger, ensure_dirs, print_status, read_json, 
    rounded, set_step_logger, step_json_path, step_csv_path, 
    write_csv, write_json
)

STEP_ID = "step_06_03_growth_solver"

def compute_growth_class(Om0, epsilon_T=0.1, z_T=3.0, use_tep=False):
    """Compute growth function using CLASS/TEP-CLASS.
    
    Returns:
        dict with growth factor D(z), growth rate f(z), sigma_8
    """
    try:
        from classy import Class
    except ImportError:
        print_status("TEP-CLASS not available, using approximation", "WARNING")
        return None
    
    params = {
        'output': 'mPk',
        'P_k_max_h/Mpc': 10.0,
        'z_max_pk': 3.0,
        'A_s': 2.1e-9,
        'n_s': 0.966,
        'h': 0.7,
        'omega_b': 0.0224,
        'omega_cdm': Om0 * 0.7**2 - 0.0224,
        'Omega_k': 0.0,
        'N_ur': 2.0328,
        'N_ncdm': 1,
        'm_ncdm': 0.06,
    }
    
    if use_tep:
        params['tep_epsilon_T'] = epsilon_T
        params['tep_z_T'] = z_T
        params['tep_n_T'] = 1.0
    
    cosmo = Class()
    cosmo.set(params)
    
    try:
        cosmo.compute()
    except Exception as e:
        print_status(f"CLASS computation failed: {e}", "WARNING")
        return None
    
    # Compute growth at various redshifts
    z_grid = np.linspace(0, 3, 30)
    growth_D = []
    growth_f = []
    
    for z in z_grid:
        try:
            # Get sigma_8 at this redshift
            sigma8_z = cosmo.sigma(8.0 / 0.7, z)
            growth_D.append(sigma8_z)
            
            # Approximate growth rate f ~ -d ln D / d ln a
            if z > 0:
                dz = 0.01
                sigma8_zp = cosmo.sigma(8.0 / 0.7, z + dz)
                f = (np.log(sigma8_zp) - np.log(sigma8_z)) / (np.log(1/(1+z+dz)) - np.log(1/(1+z)))
                growth_f.append(f)
            else:
                growth_f.append(0.55)  # Approximation at z=0
        except (ValueError, IndexError, KeyError, AttributeError) as e:
            growth_D.append(0.0)
            growth_f.append(0.0)
    
    cosmo.struct_cleanup()
    cosmo.empty()
    
    return {
        'z': z_grid,
        'D': np.array(growth_D),
        'f': np.array(growth_f),
        'sigma_8': growth_D[0] if len(growth_D) > 0 else 0.811
    }

def run():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    step022 = read_json(step_json_path("step_03_01_three_model_comparison"))
    m1_key = "M1_free_zT" if "M1_free_zT" in step022.get("models", {}) else "M1_NoLambda_zT5"
    m1 = step022['models'][m1_key]['parameters_mle']
    Om0 = 1.0  # M1_NoLambda is matter-only (no Lambda); Om0 is fixed, not fitted
    epsilon_T = m1.get('epsilon_T', 0.1)

    # Try to use CLASS for proper growth calculation
    class_available = False
    try:
        lcdm_results = compute_growth_class(0.3, use_tep=False)
        tep_results = compute_growth_class(Om0, epsilon_T=epsilon_T, use_tep=True)
        class_available = lcdm_results is not None and tep_results is not None
    except Exception as e:
        print_status(f"CLASS computation failed: {e}", "WARNING")
        lcdm_results = None
        tep_results = None

    if class_available and lcdm_results is not None and tep_results is not None:
        # Use CLASS results
        z_grid = lcdm_results['z']
        D_lcdm = lcdm_results['D']
        D_tep = tep_results['D']
        f_lcdm = lcdm_results['f']
        f_tep = tep_results['f']
        sigma_8_lcdm = lcdm_results['sigma_8']
        sigma_8_tep = tep_results['sigma_8']
        research_grade = True
        blockers = []
        claim_gate = 'open'
        forecast_only = False
    else:
        # Fallback to Carroll+Press 1992 approximation for ΛCDM growth
        print_status("Using fallback approximation for growth", "WARNING")
        z_grid = np.linspace(0, 3, 100)
        a_grid = 1.0 / (1.0 + z_grid)
        Om0_ref = 0.3
        # Growth factor: D(a) ∝ a * g(a) where g(a) ~ Omega_m(a)^0.55 / Omega_m0^0.55
        # Normalized to D(a=1)=1
        Omega_m_a = Om0_ref / (Om0_ref + (1-Om0_ref) * a_grid**3)
        D_lcdm = a_grid * (Omega_m_a / Om0_ref)**0.55
        D_tep = a_grid * (Omega_m_a / Om0)**0.55  # TEP with different Om0
        f_lcdm = np.full_like(z_grid, 0.55)
        f_tep = np.full_like(z_grid, 0.55)
        sigma_8_lcdm = 0.811
        sigma_8_tep = 0.811 * (Om0/0.3)**0.5
        research_grade = False
        blockers = ['TEP-CLASS growth computation unavailable or failed']
        claim_gate = 'blocked'
        forecast_only = True

    results = {
        'step': STEP_ID,
        'metrics': {
            'Ok0': 0.0,
            'sigma_8_tep': rounded(sigma_8_tep, 3),
            'sigma_8_lcdm': rounded(sigma_8_lcdm, 3),
            'growth_at_z1_lcdm': rounded(np.interp(1.0, z_grid, D_lcdm), 3),
            'growth_at_z1_tep': rounded(np.interp(1.0, z_grid, D_tep), 3),
            'f_sigma8_z1_lcdm': rounded(np.interp(1.0, z_grid, f_lcdm) * sigma_8_lcdm, 3),
            'f_sigma8_z1_tep': rounded(np.interp(1.0, z_grid, f_tep) * sigma_8_tep, 3),
            'class_used': class_available
        },
        'validation': {
            'forecast_only': forecast_only,
            'research_grade_growth': research_grade,
            'claim_gate': claim_gate,
            'blockers': blockers,
        },
    }
    
    write_json(step_json_path(STEP_ID), results)
    
    csv_rows = []
    for i in range(len(z_grid)):
        row = {
            "z": rounded(z_grid[i], 3),
            "D_lcdm": rounded(D_lcdm[i], 4),
            "D_tep": rounded(D_tep[i], 4)
        }
        if class_available:
            row["f_lcdm"] = rounded(f_lcdm[i], 4)
            row["f_tep"] = rounded(f_tep[i], 4)
        csv_rows.append(row)
    write_csv(step_csv_path(STEP_ID), csv_rows)
    
    print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
    return results

if __name__ == "__main__":
    run()
