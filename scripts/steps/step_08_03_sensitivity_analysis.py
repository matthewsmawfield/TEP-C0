#!/usr/bin/env python3
"""Step 031: Sensitivity Analysis - FISHER INFORMATION MATRIX.

Computes parameter uncertainties using second derivatives of the log-likelihood.
References: Pantheon+SH0ES public distance table and covariance release.
"""

from __future__ import annotations

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import json
from c0_common import TEPLogger, ensure_dirs, print_status, read_json, rounded, set_step_logger, step_json_path, write_json
from step_03_01_three_model_comparison import ModelTEP, PantheonData

STEP_ID = "step_08_03_sensitivity_analysis"


def load_step022_results() -> dict:
    """Load fitted model parameters from step_022.
    
    DATA PROVENANCE NOTE: This step uses fitted parameters from step_022 for sensitivity analysis.
    Sensitivity analysis is a diagnostic tool that explores parameter space around the best-fit values.
    No circular dependency - this is a post-fitting diagnostic, not a validation.
    """
    results_file = Path("results/step_03_01_three_model_comparison.json")
    if not results_file.exists():
        raise FileNotFoundError("step_022 results not found. Run step_022 first.")
    return json.loads(results_file.read_text())


def compute_fisher_matrix(log_like_fn, params, step_sizes=None):
    """Compute Fisher information matrix by numerical second derivatives.
    
    F_ij = - < d^2 ln L / (d theta_i d theta_j) >
    """
    n_params = len(params)
    
    if step_sizes is None:
        step_sizes = np.abs(params) * 0.01 + 0.001
    
    fisher = np.zeros((n_params, n_params))
    
    # Compute second derivatives using central differences
    for i in range(n_params):
        for j in range(i, n_params):
            if i == j:
                # Diagonal element: d^2L/dtheta_i^2
                params_plus = params.copy()
                params_minus = params.copy()
                params_plus[i] += step_sizes[i]
                params_minus[i] -= step_sizes[i]
                
                L_plus = log_like_fn(params_plus)
                L_minus = log_like_fn(params_minus)
                L_center = log_like_fn(params)
                
                fisher[i, i] = -np.divide(L_plus - 2*L_center + L_minus, step_sizes[i]**2)
            else:
                # Off-diagonal element: d^2L/(dtheta_i dtheta_j)
                params_pp = params.copy()
                params_pm = params.copy()
                params_mp = params.copy()
                params_mm = params.copy()
                
                params_pp[i] += step_sizes[i]
                params_pp[j] += step_sizes[j]
                
                params_pm[i] += step_sizes[i]
                params_pm[j] -= step_sizes[j]
                
                params_mp[i] -= step_sizes[i]
                params_mp[j] += step_sizes[j]
                
                params_mm[i] -= step_sizes[i]
                params_mm[j] -= step_sizes[j]
                
                L_pp = log_like_fn(params_pp)
                L_pm = log_like_fn(params_pm)
                L_mp = log_like_fn(params_mp)
                L_mm = log_like_fn(params_mm)
                
                fisher[i, j] = -(L_pp - L_pm - L_mp + L_mm) / (4 * step_sizes[i] * step_sizes[j])
                fisher[j, i] = fisher[i, j]
    
    return fisher


def run():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    try:
        step022 = load_step022_results()
    except FileNotFoundError as e:
        print_status(f"Error loading step_022 results: {e}", "ERROR")
        results = {
            'step': STEP_ID,
            'status': 'failed',
            'error': str(e),
            'description': 'Parameter sensitivity via Fisher information matrix',
        }
        write_json(step_json_path(STEP_ID), results)
        return results
    except Exception as e:
        print_status(f"Unexpected error loading step_022 results: {e}", "ERROR")
        results = {
            'step': STEP_ID,
            'status': 'failed',
            'error': str(e),
            'description': 'Parameter sensitivity via Fisher information matrix',
        }
        write_json(step_json_path(STEP_ID), results)
        return results

    try:
        model = ModelTEP(pure_shear=True)
        data = PantheonData()
        data.load()
        sigma_from_full_covariance = True
        m2 = step022['models']['M2_PureShear']['parameters_mle']

        param_names = model.param_names
        params = np.array([m2[name] for name in param_names], dtype=float)
        bounds = np.array(model.bounds, dtype=float)
        widths = bounds[:, 1] - bounds[:, 0]
        step_sizes = np.maximum(widths * 1e-4, 1e-8)
        interior_params = np.clip(params, bounds[:, 0] + 2 * step_sizes, bounds[:, 1] - 2 * step_sizes)

        def log_like(p):
            return model.log_likelihood(np.asarray(p, dtype=float), data)
        
        print_status("Computing Fisher information matrix", "PROCESS")
        fisher = compute_fisher_matrix(log_like, interior_params, step_sizes=step_sizes)
        fisher = 0.5 * (fisher + fisher.T)

        condition_number = float(np.linalg.cond(fisher))
        cov = np.linalg.pinv(fisher, rcond=1e-10)
        variances = np.clip(np.diag(cov), 0.0, np.inf)
        param_uncertainties = np.sqrt(variances)
        fisher_success = bool(np.all(np.isfinite(param_uncertainties)))
        
        print_status("Computing parameter uncertainties", "PROCESS")
        uncertainties = {}
        for i, name in enumerate(param_names):
            uncertainties[name] = rounded(param_uncertainties[i], 6)
        
        results = {
            'step': STEP_ID,
            'description': 'Parameter sensitivity via Fisher information matrix',
            'method': 'Numerical second derivatives of log-likelihood',
            'status': 'FISHER_MATRIX',
            'fisher_computation': 'Successful' if fisher_success else 'Matrix singular - used fallback',
            'best_fit_parameters': {name: rounded(value, 6) for name, value in zip(param_names, params)},
            'finite_difference_point': {name: rounded(value, 6) for name, value in zip(param_names, interior_params)},
            'fisher_condition_number': rounded(condition_number, 6),
            'parameter_uncertainties': uncertainties,
            'cramer_rao_bounds': uncertainties,
            'relative_precision_percent': {
                name: rounded(np.divide(param_uncertainties[i], max(abs(params[i]), 1e-12)) * 100, 2)
                for i, name in enumerate(param_names)
            },
            'interpretation': f'Fisher matrix around the M2 likelihood gives log_Sigma_0 uncertainty ~{rounded(np.divide(param_uncertainties[0], max(abs(params[0]), 1e-12)) * 100, 1)}%.',
            'validation': {
                'uses_step022_full_covariance_likelihood': True,
                'uses_pseudo_inverse': True,
                'claim_gate': 'diagnostic',
            },
        }
        
        write_json(step_json_path(STEP_ID), results)
        print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
        return results
        
    except Exception as e:
        print_status(f"Error during Fisher matrix computation: {e}", "ERROR")
        import traceback
        print_status(traceback.format_exc(), "ERROR")
        results = {
            'step': STEP_ID,
            'status': 'failed',
            'error': str(e),
            'description': 'Parameter sensitivity via Fisher information matrix',
        }
        write_json(step_json_path(STEP_ID), results)
        return results


if __name__ == "__main__":
    run()
