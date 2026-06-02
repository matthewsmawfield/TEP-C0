#!/usr/bin/env python3
"""
Diagnostic: test whether epsilon_T wants to exceed the [0, 1] prior bound.
Computes likelihood for epsilon_T = [0.5, 0.7, 0.9, 1.0, 1.2, 1.5, 2.0]
using the step_03_01 data and model setup.
"""
import sys
from pathlib import Path
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "steps"))

from core.tep_cosmology import TEPCosmology
from core.cosmology import CosmologyFLRW
from c0_common import read_json

def load_step_03_data():
    """Load the Pantheon+ data and covariance from step_03_01 results."""
    # Read the step_03_01 JSON for provenance info
    step022 = read_json(PROJECT_ROOT / "results" / "step_03_01_three_model_comparison.json")
    
    # Load actual data
    import pandas as pd
    df = pd.read_csv(PROJECT_ROOT / "data" / "raw" / "pantheon_plus_shoes.dat", sep=r'\s+', comment='#')
    z = df['zCMB'].values
    mb = df['m_b_corr'].values if 'm_b_corr' in df.columns else df['mB'].values
    
    n_raw = len(z)
    # Filter same as step_03_01
    valid = np.isfinite(z) & np.isfinite(mb) & (z > 0.001)
    z = z[valid]
    mb = mb[valid]
    
    # Load covariance (full) - .cov file has special format
    cov_flat = np.loadtxt(PROJECT_ROOT / "data" / "raw" / "Pantheon+SH0ES.cov")
    if cov_flat.size == n_raw * n_raw + 1:
        cov_full = cov_flat[1:].reshape(n_raw, n_raw)
    elif cov_flat.size == n_raw * n_raw:
        cov_full = cov_flat.reshape(n_raw, n_raw)
    else:
        raise ValueError(f"Covariance size {cov_flat.size} mismatch with data size {n_raw}")
    # Apply same mask
    cov = cov_full[np.ix_(valid, valid)]
    
    return z, mb, cov

def compute_loglike(z, mb, cov, epsilon_T, z_T=5.0, H0_ref=70.0):
    """Compute Gaussian log-likelihood for a given epsilon_T (fixed z_T)."""
    from scipy.linalg import cholesky
    
    tep = TEPCosmology(H0=H0_ref, Omega_m=1.0, epsilon_T=epsilon_T, z_T=z_T)
    mu_tep = tep.distance_modulus(z)
    
    # Nuisance M: fit by minimizing chi2 w.r.t. M
    # chi2 = (mb - mu_tep - M)^T C^{-1} (mb - mu_tep - M)
    # d(chi2)/dM = 0 => M = sum_i (C^{-1} (mb - mu_tep))_i / sum_i (C^{-1} 1)_i
    C_inv = np.linalg.inv(cov)
    residuals = mb - mu_tep
    ones = np.ones(len(z))
    M_opt = float(np.dot(ones, np.dot(C_inv, residuals)) / np.dot(ones, np.dot(C_inv, ones)))
    
    residuals_best = residuals - M_opt
    chi2 = float(np.dot(residuals_best, np.dot(C_inv, residuals_best)))
    
    # Log-likelihood (ignoring constant)
    sign, logdet = np.linalg.slogdet(cov)
    loglike = -0.5 * chi2 - 0.5 * logdet
    
    return {
        'epsilon_T': epsilon_T,
        'M': M_opt,
        'chi2': chi2,
        'loglike': loglike,
        'dof': len(z) - 2,  # epsilon_T and M
    }

def main():
    print("=" * 70)
    print("epsilon_T Prior Bound Diagnostic")
    print("=" * 70)
    
    z, mb, cov = load_step_03_data()
    print(f"Data: {len(z)} SNe, covariance shape: {cov.shape}")
    
    # Test epsilon_T values including beyond the current [0, 1] bound
    test_values = [0.0, 0.3, 0.5, 0.7, 0.9, 1.0, 1.2, 1.5, 2.0]
    
    print()
    print(f"{'epsilon_T':>10} {'M':>12} {'chi2':>12} {'loglike':>14} {'chi2/dof':>10}")
    print("-" * 70)
    
    results = []
    for eps in test_values:
        r = compute_loglike(z, mb, cov, eps)
        results.append(r)
        print(f"{eps:10.2f} {r['M']:12.4f} {r['chi2']:12.2f} {r['loglike']:14.4f} {r['chi2']/r['dof']:10.4f}")
    
    # Find maximum likelihood
    best = max(results, key=lambda x: x['loglike'])
    print()
    print(f"Maximum likelihood at epsilon_T = {best['epsilon_T']:.2f}")
    print(f"  M = {best['M']:.4f}")
    print(f"  chi2 = {best['chi2']:.2f}")
    print(f"  loglike = {best['loglike']:.4f}")
    
    # Compare to current bound [0, 1]
    at_1 = next(r for r in results if r['epsilon_T'] == 1.0)
    print()
    print(f"At epsilon_T = 1.0 (current prior upper bound):")
    print(f"  loglike = {at_1['loglike']:.4f}")
    print(f"  Delta loglike vs best = {at_1['loglike'] - best['loglike']:.4f}")
    
    # Check if likelihood increases monotonically toward the boundary
    increasing = all(results[i]['loglike'] < results[i+1]['loglike'] 
                     for i in range(len(results)-1) 
                     if results[i]['epsilon_T'] <= 1.0)
    print()
    if increasing:
        print("WARNING: Likelihood increases monotonically toward epsilon_T = 1.0")
        print("         The prior upper bound is truncating the posterior.")
        print("         Consider extending the prior to epsilon_T > 1.0.")
    else:
        print("Likelihood does NOT increase monotonically toward the bound.")
        print("The [0, 1] prior is likely sufficient.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
