#!/usr/bin/env python3
"""
Diagnostic: test unscreened TEP (z_T = inf) to see maximum possible evidence.
In this limit, gamma(z) = 1 + epsilon_T * ln(1+z) for all z.
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

# Monkey-patch TEPCosmology to disable screening
def tep_gamma_no_screening(self, z):
    """Unscreened: gamma(z) = 1 + epsilon_T * ln(1+z) everywhere."""
    z_arr = np.asarray(z, dtype=float)
    if self.epsilon_T == 0:
        return np.ones_like(z_arr)
    log_factor = np.log(1.0 + z_arr)
    enhancement = 1.0 + self.epsilon_T * log_factor
    return np.maximum(enhancement, 0.1)

TEPCosmology.tep_gamma = tep_gamma_no_screening


def load_data_and_cov():
    import pandas as pd
    df = pd.read_csv(PROJECT_ROOT / "data" / "raw" / "pantheon_plus_shoes.dat", sep=r'\s+', comment='#')
    z = df['zCMB'].values
    mb = df['m_b_corr'].values if 'm_b_corr' in df.columns else df['mB'].values
    valid = np.isfinite(z) & np.isfinite(mb) & (z > 0.001)
    z = z[valid]
    mb = mb[valid]
    n_raw = len(valid)
    cov_flat = np.loadtxt(PROJECT_ROOT / "data" / "raw" / "Pantheon+SH0ES.cov")
    if cov_flat.size == n_raw * n_raw + 1:
        cov_full = cov_flat[1:].reshape(n_raw, n_raw)
    elif cov_flat.size == n_raw * n_raw:
        cov_full = cov_flat.reshape(n_raw, n_raw)
    else:
        raise ValueError(f"Covariance size {cov_flat.size} mismatch")
    cov = cov_full[np.ix_(valid, valid)]
    return z, mb, cov


def fit_epsilon_t(z, mb, cov, H0_ref=70.0, epsilon_T_range=None):
    """Grid search epsilon_T with analytical M marginalization."""
    C_inv = np.linalg.inv(cov)
    ones = np.ones(len(z))
    
    if epsilon_T_range is None:
        epsilon_T_range = np.linspace(0.0, 1.5, 61)
    
    best_loglike = -np.inf
    best_eps = 0.0
    best_M = 0.0
    best_chi2 = 0.0
    
    for eps in epsilon_T_range:
        tep = TEPCosmology(H0=H0_ref, Omega_m=1.0, epsilon_T=eps, z_T=1e9)
        mu_tep = tep.distance_modulus(z)
        residuals = mb - mu_tep
        
        # Analytical M = offset fit
        denom = np.dot(ones, np.dot(C_inv, ones))
        M_opt = float(np.dot(ones, np.dot(C_inv, residuals)) / denom)
        residuals_best = residuals - M_opt
        chi2 = float(np.dot(residuals_best, np.dot(C_inv, residuals_best)))
        
        sign, logdet = np.linalg.slogdet(cov)
        loglike = -0.5 * chi2 - 0.5 * logdet
        
        if loglike > best_loglike:
            best_loglike = loglike
            best_eps = eps
            best_M = M_opt
            best_chi2 = chi2
    
    return {
        'epsilon_T': best_eps,
        'M': best_M,
        'chi2': best_chi2,
        'loglike': best_loglike,
        'dof': len(z) - 2,
    }


def main():
    print("=" * 70)
    print("Unscreened TEP (z_T = inf) Diagnostic")
    print("=" * 70)
    
    z, mb, cov = load_data_and_cov()
    print(f"Data: {len(z)} SNe")
    
    # LCDM baseline
    lcdm = CosmologyFLRW(H0=70.0, Om0=0.3809)
    mu_lcdm = lcdm.distance_modulus(z)
    C_inv = np.linalg.inv(cov)
    ones = np.ones(len(z))
    
    residuals_lcdm = mb - mu_lcdm
    M_lcdm = float(np.dot(ones, np.dot(C_inv, residuals_lcdm)) / np.dot(ones, np.dot(C_inv, ones)))
    chi2_lcdm = float(np.dot(residuals_lcdm - M_lcdm, np.dot(C_inv, residuals_lcdm - M_lcdm)))
    
    # Fit unscreened TEP
    print("Fitting unscreened TEP on grid...")
    result = fit_epsilon_t(z, mb, cov, epsilon_T_range=np.linspace(0.0, 1.5, 151))
    
    n = len(z)
    k_tep = 2  # epsilon_T, M
    k_lcdm = 2  # Om0, M
    
    bic_tep = result['chi2'] + k_tep * np.log(n)
    bic_lcdm = chi2_lcdm + k_lcdm * np.log(n)
    delta_bic = bic_tep - bic_lcdm
    
    print()
    print(f"{'':30} {'LCDM':>15} {'Unscreened TEP':>15}")
    print("-" * 70)
    print(f"{'epsilon_T':30} {'0.000':>15} {result['epsilon_T']:>15.4f}")
    print(f"{'M':30} {M_lcdm:>15.4f} {result['M']:>15.4f}")
    print(f"{'chi2':30} {chi2_lcdm:>15.2f} {result['chi2']:>15.2f}")
    print(f"{'chi2/dof':30} {chi2_lcdm/(n-2):>15.4f} {result['chi2']/result['dof']:>15.4f}")
    print(f"{'BIC':30} {bic_lcdm:>15.2f} {bic_tep:>15.2f}")
    print(f"{'Delta BIC (TEP - LCDM)':30} {'':15} {delta_bic:>15.2f}")
    
    # Evidence interpretation
    ln_bf = -0.5 * delta_bic
    print()
    print(f"ln(BF) = {ln_bf:.3f}, BF = {np.exp(ln_bf):.2f}")
    if delta_bic < 0:
        print("Unscreened TEP is BIC-preferred over LCDM")
    else:
        print("LCDM is BIC-preferred")
    
    # Compare to screened z_T=5 result
    step022 = read_json(PROJECT_ROOT / "results" / "step_03_01_three_model_comparison.json")
    m1_zt5 = step022['models']['M1_NoLambda_zT5']
    print()
    print("Comparison with screened z_T=5 (from step_03_01):")
    print(f"  z_T=5:    epsilon_T = {m1_zt5['parameters_mle']['epsilon_T']:.4f}, "
          f"chi2 = {m1_zt5['chi2_mle']:.2f}, BIC = {m1_zt5['bic']:.2f}")
    print(f"  z_T=inf:  epsilon_T = {result['epsilon_T']:.4f}, "
          f"chi2 = {result['chi2']:.2f}, BIC = {bic_tep:.2f}")
    print(f"  Delta chi2 (inf - 5) = {result['chi2'] - m1_zt5['chi2_mle']:.2f}")
    print(f"  Delta BIC (inf - 5)  = {bic_tep - m1_zt5['bic']:.2f}")
    
    if bic_tep < m1_zt5['bic']:
        print("\n>>> Unscreened TEP has BETTER BIC than screened z_T=5")
    else:
        print("\n>>> Screened z_T=5 is still preferred")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
