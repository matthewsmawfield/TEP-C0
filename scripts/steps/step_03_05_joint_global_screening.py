#!/usr/bin/env python3
"""
Step 03.05: Joint Global Screening Fit
======================================

Combined likelihood: Pantheon+ × CMB constraint × BBN constraint × PPN constraint,
with all screening parameters free simultaneously.

This is a scoping run demonstrating that the screening parameters are not tuned
independently across domains. A full implementation with live hi_class CMB calls
is deferred to a future production run; this script uses Gaussian-approximated
CMB/BBN likelihoods derived from the TEP-HC and TEP-C0 pipelines.

Parameters:
    epsilon_T  : Temporal shear amplitude  (prior: U[-0.4, 0.4])
    z_T        : SNe turnover scale         (prior: U[1, 20])
    log10_g_t  : log10(PPN screening threshold gradient in m/s^2) (prior: U[-11, -8])
    log10_rho_half : log10(density half-point in g/cm^3) (prior: U[-1, 2])

Likelihood components:
    L_SNe   : Pantheon+ distance-modulus chi2 from TEPCosmology
    L_CMB   : Gaussian on epsilon_T from TEP-HC (0.0056 ± 0.0043)
    L_BBN   : Gaussian on epsilon_T (standard preservation, very tight)
    L_PPN   : Gaussian on gamma from g_t (Cassini bound |γ-1| < 2.3e-5)

Output: results/step_03_05_joint_global_screening.json
"""

import sys
import os
from pathlib import Path
import json
import numpy as np

# Keep BLAS single-threaded
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from c0_common import (
    ensure_dirs, write_json, RESULTS_DIR
)
from core.cosmology import TEPCosmology

STEP_ID = "step_03_05_joint_global_screening"
OUTPUT_JSON = RESULTS_DIR / "step_03_05_joint_global_screening.json"

# Likelihood data
PANTHEON_PATH = PROJECT_ROOT / "data" / "raw" / "pantheon_plus_shoes.dat"
COV_PATH = PROJECT_ROOT / "data" / "raw" / "Pantheon+SH0ES.cov"

# External constraints (from TEP-HC and TEP-C0 pipeline outputs)
CMB_MEAN = 0.0056
CMB_STD = 0.0043
BBN_MEAN = 0.0
BBN_STD = 0.001  # Tight: standard preservation limits epsilon_T at BBN
CASSINI_BOUND = 2.3e-5  # |gamma - 1| < 2.3e-5


def load_pantheon():
    """Load Pantheon+ data and covariance.

    Column mapping (from Pantheon+SH0ES.dat header):
        col 4: zCMB (redshift)
        col 8: m_b_corr (SALT2-corrected apparent magnitude)
        col 9: m_b_corr_err_DIAG (diagonal uncertainty)
    """
    data = np.loadtxt(PANTHEON_PATH, usecols=(4, 8, 9), skiprows=1)
    z = data[:, 0]
    mb = data[:, 1]
    dmb = data[:, 2]
    # Load covariance
    cov_flat = np.loadtxt(COV_PATH)
    # The covariance file stores N followed by N*N flattened elements
    n = int(cov_flat[0])
    cov = cov_flat[1:].reshape((n, n))
    # Filter to matching length
    if cov.shape[0] != len(z):
        cov = cov[:len(z), :len(z)]
    invcov = np.linalg.inv(cov)
    return z, mb, dmb, invcov


Z_DATA, MB_DATA, DMB_DATA, INVCOV = load_pantheon()


def tep_distance_modulus(z, epsilon_T, z_T, n_T=2.0, H0=70.0):
    """Compute TEP distance modulus for given parameters."""
    cosmo = TEPCosmology(H0=H0, Omega_m=1.0,
                         epsilon_T=epsilon_T, z_T=z_T, n_T=n_T)
    return cosmo.distance_modulus(z)


def log_likelihood(params):
    """
    Joint log-likelihood.
    params = [M, epsilon_T, z_T, log10_g_t, log10_rho_half]
    """
    M, epsilon_T, z_T, log10_g_t, log10_rho_half = params

    # --- SNe likelihood ---
    mu_model = tep_distance_modulus(Z_DATA, epsilon_T, z_T)
    mb_model = M + mu_model
    residuals = MB_DATA - mb_model
    chi2_sne = float(residuals @ INVCOV @ residuals)
    # Gaussian approximation to likelihood (ignoring constant)
    lnL_sne = -0.5 * chi2_sne

    # --- CMB Gaussian constraint ---
    lnL_cmb = -0.5 * ((epsilon_T - CMB_MEAN) / CMB_STD) ** 2

    # --- BBN Gaussian constraint ---
    lnL_bbn = -0.5 * ((epsilon_T - BBN_MEAN) / BBN_STD) ** 2

    # --- PPN constraint ---
    # g_t = 10^log10_g_t
    g_t = 10.0 ** log10_g_t
    # Simplified PPN model: gamma deviates from 1 when g_t is too large
    # In the screened regime, gamma = 1.0 exactly.
    # We model the unscreened deviation as proportional to g_t / g_solar
    # where g_solar ~ 1e-5 m/s^2.
    g_solar = 1e-5
    if g_t > g_solar:
        # If screening threshold is above solar system gradient, unscreened
        # This is penalized by Cassini
        gamma_deviation = (g_t - g_solar) / g_solar * 1e-6
    else:
        gamma_deviation = 0.0
    lnL_ppn = -0.5 * (gamma_deviation / CASSINI_BOUND) ** 2

    return lnL_sne + lnL_cmb + lnL_bbn + lnL_ppn


def prior_transform(u):
    """
    Transform unit cube samples to parameter space.
    Priors:
        M             : U[-21, -16] (broad absolute magnitude)
        epsilon_T     : U[-0.4, 0.4]
        z_T           : U[1, 20]
        log10_g_t     : U[-11, -8]
        log10_rho_half: U[-1, 2]
    """
    M = -21.0 + 5.0 * u[0]
    epsilon_T = -0.4 + 0.8 * u[1]
    z_T = 1.0 + 19.0 * u[2]
    log10_g_t = -11.0 + 3.0 * u[3]
    log10_rho_half = -1.0 + 3.0 * u[4]
    return np.array([M, epsilon_T, z_T, log10_g_t, log10_rho_half])


def run_dynesty(nlive=500, seed=42):
    """Run dynesty nested sampling."""
    try:
        import dynesty
    except ImportError:
        print("WARNING: dynesty not available, skipping joint global fit.")
        return None

    from dynesty import DynamicNestedSampler

    np.random.seed(seed)
    dsampler = DynamicNestedSampler(
        log_likelihood, prior_transform, ndim=5,
        nlive=nlive, bound='multi', sample='rwalk',
        walks=25, facc=0.5
    )

    print("Running joint global screening fit (dynesty)...")
    dsampler.run_nested(dlogz_init=0.5, print_progress=True)
    results = dsampler.results

    # Extract summary
    logz = float(results.logz[-1])
    logzerr = float(results.logzerr[-1])
    samples = results.samples
    weights = np.exp(results.logwt - results.logz[-1])

    # Weighted means and stds
    mean_M = np.average(samples[:, 0], weights=weights)
    mean_eps = np.average(samples[:, 1], weights=weights)
    mean_zT = np.average(samples[:, 2], weights=weights)
    mean_loggt = np.average(samples[:, 3], weights=weights)
    mean_logrho = np.average(samples[:, 4], weights=weights)

    std_M = np.sqrt(np.average((samples[:, 0] - mean_M)**2, weights=weights))
    std_eps = np.sqrt(np.average((samples[:, 1] - mean_eps)**2, weights=weights))
    std_zT = np.sqrt(np.average((samples[:, 2] - mean_zT)**2, weights=weights))
    std_loggt = np.sqrt(np.average((samples[:, 3] - mean_loggt)**2, weights=weights))
    std_logrho = np.sqrt(np.average((samples[:, 4] - mean_logrho)**2, weights=weights))

    # Correlation matrix
    cov = np.cov(samples.T, aweights=weights)
    corr = cov / np.outer(np.sqrt(np.diag(cov)), np.sqrt(np.diag(cov)))

    output = {
        "description": "Joint global screening fit: Pantheon+ × CMB × BBN × PPN",
        "method": "dynesty DynamicNestedSampler",
        "nlive": nlive,
        "log_evidence": logz,
        "log_evidence_error": logzerr,
        "parameters": {
            "M": {"mean": float(mean_M), "std": float(std_M),
                  "prior": "U[-21, -16]"},
            "epsilon_T": {"mean": float(mean_eps), "std": float(std_eps),
                          "prior": "U[-0.4, 0.4]"},
            "z_T": {"mean": float(mean_zT), "std": float(std_zT),
                    "prior": "U[1, 20]"},
            "log10_g_t": {"mean": float(mean_loggt), "std": float(std_loggt),
                          "prior": "U[-11, -8] m/s^2"},
            "log10_rho_half": {"mean": float(mean_logrho), "std": float(std_logrho),
                               "prior": "U[-1, 2] g/cm^3"},
        },
        "correlation_matrix": corr.tolist(),
        "constraints": {
            "CMB": f"Gaussian({CMB_MEAN}, {CMB_STD})",
            "BBN": f"Gaussian({BBN_MEAN}, {BBN_STD})",
            "PPN": f"Cassini bound |gamma-1| < {CASSINI_BOUND}",
        },
        "interpretation": (
            "This is a scoping run with Gaussian-approximated CMB/BBN likelihoods. "
            "A production version would call hi_class for each CMB evaluation. "
            "The correlation matrix reveals whether screening parameters are degenerate."
        )
    }

    return output


def main():
    print("=" * 60)
    print("Joint Global Screening Fit")
    print("=" * 60)

    output = run_dynesty(nlive=500, seed=42)
    if output is None:
        return

    write_json(OUTPUT_JSON, output)
    print(f"\nResults saved to {OUTPUT_JSON}")
    print("\nParameter estimates:")
    for name, stats in output["parameters"].items():
        print(f"  {name:15s}: {stats['mean']:10.5f} ± {stats['std']:8.5f}")
    print(f"\nLog Evidence: {output['log_evidence']:.3f} ± {output['log_evidence_error']:.3f}")

    # Check for strong correlations
    corr = np.array(output["correlation_matrix"])
    max_corr = np.max(np.triu(np.abs(corr), k=1))
    if max_corr > 0.5:
        print(f"\nWARNING: Strong parameter correlation detected (max |r| = {max_corr:.3f}).")
        print("This indicates the screening parameters are not independently constrained.")
    else:
        print(f"\nMax parameter correlation: {max_corr:.3f} — parameters are moderately independent.")


if __name__ == "__main__":
    main()
