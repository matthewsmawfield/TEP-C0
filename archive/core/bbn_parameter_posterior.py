"""BBN parameter posterior estimation using MCMC.

Fits baryon-to-photon ratio (eta), neutron lifetime (tau_n), and
effective neutrino species (N_eff) against the light-element abundance registry.

References:
- Pitrou et al. 2021, Phys. Rep. (for BBN review)
- Particle Data Group 2022, neutron lifetime
- Planck 2018, eta constraint
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .bbn_network_resolver import BBNInputs, run_network


@dataclass(frozen=True)
class BBNPosteriorConfig:
    """Configuration for BBN parameter posterior estimation."""
    eta_min: float = 5.0e-10
    eta_max: float = 7.0e-10
    tau_n_min: float = 850.0
    tau_n_max: float = 900.0
    n_nu_min: float = 2.0
    n_nu_max: float = 4.0
    sigma_0: float = 0.0
    epsilon_t: float = 0.0
    h0_km_s_mpc: float = 70.0
    n_walkers: int = 32
    n_steps: int = 2000
    burn_in: int = 500


def load_bbn_registry(registry_path: Path) -> list[dict[str, Any]]:
    """Load BBN abundance registry CSV."""
    with registry_path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def compute_abundance_chi2(
    abundances: dict[str, float],
    registry_rows: list[dict[str, Any]],
) -> tuple[float, dict[str, float]]:
    """Compute chi2 against observational registry."""
    mapping = {
        "helium_4_mass_fraction_Yp": "Y_p",
        "deuterium_to_hydrogen_DH": "D_H",
        "helium_3_to_hydrogen_He3H": "He3_H",
        "lithium_7_to_hydrogen_Li7H": "Li7_H",
    }
    chi2_terms: dict[str, float] = {}
    chi2 = 0.0
    
    for row in registry_rows:
        label = mapping.get(str(row.get("quantity")))
        if not label or label not in abundances:
            continue
        
        obs_value = float(row["observed_value"])
        obs_sigma = float(row["observed_value"])
        pred_sigma = float(row.get("prediction_sigma", 1e-10))
        
        # Combined uncertainty (observational + prediction systematic)
        sigma_tot = np.sqrt(obs_sigma**2 + pred_sigma**2)
        if sigma_tot < 1e-20:
            continue
        
        residual = abundances[label] - obs_value
        term = (residual / sigma_tot) ** 2
        chi2_terms[label] = float(term)
        chi2 += term
    
    return float(chi2), chi2_terms


def bbn_log_posterior(
    theta: np.ndarray,
    registry_rows: list[dict[str, Any]],
    config: BBNPosteriorConfig,
) -> float:
    """Log posterior for BBN parameters.
    
    Parameters
    ----------
    theta : array [eta_log10, tau_n, n_nu]
        eta_log10 = log10(eta * 1e10)
        tau_n = neutron lifetime in seconds
        n_nu = effective number of neutrino species
    """
    eta_log10, tau_n, n_nu = theta
    
    # Priors (uniform within bounds)
    eta = 10 ** (eta_log10 - 10)
    if not (config.eta_min <= eta <= config.eta_max):
        return -np.inf
    if not (config.tau_n_min <= tau_n <= config.tau_n_max):
        return -np.inf
    if not (config.n_nu_min <= n_nu <= config.n_nu_max):
        return -np.inf
    
    # Run BBN network
    inputs = BBNInputs(
        eta=eta,
        h0_km_s_mpc=config.h0_km_s_mpc,
        sigma_0=config.sigma_0,
        epsilon_t=config.epsilon_t,
        n_nu=n_nu,
        tau_n_s=tau_n,
    )
    
    try:
        result = run_network(inputs)
        abundances = result["abundances"]
    except Exception:
        return -np.inf
    
    # Compute likelihood
    chi2, _ = compute_abundance_chi2(abundances, registry_rows)
    # Gaussian likelihood
    log_like = -0.5 * chi2
    
    return log_like


def run_bbn_mcmc(
    registry_rows: list[dict[str, Any]],
    config: BBNPosteriorConfig,
    seed: int = 42,
) -> dict[str, Any]:
    """Run MCMC for BBN parameter posterior.
    
    Returns
    -------
    Dictionary with posterior samples and summary statistics.
    """
    try:
        import emcee
    except ImportError:
        # Fallback to grid sampling
        return _grid_sample_bbn(registry_rows, config)
    
    np.random.seed(seed)
    
    # Initialize walkers with better spread to avoid condition number issues
    n_walkers, n_dim = config.n_walkers, 3
    
    # Bounds for initialization
    eta_log10_min = np.log10(config.eta_min * 1e10)
    eta_log10_max = np.log10(config.eta_max * 1e10)
    
    # Initialize walkers uniformly across parameter space
    # This ensures linear independence and avoids condition number warning
    pos = np.zeros((n_walkers, n_dim))
    pos[:, 0] = np.random.uniform(eta_log10_min, eta_log10_max, n_walkers)
    pos[:, 1] = np.random.uniform(config.tau_n_min, config.tau_n_max, n_walkers)
    pos[:, 2] = np.random.uniform(config.n_nu_min, config.n_nu_max, n_walkers)
    
    def log_prob(theta):
        return bbn_log_posterior(theta, registry_rows, config)
    
    sampler = emcee.EnsembleSampler(n_walkers, n_dim, log_prob)
    sampler.run_mcmc(pos, config.n_steps, progress=False)
    
    # Extract samples
    samples = sampler.get_chain(discard=config.burn_in, flat=True)
    
    # Compute Gelman-Rubin
    try:
        chain = sampler.get_chain(discard=config.burn_in, flat=False)
        r_hat = _compute_r_hat(chain)
    except Exception:
        r_hat = np.nan
    
    # Summary statistics
    eta_samples = 10 ** (samples[:, 0] - 10)
    tau_n_samples = samples[:, 1]
    n_nu_samples = samples[:, 2]
    
    return {
        "method": "emcee_mcmc",
        "n_walkers": n_walkers,
        "n_steps": config.n_steps,
        "burn_in": config.burn_in,
        "r_hat": float(r_hat),
        "converged": r_hat < 1.1 if not np.isnan(r_hat) else False,
        "eta": {
            "median": float(np.median(eta_samples)),
            "mean": float(np.mean(eta_samples)),
            "std": float(np.std(eta_samples)),
            "16th_percentile": float(np.percentile(eta_samples, 16)),
            "84th_percentile": float(np.percentile(eta_samples, 84)),
            "samples": eta_samples.tolist(),
        },
        "tau_n": {
            "median": float(np.median(tau_n_samples)),
            "mean": float(np.mean(tau_n_samples)),
            "std": float(np.std(tau_n_samples)),
            "16th_percentile": float(np.percentile(tau_n_samples, 16)),
            "84th_percentile": float(np.percentile(tau_n_samples, 84)),
            "samples": tau_n_samples.tolist(),
        },
        "n_nu": {
            "median": float(np.median(n_nu_samples)),
            "mean": float(np.mean(n_nu_samples)),
            "std": float(np.std(n_nu_samples)),
            "16th_percentile": float(np.percentile(n_nu_samples, 16)),
            "84th_percentile": float(np.percentile(n_nu_samples, 84)),
            "samples": n_nu_samples.tolist(),
        },
    }


def _grid_sample_bbn(
    registry_rows: list[dict[str, Any]],
    config: BBNPosteriorConfig,
) -> dict[str, Any]:
    """Fallback grid sampling when emcee not available."""
    eta_grid = np.linspace(config.eta_min, config.eta_max, 8)
    tau_n_grid = np.linspace(config.tau_n_min, config.tau_n_max, 8)
    n_nu_grid = np.linspace(config.n_nu_min, config.n_nu_max, 8)
    
    log_posterior = np.zeros((len(eta_grid), len(tau_n_grid), len(n_nu_grid)))
    
    for i, eta in enumerate(eta_grid):
        for j, tau_n in enumerate(tau_n_grid):
            for k, n_nu in enumerate(n_nu_grid):
                theta = np.array([np.log10(eta * 1e10), tau_n, n_nu])
                log_posterior[i, j, k] = bbn_log_posterior(theta, registry_rows, config)
    
    # Find maximum
    max_idx = np.unravel_index(np.argmax(log_posterior), log_posterior.shape)
    eta_best = eta_grid[max_idx[0]]
    tau_n_best = tau_n_grid[max_idx[1]]
    n_nu_best = n_nu_grid[max_idx[2]]
    
    # Compute approximate uncertainties from curvature
    # (simplified - full MCMC would give better errors)
    
    return {
        "method": "grid_sample",
        "eta": {
            "median": float(eta_best),
            "mean": float(eta_best),
            "std": float(eta_grid[1] - eta_grid[0]),
            "16th_percentile": float(eta_best - (eta_grid[1] - eta_grid[0])),
            "84th_percentile": float(eta_best + (eta_grid[1] - eta_grid[0])),
        },
        "tau_n": {
            "median": float(tau_n_best),
            "mean": float(tau_n_best),
            "std": float(tau_n_grid[1] - tau_n_grid[0]),
            "16th_percentile": float(tau_n_best - (tau_n_grid[1] - tau_n_grid[0])),
            "84th_percentile": float(tau_n_best + (tau_n_grid[1] - tau_n_grid[0])),
        },
        "n_nu": {
            "median": float(n_nu_best),
            "mean": float(n_nu_best),
            "std": float(n_nu_grid[1] - n_nu_grid[0]),
            "16th_percentile": float(n_nu_best - (n_nu_grid[1] - n_nu_grid[0])),
            "84th_percentile": float(n_nu_best + (n_nu_grid[1] - n_nu_grid[0])),
        },
    }


def _compute_r_hat(chains: np.ndarray) -> float:
    """Compute Gelman-Rubin R-hat statistic."""
    n_chains, n_steps, n_dim = chains.shape
    
    # Chain means
    chain_means = np.mean(chains, axis=1)
    
    # Overall mean
    overall_mean = np.mean(chain_means, axis=0)
    
    # Between-chain variance
    B = n_steps * np.var(chain_means, axis=1, ddof=1)
    
    # Within-chain variance
    W = np.mean(np.var(chains, axis=1, ddof=1), axis=0)
    
    # R-hat
    V_hat = ((n_steps - 1) / n_steps) * W + B / n_steps
    r_hat = np.sqrt(V_hat / W)
    
    return float(np.max(r_hat))


if __name__ == "__main__":
    registry = load_bbn_registry(
        Path(__file__).resolve().parents[3] / "data" / "processed" / "tep_c0_bbn_abundance_registry.csv"
    )
    
    config = BBNPosteriorConfig(
        sigma_0=0.0,
        epsilon_t=0.0,
        n_walkers=16,
        n_steps=500,
        burn_in=100,
    )
    
    result = run_bbn_mcmc(registry, config)
    print(f"eta = {result['eta']['median']:.2e} +/- {result['eta']['std']:.2e}")
    print(f"tau_n = {result['tau_n']['median']:.1f} +/- {result['tau_n']['std']:.1f} s")
    print(f"N_nu = {result['n_nu']['median']:.2f} +/- {result['n_nu']['std']:.2f}")
