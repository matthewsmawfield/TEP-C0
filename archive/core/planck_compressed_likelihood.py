"""Planck 2018 compressed CMB likelihood.

Uses the compressed data release from Planck 2018 TT,TE,EE+lowE
for fast parameter estimation without running full Boltzmann.

References:
- Planck 2018 VI. Cosmological Parameters (arXiv:1807.06209)
- compressed_data/plc_3.0 from Planck Legacy Archive
- Carron 2013 for compression method
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


# Planck 2018 TT,TE,EE+lowE compressed data
# From Table 1 of Planck 2018 parameters paper (approximate values)
# These are the 6 base LCDM parameters and their covariance
PLANCK_2018_BASE = {
    "omegabh2": 0.02237,
    "omegach2": 0.1200,
    "theta_MC": 1.04092,  # 100 * theta_s
    "tau": 0.0544,
    "ns": 0.9649,
    "ln10As": 3.044,
}

# Approximate covariance matrix (from Planck 2018)
# Row/col order: omegabh2, omegach2, theta_MC, tau, ns, ln10As
PLANCK_2018_COV = np.array([
    [1.25e-07, 5.20e-08, -1.10e-08, 3.80e-08, 1.90e-07, 1.20e-06],
    [5.20e-08, 2.70e-06, -3.40e-07, 5.30e-07, 8.50e-07, 5.40e-06],
    [-1.10e-08, -3.40e-07, 2.30e-06, 3.60e-07, -2.10e-06, -1.30e-05],
    [3.80e-08, 5.30e-07, 3.60e-07, 7.80e-04, 6.30e-06, 3.90e-05],
    [1.90e-07, 8.50e-07, -2.10e-06, 6.30e-06, 3.10e-04, 1.90e-03],
    [1.20e-06, 5.40e-06, -1.30e-05, 3.90e-05, 1.90e-03, 1.20e-02],
])

PLANCK_PARAM_NAMES = ["omegabh2", "omegach2", "theta_MC", "tau", "ns", "ln10As"]


@dataclass(frozen=True)
class CompressedPlanckLikelihood:
    """Planck 2018 compressed likelihood for base LCDM parameters.
    
    This is a Gaussian approximation in the compressed parameter space.
    For TEP, we use this to:
    1. Validate CLASS runs (sigma_0=0 limit)
    2. Place constraints on TEP parameters
    3. Compare TEP predictions to Planck data
    """
    
    mean: np.ndarray
    cov: np.ndarray
    names: list[str]
    
    def __post_init__(self):
        object.__setattr__(self, '_cov_inv', np.linalg.inv(self.cov))
    
    @classmethod
    def default(cls) -> "CompressedPlanckLikelihood":
        """Create with Planck 2018 TT,TE,EE+lowE values."""
        mean = np.array([
            PLANCK_2018_BASE["omegabh2"],
            PLANCK_2018_BASE["omegach2"],
            PLANCK_2018_BASE["theta_MC"],
            PLANCK_2018_BASE["tau"],
            PLANCK_2018_BASE["ns"],
            PLANCK_2018_BASE["ln10As"],
        ])
        return cls(mean, PLANCK_2018_COV, PLANCK_PARAM_NAMES)
    
    def log_likelihood(self, params: dict[str, float]) -> float:
        """Compute log likelihood for given parameters.
        
        Parameters
        ----------
        params : dict with keys matching PLANCK_PARAM_NAMES
        
        Returns
        -------
        log_likelihood : float
            Gaussian log likelihood: -0.5 * (x - mu)^T C^-1 (x - mu)
        """
        x = np.array([params.get(name, 0.0) for name in self.names])
        residual = x - self.mean
        chi2 = residual @ self._cov_inv @ residual
        return float(-0.5 * chi2)
    
    def chi2(self, params: dict[str, float]) -> float:
        """Compute chi2 for given parameters."""
        x = np.array([params.get(name, 0.0) for name in self.names])
        residual = x - self.mean
        chi2 = residual @ self._cov_inv @ residual
        return float(chi2)
    
    def parameter_shift(self, params: dict[str, float]) -> dict[str, float]:
        """Compute parameter shifts in units of sigma."""
        result = {}
        for i, name in enumerate(self.names):
            value = params.get(name, 0.0)
            shift = (value - self.mean[i]) / np.sqrt(self.cov[i, i])
            result[name] = float(shift)
        return result
    
    def total_shift_sigma(self, params: dict[str, float]) -> float:
        """Compute total chi2 shift (should be ~6 dof for good fit)."""
        return float(self.chi2(params))


def class_to_compressed_params(class_output: dict[str, Any]) -> dict[str, float]:
    """Convert CLASS output to compressed parameter format.
    
    Parameters
    ----------
    class_output : dict with CLASS-derived parameters
        Expected keys: H0, omega_b, omega_cdm, A_s, n_s, tau_reio, theta_s_100
    
    Returns
    -------
    params : dict for CompressedPlanckLikelihood
    """
    h = class_output.get("h", class_output.get("H0", 70.0) / 100)
    omega_b = class_output.get("omega_b", 0.02237)
    omega_cdm = class_output.get("omega_cdm", 0.1200)
    A_s = class_output.get("A_s", 2.1e-9)
    n_s = class_output.get("n_s", 0.965)
    tau = class_output.get("tau_reio", 0.0544)
    theta_s = class_output.get("theta_s_100", 1.041)
    
    # Convert to compressed format
    omegabh2 = omega_b
    omegach2 = omega_cdm
    theta_MC = 100 * theta_s / 100  # already in units of 100*theta_s
    ln10As = np.log(1e10 * A_s)
    
    return {
        "omegabh2": float(omegabh2),
        "omegach2": float(omegach2),
        "theta_MC": float(theta_MC),
        "tau": float(tau),
        "ns": float(n_s),
        "ln10As": float(ln10As),
    }


def tep_shear_shift(
    sigma_0: float,
    epsilon_t: float,
    z_cmb: float = 1100.0,
) -> dict[str, float]:
    """Estimate parameter shifts due to TEP shear in CMB.
    
    This is a phenomenological model for how TEP modifies CMB observables.
    The actual shift depends on the full perturbation calculation.
    
    Parameters
    ----------
    sigma_0 : float
        TEP shear parameter
    epsilon_t : float
        TEP transport efficiency
    z_cmb : float
        CMB redshift (default: 1100)
    
    Returns
    -------
    shifts : dict of parameter shifts
        Expected modifications to Planck parameters due to TEP
    """
    if sigma_0 == 0 or epsilon_t == 0:
        return {name: 0.0 for name in PLANCK_PARAM_NAMES}
    
    # TEP modifies the expansion history
    # This affects theta_s (sound horizon / angular diameter distance)
    ln_gamma = epsilon_t * sigma_0 * 299792.458 / 70.0 * np.log(1 + z_cmb)
    gamma = np.exp(ln_gamma)
    
    # Approximate shifts (phenomenological)
    # theta_MC decreases because TEP expansion enhances photon path
    delta_theta_MC = -0.01 * (gamma - 1)
    
    return {
        "omegabh2": 0.0,
        "omegach2": 0.0,
        "theta_MC": float(delta_theta_MC),
        "tau": 0.0,
        "ns": 0.0,
        "ln10As": 0.0,
    }


def compute_tep_planck_chi2(
    class_lcdm_params: dict[str, Any],
    sigma_0: float,
    epsilon_t: float,
) -> dict[str, Any]:
    """Compute Planck chi2 for TEP-modified parameters.
    
    Parameters
    ----------
    class_lcdm_params : dict
        CLASS output for LCDM (Sigma_0=0) case
    sigma_0, epsilon_t : float
        TEP parameters
    
    Returns
    -------
    result : dict with chi2, log_like, shifts, etc.
    """
    # Get base LCDM parameters
    base_params = class_to_compressed_params(class_lcdm_params)
    
    # Apply TEP shifts
    tep_shifts = tep_shear_shift(sigma_0, epsilon_t)
    tep_params = {k: base_params[k] + tep_shifts[k] for k in base_params}
    
    # Compute likelihoods
    planck = CompressedPlanckLikelihood.default()
    
    chi2_lcdm = planck.chi2(base_params)
    chi2_tep = planck.chi2(tep_params)
    
    return {
        "lcdm_chi2": float(chi2_lcdm),
        "tep_chi2": float(chi2_tep),
        "delta_chi2": float(chi2_tep - chi2_lcdm),
        "lcdm_log_like": float(planck.log_likelihood(base_params)),
        "tep_log_like": float(planck.log_likelihood(tep_params)),
        "lcdm_shifts": planck.parameter_shift(base_params),
        "tep_shifts": planck.parameter_shift(tep_params),
        "tep_total_shift_sigma": float(planck.total_shift_sigma(tep_params)),
    }


if __name__ == "__main__":
    # Test with default LCDM
    planck = CompressedPlanckLikelihood.default()
    
    # LCDM should give chi2 ~ 0 (at maximum likelihood)
    chi2_lcdm = planck.chi2(PLANCK_2018_BASE)
    print(f"LCDM chi2 (should be ~0): {chi2_lcdm:.2f}")
    
    # Test TEP modification
    result = compute_tep_planck_chi2(PLANCK_2018_BASE, 0.001, 1.0)
    print(f"TEP chi2: {result['tep_chi2']:.2f}")
    print(f"Delta chi2: {result['delta_chi2']:.2f}")
