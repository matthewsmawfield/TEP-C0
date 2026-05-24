"""Pantheon+ SNe Ia likelihood for Cobaya.

Integrates Pantheon+ data with Cobaya's likelihood interface for joint
parameter estimation with TEP-CLASS.

This module reuses the data loading infrastructure from step_022 to ensure
consistency with the existing pipeline.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

# Add parent to path for step_022 imports
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    from step_022_three_model_comparison import PantheonData
    HAS_STEP022 = True
except ImportError:
    HAS_STEP022 = False

# Try to import Cobaya base class
try:
    from cobaya.likelihood import Likelihood
    HAS_COBAYA = True
except ImportError:
    HAS_COBAYA = False
    # Create dummy base class for type checking
    class Likelihood:  # type: ignore
        """Dummy Likelihood base class when Cobaya not available."""
        pass


class PantheonCobaya(Likelihood):
    """Pantheon+ SNe Ia likelihood for Cobaya.
    
    Uses the PantheonData class from step_022 for consistent data handling.
    Computes chi2 from distance modulus residuals with full covariance.
    
    Attributes:
        z_sn: Redshifts of SNe
        mu_sn: Distance moduli (SH0ES-calibrated)
        inv_cov: Inverse covariance matrix
    """
    
    # Cobaya configuration
    likelihood_name = "pantheon_plus"
    
    def initialize(self):
        """Initialize Pantheon+ data.
        
        Called by Cobaya during model setup.
        """
        if not HAS_STEP022:
            raise RuntimeError("step_022_three_model_comparison not available")
        
        # Load Pantheon data using step_022 infrastructure
        self.data = PantheonData()
        
        try:
            self.data.load()
        except Exception as e:
            raise RuntimeError(f"Failed to load Pantheon+ data: {e}")
        
        # Extract arrays
        self.z_sn = self.data.z
        self.mu_sn = self.data.mb  # These are actually m_b_corr, not mu
        
        # Convert apparent magnitudes to distance moduli using absolute magnitude
        # M_B from SH0ES calibration ~ -19.25
        self.MB_fixed = -19.25
        self.mu_sn = self.mu_sn - self.MB_fixed
        
        # Use covariance from PantheonData
        self.cov = self.data.cov
        self.inv_cov = self.data.cov_cholesky_inv if hasattr(self.data, 'cov_cholesky_inv') else np.linalg.inv(self.cov)
        
        self.n_sn = len(self.z_sn)
        
        self.log.info(f"Pantheon+ loaded: {self.n_sn} SNe Ia")
        self.log.info(f"Redshift range: [{self.z_sn.min():.4f}, {self.z_sn.max():.4f}]")
    
    def get_requirements(self):
        """What we need from the theory code (CLASS)."""
        return {
            "angular_diameter_distance": {"z": list(self.z_sn)}
        }
    
    def _compute_distance_modulus(self, z: np.ndarray, Da: np.ndarray) -> np.ndarray:
        """Compute distance modulus from angular diameter distance.
        
        For TEP: mu = 5*log10(D_L) + 25 where D_L = D_A * (1+z)^2
        
        Args:
            z: Redshifts
            Da: Angular diameter distances in Mpc
            
        Returns:
            Distance moduli
        """
        # D_L = D_A * (1+z)^2
        D_L = Da * (1 + z) ** 2
        
        # Distance modulus
        mu = 5.0 * np.log10(D_L) + 25.0
        return mu
    
    def logp(self, **params_values):
        """Compute log-likelihood.
        
        Called by Cobaya during sampling.
        
        Args:
            **params_values: Parameters from theory
            
        Returns:
            Log-likelihood value
        """
        # Get angular diameter distances from theory (CLASS/TEP-CLASS)
        Da_list = []
        for z in self.z_sn:
            Da = self.provider.get_angular_diameter_distance(z)
            # Ensure scalar
            if hasattr(Da, '__len__'):
                Da = float(Da[0]) if len(Da) > 0 else float(Da)
            else:
                Da = float(Da)
            Da_list.append(Da)
        
        Da = np.array(Da_list)
        
        # Compute theoretical distance moduli
        mu_theory = self._compute_distance_modulus(self.z_sn, Da)
        
        # Residuals: data - theory
        residuals = self.mu_sn - mu_theory
        
        # Chi2 with covariance
        # Use precomputed Cholesky if available for efficiency
        if hasattr(self.data, 'cov_cholesky') and self.data.cov_cholesky is not None:
            from scipy.linalg import solve_triangular
            y = solve_triangular(self.data.cov_cholesky, residuals, lower=True)
            chi2 = np.dot(y, y)
        else:
            chi2 = residuals @ self.inv_cov @ residuals
        
        # Log-likelihood (neglecting normalization)
        return -0.5 * chi2
    
    def get_can_provide_params(self):
        """Parameters this likelihood can provide."""
        return []


# Backward compatibility: standalone test function
def test_likelihood():
    """Test the likelihood without Cobaya."""
    if not HAS_STEP022:
        print("step_022 not available")
        return
    
    data = PantheonData()
    data.load()
    
    print(f"Loaded {len(data.z)} SNe")
    print(f"Redshift range: [{data.z.min():.4f}, {data.z.max():.4f}]")
    print(f"Covariance shape: {data.cov.shape}")
    
    # Test chi2 computation with null model
    MB = -19.25
    mu_data = data.mb - MB
    residuals = mu_data - mu_data  # Null residuals
    
    from scipy.linalg import cholesky, solve_triangular
    cov_cholesky = cholesky(data.cov, lower=True)
    y = solve_triangular(cov_cholesky, residuals, lower=True)
    chi2 = np.dot(y, y)
    
    print(f"Test chi2 (null): {chi2:.2f}")


if __name__ == "__main__":
    test_likelihood()
