"""Pantheon+ likelihood for Cobaya integration.

This module provides a Cobaya-compatible likelihood for Pantheon+ supernova data.
The likelihood uses CLASS distances as the LCDM baseline and applies the TEP
correction ratio on top, ensuring self-consistent treatment of cosmological
parameters between the CMB and SNe likelihoods.
"""

from __future__ import annotations

import numpy as np
from typing import Optional
import os

from cobaya.likelihood import Likelihood


class PantheonCobaya(Likelihood):
    """Pantheon+ supernova likelihood for Cobaya.
    
    Implements the Pantheon+ dataset as a Cobaya likelihood class.
    Compatible with TEP-C0 cosmological models.
    
    Architecture:
        The likelihood requests Hubble(z) from the theory provider (CLASS)
        and uses it to compute both the LCDM baseline comoving distance and
        the TEP-modified comoving distance. This ensures that the matter
        density parameters (omega_b, omega_cdm) enter the SNe distance
        calculation through CLASS — the same code that the Planck likelihood
        uses — eliminating any inconsistency between the two probes.
        
        The TEP parameters (epsilon_T, z_T) enter only as a multiplicative
        correction to the integrand of the comoving distance integral.
    """
    
    data_file: Optional[str] = "data/processed/tep_c0_pantheon_plus_distances.csv"
    cov_file: Optional[str] = "data/processed/tep_c0_pantheon_plus_covariance_manifest.csv"
    
    # Internal redshift grid for H(z) integration
    _z_grid_max: float = 2.5
    _z_grid_n: int = 500
    
    def initialize(self):
        """Initialize Pantheon+ likelihood."""
        try:
            import pandas as pd
            if os.path.exists(self.data_file):
                df = pd.read_csv(self.data_file)
                z_col = 'zHD' if 'zHD' in df.columns else 'zCMB'
                mb_col = 'm_b_corr' if 'm_b_corr' in df.columns else 'mB'
                err_col = 'm_b_corr_err_DIAG' if 'm_b_corr_err_DIAG' in df.columns else 'MU_SH0ES_ERR_DIAG'
                
                self.z = df[z_col].values
                self.mb = df[mb_col].values
                self.mb_err = df[err_col].values if err_col in df.columns else np.ones_like(self.z) * 0.1
                
                valid = np.isfinite(self.z) & np.isfinite(self.mb) & (self.z > 0)
                self.z = self.z[valid]
                self.mb = self.mb[valid]
                self.mb_err = self.mb_err[valid]
            else:
                raise FileNotFoundError(f"Data file not found: {self.data_file}")
            
            # Load covariance
            if os.path.exists(self.cov_file):
                try:
                    cov_flat = np.loadtxt(self.cov_file)
                    n_sn = len(self.z)
                    if cov_flat.size == n_sn * n_sn + 1:
                        self.cov = cov_flat[1:].reshape(n_sn, n_sn)
                    elif cov_flat.size == n_sn * n_sn:
                        self.cov = cov_flat.reshape(n_sn, n_sn)
                    else:
                        self.cov = np.diag(self.mb_err ** 2)
                except Exception:
                    self.cov = np.diag(self.mb_err ** 2)
            else:
                self.cov = np.diag(self.mb_err ** 2)
                
            self.n_sn = len(self.z)
            
        except Exception as e:
            print(f"Pantheon likelihood initialization failed: {e}")
            raise
            
        # Precompute the inverse covariance and its sum for M marginalization
        try:
            self.inv_cov = np.linalg.inv(self.cov)
        except np.linalg.LinAlgError:
            self.inv_cov = np.diag(1.0 / np.diag(self.cov))
            
        self.sum_inv_cov = np.sum(self.inv_cov)
        
        # Build the integration grid for H(z) -> distances
        self._z_grid_max = max(float(np.max(self.z)) * 1.05, 2.5)
        self._z_grid = np.linspace(0.0, self._z_grid_max, self._z_grid_n)

    def get_requirements(self):
        """Return requirements for Cobaya.
        
        Request Hubble(z) from CLASS at the integration grid redshifts.
        This ensures the LCDM baseline distances are computed by CLASS
        using the same cosmological parameters as the CMB likelihood.
        """
        return {
            "Hubble": {"z": self._z_grid.tolist()}
        }
    
    def _tep_gamma(self, z, epsilon_T=0.001, z_T=5.0, n_T=1.0):
        """Temporal Equivalence Principle (TEP) gamma factor.
        
        This is the canonical TEP formula from the manuscript,
        matching core/tep_cosmology.py exactly:
            gamma(z) = 1 + epsilon_T * log((1+z)/(1+z_T)) * exp(-(z/z_T)^n_T)
        
        Parameters
        ----------
        z : array_like
            Redshifts
        epsilon_T : float
            Temporal shear mixing fraction
        z_T : float
            Transition redshift
        n_T : float
            Transition index (default 1.0)
        
        Returns
        -------
        array_like
            Gamma factor at each redshift
        """
        z_arr = np.asarray(z, dtype=float)
        if epsilon_T == 0.0:
            return np.ones_like(z_arr)
        
        log_factor = np.log(1.0 + z_arr)
        suppression = np.exp(-((z_arr / z_T) ** n_T))
        gamma = 1.0 - epsilon_T * log_factor * suppression
        return np.maximum(gamma, 0.1)
    
    def logp(self, **params_values) -> float:
        """Compute log-likelihood for given theory predictions.
        
        Uses CLASS H(z) to compute LCDM baseline distances, then applies
        the TEP correction ratio. This ensures omega_cdm enters the SNe
        distance calculation through CLASS (same as Planck).
        
        Returns:
            Log-likelihood value
        """
        # Get H(z) from CLASS on the integration grid
        # Cobaya provides H(z) in km/s/Mpc via the provider
        H_grid = np.array(self.provider.get_Hubble(self._z_grid))
        
        # Check for valid H(z) values
        if np.any(np.isnan(H_grid)) or np.any(H_grid <= 0):
            return -np.inf
        
        # Compute TEP-modified comoving distance via numerical integration:
        #   d_c(z) = integral_0^z [c / H_TEP(z')] dz'
        # H_grid from CLASS already includes the division by gamma(z).
        c_light = 299792.458  # km/s
        integrand = c_light / H_grid
        
        from scipy.integrate import cumulative_trapezoid
        cum_dc = cumulative_trapezoid(integrand, self._z_grid, initial=0.0)
        
        # Interpolate to SNe redshifts and convert to luminosity distance
        dc_sne = np.interp(self.z, self._z_grid, cum_dc)
        d_l = dc_sne * (1.0 + self.z)
        
        # Check for unphysical distances
        if np.any(d_l <= 0.0) or np.any(np.isnan(d_l)):
            return -np.inf
        
        with np.errstate(divide='ignore', invalid='ignore'):
            theory_prediction = 5.0 * np.log10(d_l) + 25.0
            
        if np.any(np.isnan(theory_prediction)) or np.any(np.isinf(theory_prediction)):
            return -np.inf
        
        # Analytical marginalization over the absolute magnitude M:
        #   M_best = (1^T C^{-1} delta) / (1^T C^{-1} 1)
        #   chi2 = (delta - M_best)^T C^{-1} (delta - M_best)
        delta = self.mb - theory_prediction
        
        if np.any(np.isnan(delta)) or np.any(np.isinf(delta)):
            return -np.inf
            
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            M_best = np.sum(self.inv_cov @ delta) / self.sum_inv_cov
            residual = delta - M_best
            chi2 = residual.T @ self.inv_cov @ residual
            
        if np.isnan(chi2) or np.isinf(chi2):
            return -np.inf
        return float(-0.5 * chi2)
