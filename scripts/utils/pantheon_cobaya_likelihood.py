"""Pantheon+ likelihood for Cobaya integration.

This module provides a Cobaya-compatible likelihood for Pantheon+ supernova data.
"""

from __future__ import annotations

import numpy as np
from typing import Dict, Any, Optional
import os


from cobaya.likelihood import Likelihood

class PantheonCobaya(Likelihood):
    """Pantheon+ supernova likelihood for Cobaya.
    
    Implements the Pantheon+ dataset as a Cobaya likelihood class.
    Compatible with TEP-C0 cosmological models.
    """
    
    data_file: Optional[str] = "data/processed/tep_c0_pantheon_plus_distances.csv"
    cov_file: Optional[str] = "data/processed/tep_c0_pantheon_plus_covariance_manifest.csv"
    raw_cov_file: Optional[str] = "data/raw/Pantheon+SH0ES.cov"
    
    def initialize(self):
        """Initialize Pantheon+ likelihood."""
        # Load data
        try:
            import pandas as pd
            if os.path.exists(self.data_file):
                df = pd.read_csv(self.data_file)
                # Determine correct columns based on Pantheon+ SH0ES format
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
            
            # Load covariance — try processed manifest, then raw .cov, then diagonal fallback
            cov_loaded = False
            cov_source = None
            for cov_path in [self.cov_file, self.raw_cov_file]:
                if cov_path and os.path.exists(cov_path):
                    try:
                        cov_flat = np.loadtxt(cov_path)
                        n_sn = len(self.z)
                        if cov_flat.size == n_sn * n_sn + 1:
                            self.cov = cov_flat[1:].reshape(n_sn, n_sn)
                            cov_loaded = True
                            cov_source = cov_path
                        elif cov_flat.size == n_sn * n_sn:
                            self.cov = cov_flat.reshape(n_sn, n_sn)
                            cov_loaded = True
                            cov_source = cov_path
                        
                        # Symmetrize if slightly asymmetric (floating-point roundoff in file)
                        if cov_loaded and not np.allclose(self.cov, self.cov.T):
                            self.cov = (self.cov + self.cov.T) / 2.0
                    except Exception as e:
                        print(f"[PantheonCobaya] WARNING: Failed to load covariance from {cov_path}: {e}")
                    if cov_loaded:
                        break
            if not cov_loaded:
                print(f"[PantheonCobaya] WARNING: Full covariance not found. Using diagonal approximation ({n_sn} SNe).")
                self.cov = np.diag(self.mb_err ** 2)
                cov_source = "diagonal_fallback"
            else:
                print(f"[PantheonCobaya] Loaded full covariance from {cov_source} ({n_sn} x {n_sn})")
                
            self.n_sn = len(self.z)
            
        except Exception as e:
            print(f"Pantheon likelihood initialization failed: {e}")
            raise
            
        # Precompute the inverse covariance and its sum to save time in logp
        try:
            self.inv_cov = np.linalg.inv(self.cov)
        except np.linalg.LinAlgError:
            self.inv_cov = np.diag(1.0 / np.diag(self.cov))
            
        self.sum_inv_cov = np.sum(self.inv_cov)

    def get_requirements(self):
        """Return requirements for Cobaya."""
        return {
            "angular_diameter_distance": {
                "z": self.z
            }
        }
    
    def logp(self, **params_values) -> float:
        """Compute log-likelihood for given theory predictions.
        
        Returns:
            Log-likelihood value
        """
        # Get theory prediction
        d_a = self.provider.get_angular_diameter_distance(self.z)
        
        # Check for unphysical distances before proceeding
        if d_a is None or np.any(d_a <= 0.0) or np.any(np.isnan(d_a)):
            return -np.inf
            
        d_l = d_a * (1.0 + self.z)**2
        
        with np.errstate(divide='ignore', invalid='ignore'):
            theory_prediction = 5.0 * np.log10(d_l) + 25.0
            
        if np.any(np.isnan(theory_prediction)) or np.any(np.isinf(theory_prediction)):
            return -np.inf
        
        # We need the nuisance parameter M, but since we are just doing a basic
        # placeholder for the full Pantheon+ implementation here, we'll
        # compute an analytical marginalization over M or just use a simple residual.
        # For simplicity, calculate best-fit M analytically.
        delta = self.mb - theory_prediction
        
        # Handle cases where parameters lead to extreme deltas
        if np.any(np.isnan(delta)) or np.any(np.isinf(delta)):
            print(f"Pantheon Likelihood: delta contains NaN or inf. z_len={len(self.z)}")
            return -np.inf
            
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            # Use precomputed inverse covariance
            M_best = np.sum(self.inv_cov @ delta) / self.sum_inv_cov
            residual = delta - M_best
            chi2 = residual.T @ self.inv_cov @ residual
            
        if np.isnan(chi2) or np.isinf(chi2):
            print(f"Pantheon Likelihood: chi2 is NaN or inf! chi2={chi2}")
            return -np.inf
        return float(-0.5 * chi2)
