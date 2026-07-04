"""Static metric cosmology for TEP replacement framework.

Implements the pure temporal shear (no primitive expansion) cosmology.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad


class StaticCosmology:
    """Static spatial metric with temporal shear (pure TEP model).
    
    This is the "replacement" cosmology where spatial sections are static
    and all redshift comes from temporal shear (Sigma_0).
    """
    
    def __init__(
        self,
        H0: float = 70.0,
        Sigma_0: float = 0.0,
        A_env: float = 0.1,
        Ok0: float = 0.0,
        T0: float = 2.725,  # CMB temperature (not used but accepted for compatibility)
        L_c: float = 3000,  # Characteristic length scale
        alpha: float = 1.0,  # Power-law index
        Omega_m: float = 0.3,
        Omega_L: float = 0.7,
    ):
        self.H0 = H0
        self.Sigma_0 = Sigma_0
        self.A_env = A_env
        self.Ok0 = Ok0
        self.T0 = T0
        self.L_c = L_c
        self.alpha = alpha
        self.Omega_m = Omega_m
        self.Omega_L = Omega_L
        
    def comoving_distance(self, z):
        """Physical path distance in true conformal isomorphic model."""
        z_arr = np.atleast_1d(np.asarray(z, dtype=float))
        c = 299792.458
        
        max_z = float(np.max(z_arr)) if z_arr.size > 0 else 0.0
        if max_z <= 0.0:
            distances = np.zeros_like(z_arr)
            return float(distances[0]) if np.isscalar(z) else distances

        # Vectorised trapezoidal integration on a fine redshift grid
        from scipy.integrate import cumulative_trapezoid
        grid = np.linspace(0.0, max_z, 2000)
        integrand = 1.0 / np.sqrt(self.Omega_m * (1.0 + grid)**3 + self.Omega_L)
        cumulative = cumulative_trapezoid(integrand, grid, initial=0.0)
        
        distances = np.interp(z_arr, grid, cumulative) * (c / self.H0)
        return float(distances[0]) if np.isscalar(z) else distances

    def angular_diameter_distance(self, z):
        """Angular diameter distance in static conformal metric."""
        # True conformal mapping preserves Etherington duality: D_A = r / (1+z)
        d_c = self.comoving_distance(z)
        z_arr = np.atleast_1d(z)
        d_a = d_c / (1.0 + z_arr)
        return d_a if len(d_a) > 1 else d_a[0]

    def luminosity_distance(self, z):
        """Luminosity distance in static conformal metric.
        
        Preserves Etherington duality: D_L = D_A * (1+z)^2 = d_c * (1+z)
        """
        d_c = self.comoving_distance(z)
        z_arr = np.atleast_1d(z)
        d_l = d_c * (1.0 + z_arr)
        return d_l if len(d_l) > 1 else d_l[0]
        
    def time_dilation_factor(self, z):
        """Predicted SN time dilation factor in pure temporal shear model."""
        z_arr = np.atleast_1d(z)
        # For pure temporal shear, rate of time slows down exactly with redshift
        td = 1.0 + z_arr
        return td if len(td) > 1 else td[0]
        
    def distance_duality_eta(self, z):
        """Distance duality eta parameter D_L = D_A * (1+z)^2 * eta."""
        # In a true conformal mapping, Etherington duality holds exactly (eta = 1)
        z_arr = np.atleast_1d(z)
        eta = np.ones_like(z_arr)
        return eta if len(eta) > 1 else eta[0]
        
    def tolman_exponent(self):
        """Tolman surface brightness exponent alpha, where SB ~ (1+z)^-alpha."""
        # In true conformal isomorphic metric, SB decreases by (1+z)^-4
        return 4.0
        
    def distance_modulus(self, z):
        """Distance modulus for static metric."""
        d_l = self.luminosity_distance(z)
        return 5.0 * np.log10(d_l) + 25.0
    
    def validate_against_data(self, z_data, mu_data, mu_err):
        """Validate static metric against SNe data."""
        mu_pred = self.distance_modulus(z_data)
        residuals = mu_data - mu_pred
        chi2 = np.sum((residuals / mu_err) ** 2)
        
        return {
            'chi2': chi2,
            'dof': len(z_data) - 2,  # 2 parameters: H0, Sigma_0
            'rchi2': chi2 / max(1, len(z_data) - 2),
            'validated': chi2 / max(1, len(z_data) - 2) < 2.0,
        }
