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
        """Physical path distance in pure temporal shear model."""
        z_arr = np.atleast_1d(z)
        c = 299792.458
        d = c * np.log1p(z_arr) / self.H0
        return d if len(d) > 1 else d[0]

    def angular_diameter_distance(self, z):
        """Angular diameter distance in static non-expanding metric."""
        # Without spatial expansion, objects do not artificially appear larger
        return self.comoving_distance(z)

    def luminosity_distance(self, z):
        """Luminosity distance in pure temporal shear model.
        
        Assumes photon energy redshifts by (1+z) and arrival rate 
        slows by (1+z), yielding D_L = d * (1+z).
        """
        d = self.comoving_distance(z)
        z_arr = np.atleast_1d(z)
        d_l = d * (1 + z_arr)
        return d_l if len(d_l) > 1 else d_l[0]
        
    def time_dilation_factor(self, z):
        """Predicted SN time dilation factor in pure temporal shear model."""
        z_arr = np.atleast_1d(z)
        # For pure temporal shear, rate of time slows down exactly with redshift
        td = 1.0 + z_arr
        return td if len(td) > 1 else td[0]
        
    def distance_duality_eta(self, z):
        """Distance duality eta parameter D_L = D_A * (1+z)^2 * eta."""
        # D_L = d * (1+z), D_A = d
        # So D_L / (D_A * (1+z)^2) = 1 / (1+z)
        z_arr = np.atleast_1d(z)
        eta = 1.0 / (1.0 + z_arr)
        return eta if len(eta) > 1 else eta[0]
        
    def tolman_exponent(self):
        """Tolman surface brightness exponent alpha, where SB ~ (1+z)^-alpha."""
        # In this static metric, SB decreases by (1+z)^-2
        # (one factor for energy, one for arrival rate, no solid angle effect from expansion)
        return 2.0
        
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
