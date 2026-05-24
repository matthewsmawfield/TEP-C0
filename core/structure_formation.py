"""Structure formation calculations for TEP cosmology."""

from __future__ import annotations

import numpy as np
from scipy.integrate import odeint


class StructureFormation:
    """Structure formation in TEP cosmology."""
    
    def __init__(
        self,
        H0: float = 70.0,
        Omega_m: float = 0.3,
        Omega_L: float = 0.7,
        Omega_b: float = 0.045,
        Omega_cdm: float = None,
        sigma_8: float = 0.81,
        n_s: float = 0.965,
        Sigma_0: float = 0.0,
        Om0: float = None,  # Alias for Omega_m (backward compatibility)
    ):
        self.H0 = H0
        # Handle both Omega_m and Om0 for backward compatibility
        self.Omega_m = Omega_m if Om0 is None else Om0
        self.Omega_L = Omega_L
        self.Omega_b = Omega_b
        self.Omega_cdm = Omega_cdm if Omega_cdm is not None else (self.Omega_m - Omega_b)
        self.sigma_8 = sigma_8
        self.n_s = n_s
        self.Sigma_0 = Sigma_0
        self.h = H0 / 100.0
        # Alias for backward compatibility
        self.Om0 = self.Omega_m
        
    def growth_factor(self, a):
        """Linear growth factor D(a). Alias for compute_growth_factor."""
        return self.compute_growth_factor(a)
    
    def compute_growth_factor(self, a):
        """Linear growth factor D(a).
        
        In TEP, growth is modified by temporal shear.
        For LCDM limit (Sigma_0=0), uses Carroll et al. 1992 approximation.
        """
        a_arr = np.atleast_1d(a)
        
        # LCDM growth factor (Carroll et al. 1992)
        def D_lcdm(a_val):
            if a_val <= 0:
                return 0.0
            z = 1.0 / a_val - 1.0
            Omega_mz = self.Omega_m * (1 + z)**3 / (self.Omega_m * (1 + z)**3 + self.Omega_L)
            return a_val * 2.5 * Omega_mz / (Omega_mz**(4.0/7.0) - (1 - Omega_mz) + (1 + 0.5*Omega_mz)*(1 + (1 - Omega_mz)/70.0))
        
        # TEP correction: growth suppression from temporal shear
        D_values = np.array([D_lcdm(ai) for ai in a_arr])
        
        if self.Sigma_0 > 0:
            # Temporal shear suppresses growth
            suppression = 1.0 - 0.1 * self.Sigma_0 * 1e4 * np.log(1.0 / np.maximum(a_arr, 0.001))
            D_values = D_values * np.maximum(suppression, 0.5)
        
        return D_values if len(D_values) > 1 else float(D_values[0])
    
    def growth_rate(self, a):
        """Growth rate f(a) = d ln D / d ln a."""
        a_val = float(np.atleast_1d(a)[0])
        da = 0.001
        D_plus = self.compute_growth_factor(a_val + da)
        D_minus = self.compute_growth_factor(max(a_val - da, 0.001))
        D_center = self.compute_growth_factor(a_val)
        
        if D_center <= 0 or D_plus <= 0 or D_minus <= 0:
            return 0.55
            
        d_ln_D = np.log(D_plus / D_minus) / 2
        d_ln_a = np.log((a_val + da) / max(a_val - da, 0.001)) / 2
        
        return d_ln_D / d_ln_a
    
    def fsigma8(self, z):
        """f*sigma_8(z) for RSD comparisons."""
        a = 1.0 / (1.0 + z)
        f = self.growth_rate(a)
        D = self.compute_growth_factor(a)
        sigma_8_z = self.sigma_8 * D
        return f * sigma_8_z
    
    def power_spectrum(self, k, z=0):
        """Matter power spectrum P(k) at redshift z.
        
        Simplified Eisenstein & Hu 1998 approximation.
        """
        # Transfer function (simplified)
        Gamma = self.Omega_m * self.h
        q = k / Gamma
        T = np.log(1 + 2.34*q) / (2.34*q) / (1 + 3.89*q + (16.1*q)**2 + (5.46*q)**3 + (6.71*q)**4)**0.25
        
        # Growth factor
        a = 1.0 / (1.0 + z)
        D = self.compute_growth_factor(a)
        
        # Normalization from sigma_8
        # P(k) ~ k^n_s * T(k)^2 * D(z)^2
        P = k**self.n_s * T**2 * D**2
        
        # Normalize to sigma_8 = 0.81 at z=0
        # This is a simplified normalization
        P_norm = P * (self.sigma_8 / 0.81)**2
        
        return P_norm


class TEPStructureFormation(StructureFormation):
    """Structure formation in TEP cosmology with temporal shear effects.
    
    Extends StructureFormation with TEP-specific modifications:
    - Temporal shear suppression of growth
    - Modified growth rate from TEP transport kernel
    - Sound horizon shift from modified expansion history
    """
    
    def __init__(
        self,
        H0: float = 70.0,
        Omega_m: float = 0.3,
        Omega_L: float = 0.7,
        Omega_b: float = 0.045,
        Omega_cdm: float = None,
        sigma_8: float = 0.81,
        n_s: float = 0.965,
        Sigma_0: float = 0.0,
        epsilon_T: float = 0.0,
        z_T: float = 5.0,
        Om0: float = None,  # Alias for Omega_m (backward compatibility)
    ):
        super().__init__(
            H0=H0,
            Omega_m=Omega_m,
            Omega_L=Omega_L,
            Omega_b=Omega_b,
            Omega_cdm=Omega_cdm,
            sigma_8=sigma_8,
            n_s=n_s,
            Sigma_0=Sigma_0,
            Om0=Om0,
        )
        self.epsilon_T = epsilon_T
        self.z_T = z_T
    
    def sound_horizon(self) -> float:
        """Sound horizon at drag epoch in TEP cosmology.
        
        TEP modifies expansion history, which changes the sound horizon.
        For epsilon_T > 0, the sound horizon is reduced due to faster expansion.
        """
        # Standard LCDM sound horizon (Planck 2018 reference)
        r_drag_lcdm = 147.09  # Mpc
        
        # TEP correction: expansion faster at high z for epsilon_T > 0
        # This reduces the sound horizon
        if self.epsilon_T > 0:
            # Approximate 31% reduction for typical TEP parameters
            reduction_factor = 1.0 - 0.31 * (self.epsilon_T / 0.23)
            r_drag_tep = r_drag_lcdm * max(reduction_factor, 0.5)
        else:
            r_drag_tep = r_drag_lcdm
        
        return r_drag_tep
    
    def sigma_8(self) -> float:
        """sigma_8 in TEP cosmology.
        
        TEP temporal shear suppresses growth, reducing sigma_8.
        """
        # Base sigma_8 from parent class
        sigma_8_base = self.sigma_8
        
        # TEP suppression from temporal shear
        if self.Sigma_0 > 0:
            suppression = 1.0 - 0.1 * self.Sigma_0 * 1e4
            sigma_8_tep = sigma_8_base * max(suppression, 0.5)
        else:
            sigma_8_tep = sigma_8_base
        
        return sigma_8_tep
