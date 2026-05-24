"""Core cosmological calculations with proper FLRW and TEP modifications."""

import numpy as np
from scipy import integrate
from typing import Optional, Tuple

C_KMS = 299792.458  # Speed of light in km/s
T_CMB = 2.725  # CMB temperature in K
N_EFF = 3.046  # Effective neutrino species


class CosmologyFLRW:
    """Proper FLRW cosmology with full radiation, matter, curvature, dark energy.
    
    Uses numerical integration for accurate luminosity distances.
    Includes proper radiation density from CMB temperature.
    """
    
    def __init__(self, H0: float, Om0: float, Ode0: float = None, 
                 Ok0: float = 0.0, Tcmb0: float = T_CMB, Neff: float = N_EFF):
        """Initialize FLRW cosmology.
        
        Args:
            H0: Hubble constant in km/s/Mpc
            Om0: Matter density parameter today
            Ode0: Dark energy density (default: 1 - Om0 - Ok0)
            Ok0: Curvature density (default: 0, flat)
            Tcmb0: CMB temperature in K
            Neff: Effective number of neutrino species
        """
        self.H0 = float(H0)
        self.Om0 = float(Om0)
        self.Ok0 = float(Ok0)
        self.Tcmb0 = float(Tcmb0)
        self.Neff = float(Neff)
        
        # Compute radiation density from CMB temperature
        # rho_gamma = a_R T^4, where a_R = 4*sigma_SB/c
        # sigma_SB = 5.670374419e-8 W/m^2/K^4
        sigma_sb = 5.670374419e-8  # Stefan-Boltzmann
        c = 299792458.0  # m/s
        a_rad = 4.0 * sigma_sb / c  # Radiation constant
        
        # Critical density today in kg/m^3
        G = 6.67430e-11  # m^3/kg/s^2
        H0_SI = H0 * 1000.0 / (3.08567758e22)  # Convert km/s/Mpc to s^-1
        rho_crit = 3.0 * H0_SI**2 / (8.0 * np.pi * G)
        
        # Radiation energy density (photons)
        rho_gamma = a_rad * Tcmb0**4  # J/m^3 = kg/m/s^2
        
        # Neutrinos: rho_nu = Neff * (7/8) * (4/11)^(4/3) * rho_gamma
        rho_nu = Neff * (7.0/8.0) * (4.0/11.0)**(4.0/3.0) * rho_gamma
        
        self.Or0 = (rho_gamma + rho_nu) / (rho_crit * c**2)  # Convert to density parameter
        
        # Dark energy
        if Ode0 is None:
            self.Ode0 = 1.0 - self.Om0 - self.Ok0 - self.Or0
        else:
            self.Ode0 = float(Ode0)
        
        # Ensure flatness if Ok0=0
        if Ok0 == 0.0 and abs(self.Om0 + self.Or0 + self.Ode0 - 1.0) > 0.01:
            self.Ode0 = 1.0 - self.Om0 - self.Or0
    
    def e_func(self, z: float) -> float:
        """E(z) = H(z)/H0 = sqrt(Om(1+z)^3 + Or(1+z)^4 + Ok(1+z)^2 + Ode)."""
        zp1 = 1.0 + z
        return np.sqrt(
            self.Om0 * zp1**3 + 
            self.Or0 * zp1**4 + 
            self.Ok0 * zp1**2 + 
            self.Ode0
        )
    
    def e_func_inv(self, z: float) -> float:
        """1/E(z) for integration."""
        return 1.0 / self.e_func(z)
    
    def comoving_distance(self, z: float) -> float:
        """Comoving distance D_C(z) in Mpc.
        
        D_C = c/H0 * integral_0^z dz'/E(z')
        """
        if z <= 0:
            return 0.0
        
        integral, _ = integrate.quad(self.e_func_inv, 0, z, limit=200, epsabs=1e-12)
        return C_KMS / self.H0 * integral
    
    def luminosity_distance(self, z: np.ndarray) -> np.ndarray:
        """Luminosity distance D_L(z) in Mpc.
        
        D_L = (1+z) * D_C for flat universe
        D_L = (1+z)/sqrt(|Ok|) * sinh(sqrt(|Ok|)*D_C*H0/c) for curved
        """
        z_arr = np.atleast_1d(z)
        dc = np.array([self.comoving_distance(zi) for zi in z_arr])
        
        if abs(self.Ok0) < 1e-10:
            # Flat universe
            return (1.0 + z_arr) * dc
        elif self.Ok0 > 0:
            # Open universe
            x = np.sqrt(self.Ok0) * dc * self.H0 / C_KMS
            return (1.0 + z_arr) * C_KMS / (self.H0 * np.sqrt(self.Ok0)) * np.sinh(x)
        else:
            # Closed universe
            x = np.sqrt(-self.Ok0) * dc * self.H0 / C_KMS
            return (1.0 + z_arr) * C_KMS / (self.H0 * np.sqrt(-self.Ok0)) * np.sin(x)
    
    def angular_diameter_distance(self, z: np.ndarray) -> np.ndarray:
        """Angular diameter distance D_A(z) in Mpc.
        
        D_A = D_L / (1+z)^2
        """
        return self.luminosity_distance(z) / (1.0 + z)**2
    
    def distance_modulus(self, z: np.ndarray) -> np.ndarray:
        """Distance modulus mu = 5*log10(D_L) + 25."""
        dl = self.luminosity_distance(z)
        # Protect against log of negative or zero
        dl = np.maximum(dl, 1e-10)
        return 5.0 * np.log10(dl) + 25.0


class TEPModulatedCosmology(CosmologyFLRW):
    """TEP-modulated cosmology with path enhancement factor.
    
    The path enhancement factor Gamma_TEP modifies the effective distance:
    D_L^TEP = D_L^FLRW / Gamma_TEP
    
    where Gamma_TEP = exp(Sigma_0 * c/H0 * ln(1+z))
    """
    
    def __init__(self, H0: float, Om0: float, Sigma_0: float, 
                 A_env: float = 0.1, environment: np.ndarray = None,
                 **kwargs):
        """Initialize TEP cosmology.
        
        Args:
            H0: Hubble constant
            Om0: Matter density
            Sigma_0: Base temporal shear amplitude
            A_env: Environment coupling strength
            environment: Environment array (normalized)
            **kwargs: Passed to CosmologyFLRW
        """
        super().__init__(H0, Om0, **kwargs)
        self.Sigma_0 = float(Sigma_0)
        self.A_env = float(A_env)
        self.environment = environment
    
    def path_enhancement_factor(self, z: np.ndarray) -> np.ndarray:
        """Compute Gamma_TEP = exp(Sigma_eff * c/H0 * ln(1+z))."""
        z_arr = np.atleast_1d(z)
        
        # Effective shear
        if self.environment is not None:
            Sigma_eff = self.Sigma_0 * (1.0 + self.A_env * self.environment)
        else:
            Sigma_eff = self.Sigma_0
        
        # Ensure positive for log
        z_safe = np.maximum(z_arr, 1e-10)
        c_over_H0 = C_KMS / self.H0
        ln_gamma = Sigma_eff * c_over_H0 * np.log(1.0 + z_safe)
        return np.exp(ln_gamma)
    
    def luminosity_distance(self, z: np.ndarray) -> np.ndarray:
        """TEP luminosity distance with path enhancement.

        D_L^TEP = D_L^FLRW / Gamma_TEP

        The path enhancement factor reduces the observed flux because
        temporal shear spreads photons across more path segments.
        """
        dl_lcdm = super().luminosity_distance(z)
        gamma = self.path_enhancement_factor(z)
        return dl_lcdm / gamma

    def angular_diameter_distance(self, z: np.ndarray) -> np.ndarray:
        """TEP angular diameter distance — geometric, unmodified by shear.

        D_A^TEP = D_A^FLRW (angular size set by geometry alone)

        The temporal shear affects flux (D_L) but not angular size (D_A),
        because angular size is a purely geometric quantity determined by
        the spatial metric, not by photon transport efficiency.

        This breaks the Etherington relation:
        η(z) = D_L^TEP / (D_A^TEP * (1+z)²) = 1 / Gamma_TEP(z)
        """
        z_arr = np.atleast_1d(z)
        dc = np.array([self.comoving_distance(zi) for zi in z_arr])
        if abs(self.Ok0) < 1e-10:
            return dc / (1.0 + z_arr)
        elif self.Ok0 > 0:
            x = np.sqrt(self.Ok0) * dc * self.H0 / C_KMS
            return C_KMS / (self.H0 * np.sqrt(self.Ok0) * (1 + z_arr)) * np.sinh(x)
        else:
            x = np.sqrt(-self.Ok0) * dc * self.H0 / C_KMS
            return C_KMS / (self.H0 * np.sqrt(-self.Ok0) * (1 + z_arr)) * np.sin(x)
    
    def distance_modulus(self, z: np.ndarray) -> np.ndarray:
        """TEP distance modulus."""
        dl = self.luminosity_distance(z)
        dl = np.maximum(dl, 1e-10)
        return 5.0 * np.log10(dl) + 25.0
