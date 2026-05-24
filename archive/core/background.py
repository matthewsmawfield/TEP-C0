"""TEP-Modified Background Cosmology.

Full implementation of Friedmann equations with TEP modification:
H_TEP(a) = H_LCDM(a) * Gamma_TEP(a)
"""

import numpy as np
from scipy.integrate import quad
from typing import Callable

# Physical constants
C_KMS = 299792.458  # km/s
C_MS = 2.99792458e8  # m/s
G_NEWTON = 6.67430e-11  # m^3 kg^-1 s^-2
K_BOLTZMANN = 1.380649e-23  # J/K
H_PLANCK = 6.62607015e-34  # J s
SIGMA_SB = 5.670374419e-8  # W m^-2 K^-4
M_PROTON = 1.67262192369e-27  # kg
T_CMB = 2.725  # K today


class TEPBackground:
    """Background cosmology with TEP modification.
    
    Implements full Friedmann equations:
        H^2 = (8πG/3) * rho_total
    
    with TEP modification:
        H_TEP(z) = H_LCDM(z) * exp(Sigma_0 * c/H0 * ln(1+z))
    """
    
    def __init__(self, H0: float, Omega_b: float, Omega_cdm: float,
                 Omega_Lambda: float, Omega_k: float = 0.0,
                 T_cmb: float = T_CMB, N_eff: float = 3.046,
                 Sigma_0: float = 0.0):
        """
        Args:
            H0: Hubble constant in km/s/Mpc
            Omega_b: Baryon density parameter
            Omega_cdm: Cold dark matter density parameter
            Omega_Lambda: Dark energy density parameter
            Omega_k: Curvature density parameter
            T_cmb: CMB temperature today (K)
            N_eff: Effective number of neutrino species
            Sigma_0: TEP shear amplitude (0 for standard)
        """
        self.H0 = H0
        self.Omega_b = Omega_b
        self.Omega_cdm = Omega_cdm
        self.Omega_Lambda = Omega_Lambda
        self.Omega_k = Omega_k
        self.T_cmb = T_cmb
        self.N_eff = N_eff
        self.Sigma_0 = Sigma_0
        
        # Critical density today
        H0_SI = H0 * 1000 / 3.08567758e22  # Convert to s^-1
        self.rho_crit = 3 * H0_SI**2 / (8 * np.pi * G_NEWTON)
        
        # Compute radiation density from CMB
        # a = 4*sigma/c = 7.5657e-16 J/m^3/K^4
        # Convert to kg/m^3 by dividing by c^2
        a_rad = 4 * SIGMA_SB / C_MS / C_MS**2  # kg/m^3/K^4
        rho_gamma = a_rad * T_cmb**4  # kg/m^3
        
        # Photon density parameter
        self.Omega_gamma = rho_gamma / self.rho_crit
        
        # Neutrino density parameter
        # rho_nu = N_eff * (7/8) * (4/11)^(4/3) * rho_gamma
        self.Omega_nu = self.N_eff * (7.0/8.0) * (4.0/11.0)**(4.0/3.0) * self.Omega_gamma
        
        # Total radiation
        self.Omega_r = self.Omega_gamma + self.Omega_nu
        
        # Ensure flatness if Omega_k = 0
        if Omega_k == 0.0:
            total = self.Omega_b + self.Omega_cdm + self.Omega_r + self.Omega_Lambda
            if abs(total - 1.0) > 1e-6:
                self.Omega_Lambda = 1.0 - (self.Omega_b + self.Omega_cdm + self.Omega_r)
    
    def tep_gamma(self, z: float) -> float:
        """TEP path enhancement factor Gamma_TEP(z)."""
        if self.Sigma_0 == 0:
            return 1.0
        ln_gamma = self.Sigma_0 * C_KMS / self.H0 * np.log(1.0 + z)
        return np.exp(ln_gamma)
    
    def H(self, z: float) -> float:
        """Hubble parameter in km/s/Mpc with TEP modification."""
        zp1 = 1.0 + z
        
        # Standard LCDM Hubble parameter
        H_std = self.H0 * np.sqrt(
            self.Omega_r * zp1**4 +
            (self.Omega_b + self.Omega_cdm) * zp1**3 +
            self.Omega_k * zp1**2 +
            self.Omega_Lambda
        )
        
        # TEP modification
        return H_std * self.tep_gamma(z)
    
    def H_SI(self, z: float) -> float:
        """Hubble parameter in SI units (s^-1)."""
        return self.H(z) * 1000.0 / 3.08567758e22
    
    def dH_dz(self, z: float) -> float:
        """Derivative of Hubble parameter with respect to z."""
        eps = 1e-8
        return (self.H(z + eps) - self.H(z - eps)) / (2 * eps)
    
    def conformal_H(self, z: float) -> float:
        """Conformal Hubble parameter aH in km/s/Mpc."""
        a = 1.0 / (1.0 + z)
        return self.H(z) * a
    
    def conformal_H_SI(self, z: float) -> float:
        """Conformal Hubble parameter aH in SI (s^-1)."""
        return self.H_SI(z) / (1.0 + z)
    
    def rho_tot(self, z: float) -> float:
        """Total energy density in kg/m^3."""
        H_SI = self.H_SI(z)
        return 3 * H_SI**2 / (8 * np.pi * G_NEWTON)
    
    def equation_of_state(self, z: float) -> float:
        """Total equation of state w = p/rho."""
        zp1 = 1.0 + z
        
        # Energy densities
        rho_r = self.Omega_r * zp1**4 * self.rho_crit
        rho_m = (self.Omega_b + self.Omega_cdm) * zp1**3 * self.rho_crit
        rho_lambda = self.Omega_Lambda * self.rho_crit
        
        # Pressures
        p_r = rho_r / 3.0
        p_m = 0.0
        p_lambda = -rho_lambda
        
        rho_tot = rho_r + rho_m + rho_lambda
        p_tot = p_r + p_m + p_lambda
        
        return p_tot / rho_tot if rho_tot > 0 else 0.0
    
    def comoving_distance(self, z: float) -> float:
        """Comoving distance D_C(z) in Mpc."""
        def integrand(zp):
            return 1.0 / self.H(zp)
        
        integral, _ = quad(integrand, 0, z, limit=200, epsabs=1e-12)
        return C_KMS * integral
    
    def angular_diameter_distance(self, z: float) -> float:
        """Angular diameter distance D_A(z) in Mpc."""
        dc = self.comoving_distance(z)
        
        if abs(self.Omega_k) < 1e-10:
            # Flat universe
            return dc / (1.0 + z)
        elif self.Omega_k > 0:
            # Open universe
            x = np.sqrt(self.Omega_k) * dc * self.H0 / C_KMS
            return C_KMS / (self.H0 * np.sqrt(self.Omega_k) * (1+z)) * np.sinh(x)
        else:
            # Closed universe  
            x = np.sqrt(-self.Omega_k) * dc * self.H0 / C_KMS
            return C_KMS / (self.H0 * np.sqrt(-self.Omega_k) * (1+z)) * np.sin(x)
    
    def luminosity_distance(self, z: float) -> float:
        """Luminosity distance D_L(z) in Mpc."""
        return self.angular_diameter_distance(z) * (1.0 + z)**2
    
    def time_since_big_bang(self, z: float) -> float:
        """Cosmic time since big bang in seconds."""
        def integrand(zp):
            return 1.0 / ((1.0 + zp) * self.H_SI(zp))
        
        integral, _ = quad(integrand, z, np.inf, limit=200, epsabs=1e-12)
        return integral
    
    def sound_speed_baryon(self, z: float) -> float:
        """Baryon sound speed c_s in m/s."""
        a = 1.0 / (1.0 + z)
        R = 0.75 * self.Omega_b / self.Omega_gamma * (1.0 / a)
        return C_MS / np.sqrt(3.0 * (1.0 + R))
    
    def sound_horizon(self, z: float) -> float:
        """Sound horizon at redshift z in Mpc."""
        def integrand(zp):
            cs = self.sound_speed_baryon(zp)
            H = self.H_SI(zp)
            if H <= 0:
                return 0.0
            return cs / H
        
        # Integrate from z to high z (not infinity to avoid numerical issues)
        z_max = max(1e8, z * 1000)
        try:
            integral, err = quad(integrand, z, z_max, limit=200, epsabs=1e-10, epsrel=1e-6)
            if abs(integral) < 1e-30 or err > abs(integral):
                # Fallback: approximate integral
                # In radiation domination: r_s ~ c_s / H at z
                cs = self.sound_speed_baryon(z)
                H = self.H_SI(z)
                integral = cs / H if H > 0 else 0.0
        except (ValueError, RuntimeError, ZeroDivisionError, OverflowError):
            integral = 0.0
        
        return integral / 3.08567758e22  # Convert to Mpc
