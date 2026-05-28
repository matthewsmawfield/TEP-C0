"""Core cosmology classes for TEP pipeline.

Provides standard FLRW and TEP-modulated cosmology classes.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad, quad_vec

C_KMS = 299792.458  # Speed of light in km/s

# Numba JIT compilation for performance optimization on M4 Pro
try:
    from numba import jit, prange
    HAS_NUMBA = True
except ImportError:
    # Fallback if numba not installed
    def jit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator
    def prange(x):
        return range(x)
    HAS_NUMBA = False


# JIT-compiled helper functions for performance
@jit(nopython=True)
def _compute_integrand_grid(zp, H0, Om0, Ok0, Ode0, c):
    """Compute integrand for comoving distance on a grid (JIT-compiled)."""
    zp1 = 1.0 + zp
    return c / (H0 * np.sqrt(Om0 * zp1**3 + Ok0 * zp1**2 + Ode0))


@jit(nopython=True)
def _compute_integrand_grid_wcdm(zp, H0, Om0, Ok0, Ode0, w, c):
    """Compute integrand for comoving distance on a grid for wCDM."""
    zp1 = 1.0 + zp
    return c / (H0 * np.sqrt(Om0 * zp1**3 + Ok0 * zp1**2 + Ode0 * (zp1**(3.0 * (1.0 + w)))))


@jit(nopython=True)
def _compute_integrand_grid_cpl(zp, H0, Om0, Ok0, Ode0, w0, wa, c):
    """Compute integrand for comoving distance on a grid for CPL."""
    zp1 = 1.0 + zp
    return c / (H0 * np.sqrt(Om0 * zp1**3 + Ok0 * zp1**2 + Ode0 * (zp1**(3.0 * (1.0 + w0 + wa))) * np.exp(-3.0 * wa * zp / zp1)))


@jit(nopython=True)
def _compute_cumulative_integral(integrand, dz):
    """Compute cumulative integral using the trapezoidal rule."""
    distances = np.zeros_like(integrand)
    for i in range(1, integrand.size):
        distances[i] = distances[i - 1] + 0.5 * (integrand[i] + integrand[i - 1]) * dz
    return distances


@jit(nopython=True)
def _compute_luminosity_distance(distances, z):
    """Compute luminosity distance from comoving distance (JIT-compiled)."""
    return distances * (1.0 + np.asarray(z))


@jit(nopython=True)
def _compute_distance_modulus(d_l):
    """Compute distance modulus from luminosity distance (JIT-compiled)."""
    return 5.0 * np.log10(d_l) + 25.0


@jit(nopython=True)
def _compute_tep_correction(Sigma_eff, z_arr):
    """Compute TEP temporal shear correction (JIT-compiled)."""
    ln_1pz = np.log(1.0 + z_arr)
    return -2.5 * Sigma_eff * ln_1pz / np.log(10.0)


class CosmologyFLRW:
    """Standard FLRW cosmology for comparison."""
    
    def __init__(self, H0: float = 70.0, Om0: float = 0.3, Ok0: float = 0.0, Ode0: float = None):
        self.H0 = H0
        self.Om0 = Om0
        self.Ok0 = Ok0
        # Ode0 can be specified explicitly or computed from flatness
        if Ode0 is not None:
            self.Ode0 = Ode0
        else:
            self.Ode0 = 1.0 - Om0 - Ok0
        
        # Precompute distance lookup table for performance
        self._distance_cache = {}
        self._cache_key = (H0, Om0, Ok0, self.Ode0)
        
    @jit(nopython=True)
    def e_func(self, z):
        """Dimensionless Hubble parameter (JIT-compiled)."""
        zp1 = 1.0 + np.asarray(z)
        return np.sqrt(self.Om0 * zp1**3 + self.Ok0 * zp1**2 + self.Ode0)
    
    def comoving_distance(self, z):
        """Comoving distance in Mpc using fixed quadrature grid for speed."""
        z_arr = np.atleast_1d(z)
        c = 299792.458  # km/s
        
        # Use fixed quadrature grid for much faster integration
        # Instead of adaptive quad, use Simpson's rule with fixed grid
        n_points = 1000  # Keep low-z Pantheon distances accurate while staying fast.
        max_z = np.max(z_arr)
        
        # Create integration grid
        zp = np.linspace(0, max_z, n_points)
        dz = zp[1] - zp[0]
        
        # Compute integrand at all grid points using JIT-compiled function
        integrand = _compute_integrand_grid(zp, self.H0, self.Om0, self.Ok0, self.Ode0, c)
        
        # Simpson's rule integration using JIT-compiled function
        integral = _compute_cumulative_integral(integrand, dz)
        
        # Interpolate to desired redshifts
        distances = np.interp(z_arr, zp, integral)
        
        return distances if len(distances) > 1 else distances[0]
    
    def angular_diameter_distance(self, z):
        """Angular diameter distance in Mpc."""
        d_c = self.comoving_distance(z)
        return d_c / (1.0 + np.asarray(z))
    
    def luminosity_distance(self, z):
        """Luminosity distance in Mpc using JIT-compiled function."""
        d_c = self.comoving_distance(z)
        return _compute_luminosity_distance(d_c, z)
    
    def distance_modulus(self, z):
        """Distance modulus (m - M) using JIT-compiled function."""
        d_l = self.luminosity_distance(z)
        return _compute_distance_modulus(d_l)


class TEPModulatedCosmology(CosmologyFLRW):
    """TEP-modulated cosmology with temporal shear."""
    
    def __init__(
        self,
        H0: float = 70.0,
        Om0: float = 0.3,
        Sigma_0: float = 0.0,
        A_env: float = 0.1,
        environment: np.ndarray = None,
        Ok0: float = 0.0,
        epsilon_T: float = None,
        z_T: float = 5.0,
    ):
        super().__init__(H0=H0, Om0=Om0, Ok0=Ok0)
        self.Sigma_0 = Sigma_0
        self.A_env = A_env
        self.environment = environment if environment is not None else np.array([0.0])
        self.epsilon_T = epsilon_T if epsilon_T is not None else Sigma_0 * 1000  # Approximate conversion
        self.z_T = z_T
    
    def _effective_sigma(self, z):
        """Effective temporal shear including environment modulation."""
        z_arr = np.atleast_1d(z)
        n_points = len(z_arr)
        
        # Extend environment array if needed
        if len(self.environment) < n_points:
            env = np.tile(self.environment, (n_points // len(self.environment) + 1))[:n_points]
        else:
            env = self.environment[:n_points]
        
        Sigma_eff = self.Sigma_0 * (1.0 + self.A_env * env)
        return Sigma_eff
    
    def distance_modulus(self, z):
        """TEP distance modulus with temporal shear correction (JIT-optimized)."""
        z_arr = np.atleast_1d(z)
        
        # Base LCDM distance modulus
        mu_lcdm = super().distance_modulus(z_arr)
        
        # TEP correction from temporal shear using JIT-compiled function
        Sigma_eff = self._effective_sigma(z_arr)
        correction = _compute_tep_correction(Sigma_eff, z_arr)
        
        return mu_lcdm + correction
    
    def path_enhancement_factor(self, z):
        """TEP path enhancement factor for time dilation predictions.
        
        For TEPModulatedCosmology the temporal shear modifies spatial distances
        but the time dilation factor retains the standard (1+z) form.
        """
        z_arr = np.atleast_1d(z)
        return 1.0 + z_arr


class wCDMCosmology(CosmologyFLRW):
    """wCDM cosmology with constant dark energy equation of state."""
    
    def __init__(self, H0: float = 70.0, Om0: float = 0.3, Ok0: float = 0.0, Ode0: float = None, w: float = -1.0):
        super().__init__(H0=H0, Om0=Om0, Ok0=Ok0, Ode0=Ode0)
        self.w = w
        
    def e_func(self, z):
        zp1 = 1.0 + np.asarray(z)
        return np.sqrt(self.Om0 * zp1**3 + self.Ok0 * zp1**2 + self.Ode0 * (zp1**(3.0 * (1.0 + self.w))))
        
    def comoving_distance(self, z):
        z_arr = np.atleast_1d(z)
        c = 299792.458
        n_points = 1000
        max_z = np.max(z_arr)
        if max_z == 0:
            return np.zeros_like(z_arr)
        zp = np.linspace(0, max_z, n_points)
        dz = zp[1] - zp[0]
        integrand = _compute_integrand_grid_wcdm(zp, self.H0, self.Om0, self.Ok0, self.Ode0, self.w, c)
        integral = _compute_cumulative_integral(integrand, dz)
        distances = np.interp(z_arr, zp, integral)
        return distances if len(distances) > 1 else distances[0]


class CPLCosmology(CosmologyFLRW):
    """CPL (w0 wa) cosmology with evolving dark energy equation of state."""
    
    def __init__(self, H0: float = 70.0, Om0: float = 0.3, Ok0: float = 0.0, Ode0: float = None, w0: float = -1.0, wa: float = 0.0):
        super().__init__(H0=H0, Om0=Om0, Ok0=Ok0, Ode0=Ode0)
        self.w0 = w0
        self.wa = wa
        
    def e_func(self, z):
        zp = np.asarray(z)
        zp1 = 1.0 + zp
        return np.sqrt(self.Om0 * zp1**3 + self.Ok0 * zp1**2 + self.Ode0 * (zp1**(3.0 * (1.0 + self.w0 + self.wa))) * np.exp(-3.0 * self.wa * zp / zp1))
        
    def comoving_distance(self, z):
        z_arr = np.atleast_1d(z)
        c = 299792.458
        n_points = 1000
        max_z = np.max(z_arr)
        if max_z == 0:
            return np.zeros_like(z_arr)
        zp = np.linspace(0, max_z, n_points)
        dz = zp[1] - zp[0]
        integrand = _compute_integrand_grid_cpl(zp, self.H0, self.Om0, self.Ok0, self.Ode0, self.w0, self.wa, c)
        integral = _compute_cumulative_integral(integrand, dz)
        distances = np.interp(z_arr, zp, integral)
        return distances if len(distances) > 1 else distances[0]

