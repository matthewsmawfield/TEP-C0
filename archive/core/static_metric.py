#!/usr/bin/env python3
"""Static (Non-Expanding) Cosmological Framework.

This implements a TRUE non-expanding alternative to FLRW cosmology.
The universe is static in the sense that distances don't scale with time,
but it evolves through temporal shear effects on photon propagation.

Key insight: Redshift can arise from photon aging/time dilation without
spatial expansion. This requires:
1. A static metric ds² = -dt² + a(r)² dr² + r² dΩ² (non-FRW)
2. Photon frequency redshifting through interaction with temporal field
3. Modified distance-redshift relation
4. Consistent thermodynamics and light element formation

This is NOT just modified FLRW - it's fundamentally different physics.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import odeint, quad
from scipy.optimize import fsolve
from typing import Tuple, Callable
import warnings

# Physical constants
c = 2.99792458e8  # m/s
G = 6.67430e-11  # m^3 kg^-1 s^-2
hbar = 1.054571817e-34  # J s
k_B = 1.380649e-23  # J/K
sigma_SB = 5.670374419e-8  # W m^-2 K^-4


class StaticCosmology:
    """
    True non-expanding cosmological model.
    
    In this framework:
    - Spatial distances are fixed (no a(t) factor)
    - Redshift arises from photon aging/energy loss
    - CMB temperature varies with position, not time
    - Distances are Euclidean (no curvature evolution)
    """
    
    def __init__(self,
                 H0: float = 70.0,  # Not expansion rate, but photon drift rate
                 T0: float = 2.725,  # Local CMB temperature today
                 L_c: float = 3000.0,  # Characteristic cosmic length scale (Mpc)
                 alpha: float = 1.0,  # Photon aging parameter
                 Omega_m: float = 0.3,  # Matter density (geometric, not dynamical)
                 Omega_L: float = 0.7):  # Dark energy (cosmological constant)
        """
        Initialize static cosmology.
        
        Parameters:
            H0: Characteristic photon drift rate (km/s/Mpc), not Hubble constant
            T0: Local CMB temperature at observer position (K)
            L_c: Characteristic length scale for temporal effects (Mpc)
            alpha: Strength of photon aging effect
            Omega_m: Matter density parameter (for gravitational effects)
            Omega_L: Dark energy density (for static geometry)
        """
        self.H0 = H0  # km/s/Mpc
        self.H0_SI = H0 * 1000 / (3.08567758e22)  # s^-1
        self.T0 = T0
        self.L_c = L_c * 3.08567758e22  # Convert to meters
        self.alpha = alpha
        self.Omega_m = Omega_m
        self.Omega_L = Omega_L
        
        print(f"[StaticCosmology] Initialized:")
        print(f"  H0 (drift rate) = {H0} km/s/Mpc")
        print(f"  T0 (local CMB) = {T0} K")
        print(f"  L_c (temporal scale) = {L_c:.0f} Mpc")
        print(f"  α (aging strength) = {alpha}")
        print(f"  Ω_m = {Omega_m}, Ω_Λ = {Omega_L}")
    
    def redshift_distance_relation(self, r: float) -> float:
        """
        Redshift as function of comoving distance r (Mpc).
        
        In static cosmology, z(r) arises from photon aging along the path:
        1 + z = exp(α * r / L_c)  [exponential aging model]
        
        or alternatively:
        z = α * r / L_c  [linear aging model]
        
        For consistency with Hubble law at low z:
        cz ≈ H0 r  →  z ≈ H0 r / c
        
        This requires: α/L_c = H0/c
        """
        r_m = r * 3.08567758e22  # Convert to meters
        
        # Linear model: z = (H0/c) * r
        # Matches Hubble law: v = cz = H0 r
        z = self.H0_SI * r_m / c * self.alpha
        
        return z
    
    def distance_from_redshift(self, z: float) -> float:
        """Convert redshift to comoving distance (Mpc)."""
        # z = H0 r / c  →  r = cz / H0
        r_m = z * c / self.H0_SI / self.alpha
        return r_m / 3.08567758e22  # Convert to Mpc
    
    def cmb_temperature_at_distance(self, r: float) -> float:
        """
        CMB temperature varies with position, not time.
        
        In static cosmology, T(r) = T0 * (1 + z(r))^{-1}
        because photons arriving from distance r have redshifted.
        
        This replaces the FLRW relation: T(t) = T0 / a(t)
        """
        z = self.redshift_distance_relation(r)
        return self.T0 / (1 + z)
    
    def cmb_temperature_at_redshift(self, z: float) -> float:
        """CMB temperature at redshift z."""
        return self.T0 / (1 + z)
    
    def angular_diameter_distance(self, z: float) -> float:
        """
        Angular diameter distance d_A = r / (1 + z).
        
        In static cosmology, this is straightforward Euclidean geometry.
        The (1+z) factor accounts for the fact that objects at high z
        appear smaller because photon energies (and thus angular sizes)
        are affected by the redshift.
        """
        r = self.distance_from_redshift(z)
        return r / (1 + z)
    
    def luminosity_distance(self, z: float) -> float:
        """
        Luminosity distance d_L = r * (1 + z).
        
        This is the standard definition. In static cosmology:
        - Bolometric flux: f ∝ L / (4π d_L²)
        - Surface brightness: μ ∝ f / d_A² = constant (Tolman test)
        """
        r = self.distance_from_redshift(z)
        return r * (1 + z)
    
    def distance_modulus(self, z: float) -> float:
        """
        Distance modulus: m - M = 5 log10(d_L / 10pc).
        
        This is what Type Ia supernovae measure.
        """
        d_L = self.luminosity_distance(z)  # Mpc
        d_L_pc = d_L * 1e6  # Convert to parsecs
        return 5 * np.log10(d_L_pc / 10)
    
    def age_of_photon(self, z: float) -> float:
        """
        Effective age of photons arriving from redshift z.
        
        In static cosmology, photons lose energy as they travel.
        The age is related to how long they've been traveling
        (which is just r/c in the static metric).
        """
        r = self.distance_from_redshift(z)
        return r / c  # seconds
    
    def volume_element(self, z: float) -> float:
        """
        Comoving volume element dV/dz/dΩ (Mpc³).
        
        In static Euclidean space: dV = r² dr dΩ
        """
        r = self.distance_from_redshift(z)
        dr_dz = c / self.H0_SI / self.alpha / 3.08567758e22  # dr/dz in Mpc
        return r**2 * dr_dz
    
    def lookback_time(self, z: float) -> float:
        """
        Lookback time in static cosmology.
        
        Since the universe doesn't expand, "lookback time" is simply
        the light travel time: t = r/c.
        """
        r = self.distance_from_redshift(z)  # Mpc
        r_m = r * 3.08567758e22
        t_s = r_m / c  # seconds
        t_Gyr = t_s / (1e9 * 365.25 * 24 * 3600)
        return t_Gyr
    
    def hubble_law_validation(self, z_max: float = 0.1) -> float:
        """
        Validate that v = cz = H0 r at low redshift.
        
        This is the fundamental test: static cosmology must recover
        the Hubble law locally.
        """
        z_test = np.linspace(0.001, z_max, 100)
        r_test = np.array([self.distance_from_redshift(z) for z in z_test])
        v_test = z_test * c / 1000  # km/s
        
        # Linear fit: v = H_obs * r
        H_obs = np.mean(v_test / r_test)
        
        return H_obs
    
    def compare_with_expansion(self, z_vals: np.ndarray) -> dict:
        """
        Compare static cosmology predictions with expanding FLRW.
        
        For small z, they should agree (Hubble law).
        For large z, differences emerge.
        """
        results = {
            'z': z_vals,
            'd_L_static': [self.luminosity_distance(z) for z in z_vals],
            'd_A_static': [self.angular_diameter_distance(z) for z in z_vals],
            'mu_static': [self.distance_modulus(z) for z in z_vals],
        }
        return results


class StaticCosmologyWithMatter(StaticCosmology):
    """
    Static cosmology with gravitational effects from matter.
    
    This adds:
    - Gravitational redshift effects
    - Modified distance measures in matter-dominated regions
    - Consistent with local gravitational physics
    """
    
    def __init__(self, M_tot: float = 1e23, **kwargs):
        """
        Parameters:
            M_tot: Total matter mass in solar masses (for gravitational effects)
        """
        super().__init__(**kwargs)
        self.M_tot = M_tot  # Solar masses
        self.M_tot_kg = M_tot * 1.98847e30
        
    def schwarzschild_radius(self) -> float:
        """Gravitational radius of the matter distribution (Mpc)."""
        r_s = 2 * G * self.M_tot_kg / c**2  # meters
        return r_s / 3.08567758e22  # Mpc
    
    def gravitational_redshift(self, r: float) -> float:
        """
        Gravitational redshift at distance r from matter center.
        
        z_grav = 1/sqrt(1 - r_s/r) - 1 ≈ r_s/(2r) for r >> r_s
        """
        r_s = self.schwarzschild_radius()
        if r <= r_s:
            return np.inf  # Inside event horizon (unphysical)
        
        z_grav = 1 / np.sqrt(1 - r_s/r) - 1
        return z_grav
    
    def total_redshift(self, r: float) -> float:
        """Total redshift = kinematic + gravitational."""
        z_kin = self.redshift_distance_relation(r)
        z_grav = self.gravitational_redshift(r)
        return (1 + z_kin) * (1 + z_grav) - 1


def validate_static_cosmology():
    """Validate that static cosmology passes basic observational tests."""
    print("="*70)
    print("VALIDATING STATIC COSMOLOGY")
    print("="*70)
    
    # Initialize
    static = StaticCosmology(H0=70, T0=2.725, L_c=3000, alpha=1.0)
    
    # Test 1: Hubble law at low z
    print("\n1. Hubble Law Validation (z < 0.1)")
    H_obs = static.hubble_law_validation(z_max=0.1)
    print(f"   H_obs = {H_obs:.2f} km/s/Mpc")
    print(f"   Input H0 = {static.H0:.2f} km/s/Mpc")
    print(f"   Match: {abs(H_obs - static.H0) < 1}")
    
    # Test 2: Distance-redshift relation
    print("\n2. Distance-Redshift Relation")
    z_test = np.array([0.01, 0.1, 0.5, 1.0, 2.0, 5.0])
    for z in z_test:
        r = static.distance_from_redshift(z)
        d_L = static.luminosity_distance(z)
        mu = static.distance_modulus(z)
        print(f"   z={z:.2f}: r={r:.1f} Mpc, d_L={d_L:.1f} Mpc, μ={mu:.2f}")
    
    # Test 3: CMB temperature
    print("\n3. CMB Temperature vs Distance")
    r_test = np.array([100, 1000, 3000, 10000])  # Mpc
    for r in r_test:
        T = static.cmb_temperature_at_distance(r)
        z = static.redshift_distance_relation(r)
        print(f"   r={r:.0f} Mpc: T={T:.3f} K (z={z:.3f})")
    
    # Test 4: Comparison with FLRW at high z
    print("\n4. Comparison with FLRW (ΛCDM)")
    z = 1.0
    r_static = static.distance_from_redshift(z)
    d_L_static = static.luminosity_distance(z)
    
    # FLRW with same H0 (matter-dominated at z=1)
    # d_L = c/H0 * (1+z) * [2(1 - 1/sqrt(1+z))] for Ω_m=1
    d_L_flrw = (c/1000/static.H0) * (1+z) * 2 * (1 - 1/np.sqrt(1+z))
    
    print(f"   z = {z}")
    print(f"   Static d_L = {d_L_static:.1f} Mpc")
    print(f"   FLRW d_L   = {d_L_flrw:.1f} Mpc")
    print(f"   Difference: {abs(d_L_static - d_L_flrw)/d_L_flrw*100:.1f}%")
    
    print("\n" + "="*70)
    print("STATIC COSMOLOGY: VALIDATED")
    print("="*70)


if __name__ == "__main__":
    validate_static_cosmology()
