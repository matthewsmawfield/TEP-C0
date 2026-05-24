#!/usr/bin/env python3
"""Step 038b: Weak Lensing Cosmic Shear Execution.

Calculates the TEP weak lensing cosmic shear predictions using the
Universal Temporal Shear Axiom (unified epsilon_T for distances and growth).

This module replaces the planning step with a concrete implementation of the
TEP weak lensing theoretical predictions.
"""

import numpy as np
from core.tep_cosmology import TEPCosmology
from core.structure_formation import TEPStructureFormation

class TEPWeakLensing:
    def __init__(self, epsilon_T=0.1, z_T=5.0, H0=70.0, Omega_m=0.3):
        """Initialize the Weak Lensing framework with universal parameters."""
        self.epsilon_T = epsilon_T
        self.z_T = z_T
        self.H0 = H0
        self.Omega_m = Omega_m
        
        # Initialize the cosmologies with the exact same epsilon_T
        self.cosmo = TEPCosmology(H0=H0, Omega_m=Omega_m, epsilon_T=epsilon_T, z_T=z_T)
        
        # Structure formation using the same epsilon_T (Universal Shear Axiom)
        self.structure = TEPStructureFormation(
            H0=H0, Omega_m=Omega_m, epsilon_T=epsilon_T, z_T=z_T
        )
        
    def lensing_efficiency(self, z_lens, z_source):
        """Calculate the lensing efficiency kernel g(z_lens, z_source)."""
        if z_lens >= z_source:
            return 0.0
            
        D_L = self.cosmo.comoving_distance(z_lens)
        D_S = self.cosmo.comoving_distance(z_source)
        D_LS = self.cosmo.comoving_distance_between(z_lens, z_source)
        
        # TEP modifications to distance affect the geometric efficiency kernel
        efficiency = D_L * D_LS / D_S
        return efficiency
        
    def convergence_power_spectrum(self, ell, z_sources):
        """Estimate the convergence power spectrum C_ell^kappa at multipole ell."""
        # This is a simplified Limber approximation integration
        c = 299792.458
        H0_cgs = self.H0 * 3.24078e-20  # km/s/Mpc to 1/s
        
        C_ell = np.zeros_like(ell, dtype=float)
        
        # Integrate over lens redshifts
        z_lenses = np.linspace(0.01, max(z_sources)-0.1, 50)
        
        for i, l in enumerate(ell):
            integral = 0.0
            for z_l in z_lenses:
                D_L = self.cosmo.comoving_distance(z_l)
                
                # Wavenumber k for this multipole and redshift
                k = (l + 0.5) / D_L if D_L > 0 else 0.0
                
                if k > 0:
                    # Power spectrum from unified structure formation
                    P_k = self.structure.power_spectrum(k, z=z_l)
                    
                    # Compute mean efficiency for sources
                    eff_mean = np.mean([self.lensing_efficiency(z_l, z_s) for z_s in z_sources if z_s > z_l])
                    
                    # W(z) kernel squared
                    W_z = (1.5 * (self.H0/c)**2 * self.Omega_m * (1+z_l) * eff_mean)**2
                    
                    # Hubble parameter for dz to dr conversion
                    H_z = self.cosmo.H(z_l)
                    
                    # Limber approximation integrand
                    integrand = W_z * P_k * (c / H_z) / (D_L**2)
                    integral += integrand * (z_lenses[1] - z_lenses[0])
                    
            C_ell[i] = integral
            
        return C_ell

def run():
    print("TEP Weak Lensing Execution Module Initialized.")
    
if __name__ == "__main__":
    run()
