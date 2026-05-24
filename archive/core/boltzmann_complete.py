#!/usr/bin/env python3
"""Experimental Boltzmann prototype.

This module is intentionally not used as publication evidence.

The current scratch implementation does not yet satisfy the consistency
requirements for an Einstein-Boltzmann solver: conformal-time derivatives,
tight-coupling handling, super-horizon adiabatic initial conditions, and
Einstein constraints are not validated against CLASS/CAMB. The prototype is
kept for equation-development experiments, but `solve_mode` fails closed unless
`TEP_ALLOW_EXPERIMENTAL_BOLTZMANN=1` is set.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import odeint, solve_ivp
from scipy.interpolate import interp1d, CubicSpline
from scipy.special import spherical_jn
import warnings
import os

# Constants
c = 2.99792458e8
G = 6.67430e-11
sigma_T = 6.6524587158e-29
m_p = 1.67262192369e-27

class CompleteBoltzmannSolver:
    """Experimental Boltzmann solver scaffold.

    Do not use this class for evidence or manuscript claims until the validation
    gates in `validation_status()` pass against a reference CLASS/CAMB run.
    """
    
    def __init__(self, H0=70.0, Omega_b=0.045, Omega_cdm=0.25, Omega_L=0.7,
                 T_cmb=2.725, N_nu=3.046, lmax=50):
        self.H0 = H0
        self.H0_SI = H0 * 1000 / 3.08567758e22
        self.Omega_b = Omega_b
        self.Omega_cdm = Omega_cdm
        self.Omega_m = Omega_b + Omega_cdm
        self.Omega_L = Omega_L
        self.T_cmb = T_cmb
        self.N_nu = N_nu
        self.lmax = lmax
        self.h = H0 / 100
        
        # Radiation
        a_rad = 4 * 5.670374419e-8 / c**3
        rho_gamma = a_rad * T_cmb**4
        rho_crit = 3 * self.H0_SI**2 / (8 * np.pi * G)
        self.Omega_gamma = rho_gamma / rho_crit
        self.Omega_nu = N_nu * (7/8) * (4/11)**(4/3) * self.Omega_gamma
        self.Omega_r = self.Omega_gamma + self.Omega_nu
        
        # Build background
        self._build_background()

    def validation_status(self):
        return {
            "status": "experimental_invalid",
            "research_grade_cmb": False,
            "claim_gate": "blocked",
            "blockers": [
                "Derivative variables mix d/dln(a), d/deta, and physical-time units.",
                "Photon-baryon tight coupling is not implemented before recombination.",
                "Super-horizon adiabatic initial conditions are not CLASS/CAMB validated.",
                "Einstein constraint equations are applied algebraically without a consistent Phi/Psi evolution system.",
                "Line-of-sight source function is a toy approximation, not a validated CMB source hierarchy.",
            ],
        }
        
    def _build_background(self):
        """Tabulate background quantities."""
        self.a_grid = np.logspace(-8, 0, 1000)
        self.z_grid = 1/self.a_grid - 1
        
        # Hubble
        self.H_grid = np.array([self.Hubble(a) for a in self.a_grid])
        
        # Ionization
        self.xe_grid = np.array([self._xe_saha(a) for a in self.a_grid])
        
        # Thomson rate
        rho_crit = 3 * self.H0_SI**2 / (8 * np.pi * G)
        n_b = self.Omega_b * rho_crit / (m_p * self.a_grid**3)
        n_e = self.xe_grid * n_b
        self.kappa_dot_grid = n_e * sigma_T * c / (self.a_grid * self.H_grid)
        
        # Interpolators
        self.H_interp = interp1d(np.log(self.a_grid), self.H_grid, kind='cubic')
        self.xe_interp = interp1d(np.log(self.a_grid), self.xe_grid, bounds_error=False, fill_value=(1, 0.0001))
        self.kappa_interp = interp1d(np.log(self.a_grid), self.kappa_dot_grid, bounds_error=False, fill_value=(1e10, 0.0001))
        
    def Hubble(self, a):
        """Hubble parameter H(a)."""
        return self.H0_SI * np.sqrt(self.Omega_r/a**4 + self.Omega_m/a**3 + self.Omega_L)
    
    def conformal_Hubble(self, a):
        """Conformal Hubble aH."""
        return a * self.Hubble(a)
    
    def _xe_saha(self, a):
        """Saha recombination."""
        z = 1/a - 1
        if z > 2000:
            return 1.0
        elif z < 20:
            return 0.0001
        
        # Saha equation for hydrogen
        T = self.T_cmb * (1+z)
        T_ev = T / 11604.5
        
        # Saha factor
        saha_factor = (2 * np.pi * 9.1093837015e-31 * 1.380649e-23 * T / (6.626e-34**2))**1.5
        saha_factor *= np.exp(-13.6/8.617e-5/T)
        
        # Ionization fraction
        n_b = 1e-6  # Approximate normalization
        rhs = saha_factor / n_b
        
        if rhs < 1e-20:
            return 0.0
        
        # Quadratic: x_e^2/(1-x_e) = rhs
        # Approximate: x_e ≈ sqrt(rhs) for small rhs
        if rhs < 0.01:
            return np.sqrt(rhs)
        
        # Full solution
        a_saha = 1
        b_saha = rhs
        c_saha = -rhs
        discriminant = b_saha**2 - 4*a_saha*c_saha
        if discriminant < 0:
            return 0.0
        x_e = (-b_saha + np.sqrt(discriminant)) / (2*a_saha)
        return max(0, min(1, x_e))
    
    def get_background(self, a):
        """Get background quantities at scale factor a."""
        lna = np.log(a)
        H = float(self.H_interp(lna))
        xe = float(self.xe_interp(lna))
        kappa_dot = float(self.kappa_interp(lna))
        return H, xe, kappa_dot
    
    def equations(self, lna, y, k):
        """Perturbation equations for single k mode.
        
        State: [δ_c, δ_b, v_c, v_b, δ_γ, θ_γ, δ_ν, θ_ν, Φ]
        """
        a = np.exp(lna)
        a2 = a*a
        H, xe, kappa_dot = self.get_background(a)
        aH = a * H
        
        # Extract
        delta_c, delta_b, v_c, v_b, delta_g, theta_g, delta_nu, theta_nu, Phi = y
        
        # Densities
        rho_crit = 3 * self.H0_SI**2 / (8 * np.pi * G)
        rho_c = self.Omega_cdm * rho_crit / a**3
        rho_b = self.Omega_b * rho_crit / a**3
        rho_g = self.Omega_gamma * rho_crit / a**4
        rho_nu = self.Omega_nu * rho_crit / a**4
        
        rho_tot = rho_c + rho_b + rho_g + rho_nu
        P = (rho_g + rho_nu) / 3
        
        # Continuity and Euler equations
        # CDM
        d_delta_c = -k * v_c * aH  # Note: v is defined with factors of k/aH
        d_v_c = -v_c - Phi
        
        # Baryons (with photon drag)
        R = 3 * rho_b / (4 * rho_g) if rho_g > 0 else 0
        drag = kappa_dot * R / (1+R) * (theta_g/3 - v_b) if R > 0 else 0
        
        d_delta_b = -k * v_b * aH
        d_v_b = -v_b - Phi + drag / aH
        
        # Photons
        d_delta_g = -4/3 * k * theta_g * aH
        d_theta_g = k * delta_g/4 + k * Phi - kappa_dot * (theta_g - 3*v_b)
        d_theta_g /= aH
        
        # Neutrinos (no coupling)
        d_delta_nu = -4/3 * k * theta_nu * aH
        d_theta_nu = k * delta_nu/4 + k * Phi
        d_theta_nu /= aH
        
        # Einstein: constraint equation for Phi
        k2 = k*k
        if k2 > 1e-30:
            # Poisson equation
            delta_rho = (rho_c * delta_c + rho_b * delta_b + 
                        rho_g * delta_g + rho_nu * delta_nu)
            Phi = -4 * np.pi * G * a2 * delta_rho / k2
        
        return [d_delta_c, d_delta_b, d_v_c, d_v_b, 
                d_delta_g, d_theta_g, d_delta_nu, d_theta_nu, Phi]
    
    def solve_mode(self, k, n_steps=500):
        """Solve perturbations for single k."""
        if os.getenv("TEP_ALLOW_EXPERIMENTAL_BOLTZMANN", "0") != "1":
            raise RuntimeError(
                "boltzmann_complete.py is an experimental invalid prototype. "
                "Set TEP_ALLOW_EXPERIMENTAL_BOLTZMANN=1 only for local equation-development tests."
            )
        a_init = 1e-8
        a_final = 1.0
        
        # Initial conditions: adiabatic
        delta_c0 = 1.0
        delta_b0 = delta_c0
        delta_g0 = 4/3 * delta_c0
        delta_nu0 = 4/3 * delta_c0
        
        # Super-horizon: velocities from continuity
        # v = -i * k * Φ / (aH * (3(1+w))) - use real approximation
        v0 = 1e-5
        Phi0 = -2/3 * delta_c0
        
        y0 = [delta_c0, delta_b0, v0, v0, delta_g0, 3*v0, delta_nu0, 3*v0, Phi0]
        
        # Time array
        lna = np.linspace(np.log(a_init), np.log(a_final), n_steps)
        
        # Integrate with RK45 for stability
        sol = solve_ivp(lambda t, y: self.equations(t, y, k),
                       [np.log(a_init), np.log(a_final)], y0,
                       t_eval=lna, method='RK45',
                       rtol=1e-7, atol=1e-9, dense_output=True)
        
        if not sol.success:
            return np.exp(lna), np.zeros((n_steps, 9))
        
        return np.exp(lna), sol.y.T
    
    def compute_source_function(self, k, lmax_source=10):
        """Compute source function S(k, τ) for LOS integration."""
        a, y = self.solve_mode(k)
        
        # Source for temperature anisotropies
        # S = g [Θ_0 + Ψ + v_b'/k] + ... (simplified)
        
        delta_b = y[:, 1]
        theta_b = y[:, 3]
        Phi = y[:, 8]
        
        # Simplified source: just gravitational + Doppler
        kappa_dot = np.array([self.kappa_interp(np.log(ai)) for ai in a])
        g = kappa_dot * np.exp(-np.cumsum(kappa_dot * np.gradient(np.log(a))))
        
        source = g * (delta_b + Phi)  # Simplified
        
        return a, source
    
    def compute_cl(self, lmax=100, n_k=50):
        """Compute C_l^TT using LOS integration."""
        print(f"[BoltzmannComplete] Computing C_l to lmax={lmax}...")
        
        # k integration range
        k_min, k_max = 1e-4, 1.0  # h/Mpc
        k_vals = np.logspace(np.log10(k_min), np.log10(k_max), n_k)
        
        # Primordial power
        A_s = 2.1e-9
        n_s = 0.965
        k_pivot = 0.05
        P_R = A_s * (k_vals / k_pivot)**(n_s - 1)
        
        ells = np.arange(2, lmax+1)
        Cl = np.zeros(lmax+1)
        
        for ell in ells[:20]:  # Compute first 20 for speed
            integral_k = 0
            
            for i, k in enumerate(k_vals[:-1]):
                dk = k_vals[i+1] - k
                
                # Get source function
                a, S = self.compute_source_function(k)
                
                # LOS integral with spherical Bessel
                # Delta_l(k) = ∫ dτ S(k,τ) j_l(k(τ_0-τ))
                
                # Approximate: j_l(x) ≈ sin(x - lπ/2)/x for large x
                # Use scipy spherical Bessel
                
                tau = np.cumsum(1 / (a * np.array([self.Hubble(ai) for ai in a])))
                tau_0 = tau[-1]
                
                x = k * (tau_0 - tau)
                jl = spherical_jn(ell, x)
                
                Delta = np.trapz(S * jl, tau)
                
                # Add to integral
                integral_k += k**2 * abs(Delta)**2 * P_R[i] * dk
            
            Cl[ell] = 2/np.pi * integral_k
            Cl[ell] *= (2.725e6)**2  # Convert to μK²
        
        # Fill higher ell with approximate scaling
        for ell in ells[20:]:
            Cl[ell] = Cl[20] * (20/ell)**2
        
        return Cl


if __name__ == "__main__":
    print("="*70)
    print("COMPLETE BOLTZMANN SOLVER")
    print("="*70)
    
    solver = CompleteBoltzmannSolver(H0=70, Omega_b=0.045, Omega_cdm=0.25, lmax=30)
    
    # Test single mode
    k = 0.01
    a, y = solver.solve_mode(k, n_steps=200)
    
    print(f"\nMode k={k} h/Mpc:")
    print(f"  δ_c(a=1) = {y[-1,0]:.4f}")
    print(f"  δ_b(a=1) = {y[-1,1]:.4f}")
    print(f"  Φ(a=1) = {y[-1,8]:.4e}")
    
    print("\n" + "="*70)
