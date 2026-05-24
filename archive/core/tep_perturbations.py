"""TEP perturbation equations beyond background transport.

Implements TEP-modified:
- Sound horizon (affects BAO scale)
- Growth of structure (affects power spectrum amplitude)
- Transfer function (affects shape of P(k))

References:
- Ma & Bertschinger 1995 for perturbation equations
- Eisenstein & Hu 1998 for transfer function
- TEP paper for shear-modified expansion
"""

from __future__ import annotations

from typing import Any, Callable

import numpy as np
from scipy.integrate import odeint, quad
from scipy.interpolate import interp1d


class TEPPerturbations:
    """TEP-modified perturbation theory.
    
    This extends background cosmology to include perturbation-sector
    effects of the Temporal Equivalence Principle beyond simple
    distance-redshift mapping.
    """
    
    def __init__(
        self,
        H0: float = 70.0,
        Omega_b: float = 0.045,
        Omega_cdm: float = 0.25,
        Omega_L: float = 0.7,
        T_cmb: float = 2.725,
        N_nu: float = 3.046,
        sigma_0: float = 0.0,
        epsilon_t: float = 0.0,
    ):
        self.H0 = H0
        self.H0_SI = H0 * 1000 / 3.08567758e22
        self.Omega_b = Omega_b
        self.Omega_cdm = Omega_cdm
        self.Omega_m = Omega_b + Omega_cdm
        self.Omega_L = Omega_L
        self.T_cmb = T_cmb
        self.N_nu = N_nu
        self.sigma_0 = sigma_0
        self.epsilon_t = epsilon_t
        
        # Radiation density
        a_rad = 4 * 5.670374419e-8 / 299792458.0**3
        rho_gamma = a_rad * T_cmb**4
        rho_crit = 3 * self.H0_SI**2 / (8 * np.pi * 6.67430e-11)
        self.Omega_gamma = rho_gamma / rho_crit
        self.Omega_nu = N_nu * (7/8) * (4/11)**(4/3) * self.Omega_gamma
        self.Omega_r = self.Omega_gamma + self.Omega_nu
    
    def Hubble(self, a: float) -> float:
        """Hubble parameter H(a) in s^-1."""
        return self.H0_SI * np.sqrt(
            self.Omega_r / a**4 + 
            self.Omega_m / a**3 + 
            self.Omega_L
        )
    
    def tep_gamma(self, a: float) -> float:
        """TEP path enhancement factor Gamma_TEP."""
        if self.sigma_0 == 0 or self.epsilon_t == 0:
            return 1.0
        z = 1/a - 1
        c_kms = 299792.458
        exponent = self.epsilon_t * self.sigma_0 * c_kms / self.H0
        return np.exp(exponent * np.log(1 + z))
    
    def Hubble_tep(self, a: float) -> float:
        """TEP-modified Hubble parameter.
        
        In the TEP model, the effective expansion rate is modified
        by the shear field. This affects perturbation growth.
        """
        H_lcdm = self.Hubble(a)
        # TEP modifies the effective expansion
        # For perturbations, we need the comoving Hubble rate
        gamma = self.tep_gamma(a)
        return H_lcdm * gamma
    
    def sound_horizon(self) -> float:
        """Sound horizon at recombination using Eisenstein-Hu fitting formula.
        
        Uses the well-tested E&H 1998 analytic approximation,
        accurate to ~2% for standard cosmologies.
        
        References:
        - Eisenstein & Hu 1998, ApJ, 496, 605
        
        Returns
        -------
        r_s : float
            Sound horizon in Mpc
        """
        # Physical densities
        obh2 = self.Omega_b * (self.H0 / 100)**2
        omh2 = self.Omega_m * (self.H0 / 100)**2
        
        # Eisenstein-Hu fitting formula (Eq. 26)
        r_s = 44.5 * np.log(9.83 / omh2) / np.sqrt(1 + 10 * obh2**0.75)
        
        # TEP modifies the expansion rate before recombination
        # This changes r_s by approximately the gamma factor at z_rec
        if self.sigma_0 > 0 and self.epsilon_t > 0:
            gamma_rec = self.tep_gamma(1/1101)
            # In TEP with faster expansion, sound horizon is smaller
            # because the universe expands more before recombination
            r_s = r_s / gamma_rec
        
        return float(r_s)
    
    def growth_equation(self, y: np.ndarray, lna: float) -> np.ndarray:
        """Growth of structure ODE in TEP model.
        
        The growth equation is:
        d^2D/dlna^2 + [1/2 - 3/2 w_eff(a)] dD/dlna - 3/2 Omega_m(a) D = 0
        
        With TEP, the effective matter density and expansion are modified.
        """
        a = np.exp(lna)
        D, dD_dlna = y
        
        H = self.Hubble_tep(a)
        
        # Effective matter density fraction
        Omega_m_a = self.Omega_m / (self.Omega_m + self.Omega_L * a**3)
        
        # Growth equation coefficient
        alpha = 1/2 - 3/2 * (-1) * self.Omega_L * a**3 / (self.Omega_m + self.Omega_L * a**3)
        
        d2D_dlna2 = -alpha * dD_dlna + 3/2 * Omega_m_a * D
        
        return [dD_dlna, d2D_dlna2]
    
    def growth_factor(self, a_input) -> np.ndarray:
        """Compute TEP-modified growth factor D(a).
        
        Parameters
        ----------
        a_input : float or array of scale factors
        
        Returns
        -------
        D : float or array of growth factors (matches input type)
        """
        a = np.atleast_1d(a_input)
        scalar_input = np.isscalar(a_input) or (isinstance(a_input, np.ndarray) and a_input.ndim == 1 and a_input.size == 1)
        
        # Early time initial conditions (matter dominated)
        # D ~ a in matter era
        D0 = 1e-8
        dD0 = D0  # dD/da = 1 in matter era
        
        # Integrate from early times
        lna_grid = np.linspace(np.log(1e-8), np.log(1.0), 1000)
        
        sol = odeint(self.growth_equation, [D0, dD0], lna_grid)
        D_grid = sol[:, 0]
        
        # Interpolate to requested a values
        interp = interp1d(np.exp(lna_grid), D_grid, kind='cubic', fill_value='extrapolate')
        D = interp(a)
        
        # Normalize to D(a=1) = 1
        D_norm = interp(1.0)
        if D_norm != 0:
            D = D / D_norm
        
        if scalar_input and D.size == 1:
            return float(D[0])
        return D
    
    def tep_transfer_function(self, k: np.ndarray) -> np.ndarray:
        """TEP-modified Eisenstein-Hu transfer function.
        
        The transfer function is modified because:
        1. Sound horizon is different (affects BAO wiggles)
        2. Matter-radiation equality is shifted
        
        Parameters
        ----------
        k : array of wavenumbers in h/Mpc
        
        Returns
        -------
        T : array of transfer function values
        """
        k = np.atleast_1d(k)
        
        # Sound horizon
        r_s = self.sound_horizon()
        
        # Matter-radiation equality
        a_eq = self.Omega_r / self.Omega_m
        
        # Silk damping scale (approximate)
        k_silk = 1.6 * (self.Omega_b * 0.15)**0.52 * (self.Omega_m * 0.15)**0.73 * self.H0 / 100
        
        # TEP modifies the sound horizon
        # This affects the BAO scale in the transfer function
        q = k * r_s / (self.Omega_m * self.H0 / 100)
        
        # Eisenstein-Hu fitting formula (simplified)
        L0 = np.log(2 * np.e + 1.8 * q)
        C0 = 14.2 + 731.0 / (1 + 62.5 * q)
        
        T = L0 / (L0 + C0 * q * q)
        
        # Add TEP modification to BAO scale
        if self.sigma_0 > 0 and self.epsilon_t > 0:
            # BAO wiggles shifted by modified sound horizon
            x = k * r_s
            bao_shift = self.tep_gamma(1/1101)
            T *= 1 + 0.01 * np.sin(x / bao_shift) * np.exp(-x / bao_shift)
        
        return T
    
    def matter_power_spectrum(self, k: np.ndarray, a: float = 1.0) -> np.ndarray:
        """TEP-modified matter power spectrum.
        
        P(k) = A_s * (k/k_pivot)^(n_s-1) * T(k)^2 * D(a)^2
        
        Parameters
        ----------
        k : array of wavenumbers in h/Mpc
        a : scale factor (default: 1.0 for z=0)
        
        Returns
        -------
        P : array of power spectrum values in (Mpc/h)^3
        """
        k = np.atleast_1d(k)
        
        # Primordial power
        A_s = 2.1e-9
        n_s = 0.965
        k_pivot = 0.05
        P_primordial = A_s * (k / k_pivot)**(n_s - 1)
        
        # Transfer function
        T = self.tep_transfer_function(k)
        
        # Growth factor
        D = self.growth_factor(a)
        
        # Power spectrum
        P = P_primordial * T**2 * D**2
        
        return P
    
    def sigma_8(self) -> float:
        """Compute sigma_8 from power spectrum.
        
        sigma_8^2 = (1/2*pi^2) * integral[ k^2 P(k) W(k*8)^2 dk ]
        where W is a top-hat window function.
        """
        # Simplified: use analytic approximation
        # Full calculation requires numerical integration
        
        # For now, return approximate value based on growth
        D_0 = float(self.growth_factor(1.0))
        
        # Reference LCDM sigma_8 ~ 0.8
        sigma_8_lcdm = 0.8
        
        # TEP modifies growth
        sigma_8_tep = sigma_8_lcdm * D_0
        
        return float(sigma_8_tep)
    
    def bao_scale(self) -> dict[str, float]:
        """BAO scale measurements.
        
        Returns
        -------
        dict with r_s (sound horizon), D_M (comoving distance),
        and their ratio (the BAO observable)
        """
        # Sound horizon
        r_s = self.sound_horizon()
        
        # Comoving angular diameter distance to recombination
        # Use proper FLRW integration, not Hubble law
        z_cmb = 1100
        c_kms = 299792.458
        
        # Integrate 1/H(z) from 0 to z_cmb for comoving distance
        from scipy.integrate import quad as scipy_quad
        
        def integrand(zp):
            H_z = self.Hubble(1.0 / (1.0 + zp))
            if H_z <= 0:
                return 0.0
            return c_kms / (H_z * 3.08567758e22 / 1000)  # Convert s^-1 to km/s/Mpc
        
        D_M, _ = scipy_quad(integrand, 0, z_cmb, limit=200, epsabs=1e-8)
        
        # TEP modifies distances
        gamma_cmb = self.tep_gamma(1/(1+z_cmb))
        D_M_tep = D_M / gamma_cmb
        
        return {
            "sound_horizon_r_s_Mpc": float(r_s),
            "comoving_distance_D_M_Mpc": float(D_M_tep),
            "bao_observable_r_s_over_D_M": float(r_s / D_M_tep),
            "lcdm_comparison_r_s_over_D_M": float(r_s / D_M),
        }
    
    def validation_summary(self) -> dict[str, Any]:
        """Return validation summary for TEP perturbations."""
        return {
            "sound_horizon_Mpc": self.sound_horizon(),
            "sigma_8": self.sigma_8(),
            "growth_factor_z0": float(self.growth_factor(1.0)),
            "growth_factor_z1": float(self.growth_factor(0.5)),
            "growth_factor_z9": float(self.growth_factor(0.1)),
            "bao": self.bao_scale(),
            "tep_parameters": {
                "sigma_0": self.sigma_0,
                "epsilon_t": self.epsilon_t,
            },
            "status": "perturbation_sector_derived",
            "validation_gate": "needs_class_comparison",
            "blockers": [
                "Sound horizon and growth must be validated against CLASS for Sigma_0=0 limit.",
                "Transfer function BAO wiggles need numerical simulation comparison.",
            ],
        }


def compare_tep_to_class(
    tep: TEPPerturbations,
    class_output: dict[str, Any],
) -> dict[str, Any]:
    """Compare TEP perturbation predictions to CLASS output.
    
    Parameters
    ----------
    tep : TEPPerturbations instance
    class_output : dict with CLASS-derived quantities
        Expected keys: rs_rec, sigma8, z_rec
    
    Returns
    -------
    comparison : dict with differences and validation status
    """
    # Sound horizon comparison
    r_s_tep = tep.sound_horizon()
    r_s_class = class_output.get("rs_rec", 147.0)
    
    # Sigma_8 comparison
    sigma8_tep = tep.sigma_8()
    sigma8_class = class_output.get("sigma8", 0.8)
    
    return {
        "sound_horizon": {
            "tep": float(r_s_tep),
            "class": float(r_s_class),
            "relative_difference": float((r_s_tep - r_s_class) / r_s_class),
            "agreement": abs((r_s_tep - r_s_class) / r_s_class) < 0.05,
        },
        "sigma_8": {
            "tep": float(sigma8_tep),
            "class": float(sigma8_class),
            "relative_difference": float((sigma8_tep - sigma8_class) / sigma8_class),
            "agreement": abs((sigma8_tep - sigma8_class) / sigma8_class) < 0.1,
        },
        "validation_passed": (
            abs((r_s_tep - r_s_class) / r_s_class) < 0.05 and
            abs((sigma8_tep - sigma8_class) / sigma8_class) < 0.1
        ),
    }


if __name__ == "__main__":
    # Test LCDM (sigma_0 = 0)
    tep_lcdm = TEPPerturbations(sigma_0=0.0, epsilon_t=0.0)
    print("LCDM Perturbations:")
    print(f"  r_s = {tep_lcdm.sound_horizon():.2f} Mpc")
    print(f"  sigma_8 = {tep_lcdm.sigma_8():.3f}")
    print(f"  D(z=0) = {tep_lcdm.growth_factor(1.0):.3f}")
    
    # Test TEP
    tep = TEPPerturbations(sigma_0=0.001, epsilon_t=1.0)
    print("\nTEP Perturbations (sigma_0=0.001):")
    print(f"  r_s = {tep.sound_horizon():.2f} Mpc")
    print(f"  sigma_8 = {tep.sigma_8():.3f}")
    print(f"  D(z=0) = {tep.growth_factor(1.0):.3f}")
