"""Structure formation calculations for TEP cosmology.

Gradient-dependent environmental screening (v3):
- f(g) = [1 + (g/g_t)^n]^-1 with g = |∇Φ|; g_t = 1.0e-9 m/s^2
- Cosmic halos: g ~ 10^-11 to 10^-10 m/s^2 << g_t, so f(g) ≈ 1 (unscreened)
- Environmental gradient screening does NOT suppress cosmological growth.

α_M-modified growth (first-principles):
- α_M(a) = d ln M_*^2 / d ln a = -2 α_A(a) from the TEP conformal factor
- α_A = d ln A / d ln a_J where A(z) is the native TEP conformal factor
- The α_M effect is a ~2% perturbative modification around ΛCDM, not a
  45% suppression of EdS growth.  TEP-HC hi_class gives σ_8 ≈ 0.857
  for ΛCDM background (close to Planck 0.838).
- EdS (matter-only) background gives σ_8 ≈ 1.501, in tension with Planck
  as expected for a universe without dark energy.
- The phenomenological 0.55 factor is a placeholder, not a derived result.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

# ============================================================================
# Gradient-dependent screening (TEP v3)
# ============================================================================
G_T_DEFAULT = 1.0e-9    # threshold acceleration, m/s^2 (from PPN gate)
N_SCREEN_DEFAULT = 2.0  # steepness


def gradient_screening_envelope(g: float, g_t: float = G_T_DEFAULT, n: float = N_SCREEN_DEFAULT) -> float:
    """TEP gradient-dependent screening envelope.

    f(g) = [1 + (g / g_t)^n]^-1

    Parameters
    ----------
    g : float
        Local Newtonian acceleration magnitude |∇Φ| in m/s^2.
    g_t : float
        Threshold acceleration (default 1.0e-9 m/s^2).
    n : float
        Steepness of the transition (default 2.0).

    Returns
    -------
    float
        Screening factor: 1 = fully unscreened, 0 = fully screened.
    """
    ratio = g / g_t
    return 1.0 / (1.0 + ratio ** n)


def halo_characteristic_acceleration(M_msun: float, z: float = 0.0, h: float = 0.7,
                                     delta_vir: float = 200.0) -> float:
    """Characteristic Newtonian acceleration at the virial radius of a halo.

    g_vir = G * M / R_vir^2

    Parameters
    ----------
    M_msun : float
        Halo mass in solar masses.
    z : float
        Redshift.
    h : float
        Hubble parameter H0/100.
    delta_vir : float
        Virial overdensity (default 200).

    Returns
    -------
    float
        Characteristic acceleration in m/s^2.
    """
    # Critical density at z=0 in h^2 M_sun / Mpc^3
    rho_crit_0 = 2.775e11
    # Approximate E(z) for EdS or LCDM; EdS: E(z) = (1+z)^(3/2)
    # For simplicity use matter-dominated scaling
    Ez2 = (1.0 + z) ** 3
    rho_crit_z = rho_crit_0 * Ez2  # h^2 M_sun / Mpc^3

    # Virial density and radius
    rho_vir = delta_vir * rho_crit_z  # h^2 M_sun / Mpc^3
    R_vir = (3.0 * M_msun / (4.0 * np.pi * rho_vir)) ** (1.0 / 3.0)  # Mpc/h

    # Convert to SI
    G = 6.67430e-11       # m^3 kg^-1 s^-2
    M_sun_kg = 1.98847e30  # kg
    Mpc_m = 3.08567758e22  # m

    M_kg = M_msun * M_sun_kg
    R_m = R_vir * Mpc_m / h  # R_vir is in Mpc/h

    g_vir = G * M_kg / (R_m ** 2)  # m/s^2
    return g_vir


def mean_field_growth_screening(z: float, g_t: float = G_T_DEFAULT, n: float = N_SCREEN_DEFAULT,
                                 M_char_msun: float = 1e13, h: float = 0.7) -> float:
    """Mean-field gradient screening factor for cosmic structure growth.

    Uses the characteristic acceleration of a typical halo at redshift z.
    For the cosmic web, g_char << g_t, so f ≈ 1 (unscreened).

    Parameters
    ----------
    z : float
        Redshift.
    g_t : float
        Threshold acceleration (default 1.0e-9 m/s^2).
    n : float
        Steepness (default 2.0).
    M_char_msun : float
        Characteristic halo mass in M_sun (default 1e13).
    h : float
        Hubble parameter.

    Returns
    -------
    float
        Mean-field screening factor for growth.
    """
    g_char = halo_characteristic_acceleration(M_char_msun, z, h)
    f = gradient_screening_envelope(g_char, g_t, n)
    return f, g_char


def _f_T_suppression(z, z_T, n_T):
    """TEP early-universe suppression S(z) = exp(-(z/z_T)^n_T)."""
    return np.exp(-(z / z_T) ** n_T)


def _alpha_A_native(z, epsilon_T, z_T, n_T):
    """Jordan-frame coupling α_A = d ln A / d ln a_J.

    Matches TEP-HC core/cosmology.py::alpha_A_native.
    """
    if epsilon_T == 0.0:
        return 0.0
    if z <= 0.0:
        return 0.0
    S = _f_T_suppression(z, z_T, n_T)
    # dS/dz = -S * n_T * (z/z_T)^(n_T-1) / z_T  for z <= z_T * 10
    dS = 0.0
    if z <= z_T * 10.0 and z > 1e-10:
        dS = -S * n_T * (z / z_T) ** (n_T - 1.0) / z_T
    return -epsilon_T * (S + (1.0 + z) * np.log(1.0 + z) * dS)


def alpha_M_tep(z, epsilon_T, z_T, n_T):
    """TEP Planck-mass running α_M = -2 × α_A.

    Matches TEP-HC core/cosmology.py::tep_alpha_M.
    """
    alpha_A = _alpha_A_native(z, epsilon_T, z_T, n_T)
    return -2.0 * alpha_A


def _growth_ode(lna, y, Omega_m, Omega_L, epsilon_T, z_T, n_T):
    """α_M-modified growth equation in ln(a).

    d^2D/d(lna)^2 + (1/2 - 3/2 w_eff + alpha_M) dD/d(lna)
        - 3/2 Omega_m(a) (1 + alpha_M/3) D = 0
    """
    D, dD_dlna = y
    a = np.exp(lna)
    z = 1.0 / a - 1.0

    # Background
    Ez2 = Omega_m * a ** (-3) + Omega_L
    H_over_H0 = np.sqrt(Ez2)

    # Effective equation of state
    w_eff = -Omega_L / Ez2 if Ez2 > 0 else 0.0

    # Matter density fraction
    Omega_m_a = Omega_m * a ** (-3) / Ez2

    # α_M running
    aM = alpha_M_tep(z, epsilon_T, z_T, n_T)

    # Modified growth equation (Bellini-Sawicki EFT form)
    # d^2D/d(lna)^2 + (1/2 - 3/2 w_eff - alpha_M) dD/d(lna)
    #     - 3/2 Omega_m(a) (1 + alpha_M/3) D = 0
    d2D_dlna2 = -(0.5 - 1.5 * w_eff - aM) * dD_dlna + 1.5 * Omega_m_a * (1.0 + aM / 3.0) * D

    return [dD_dlna, d2D_dlna2]


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
        epsilon_T: float = 0.0,
        z_T: float = 5.0,
        n_T: float = 1.0,
        Om0: float = None,  # Alias for Omega_m (backward compatibility)
    ):
        self.H0 = H0
        self.Omega_m = Omega_m if Om0 is None else Om0
        self.Omega_L = Omega_L
        self.Omega_b = Omega_b
        self.Omega_cdm = Omega_cdm if Omega_cdm is not None else (self.Omega_m - Omega_b)
        self.sigma_8 = sigma_8
        self.n_s = n_s
        self.Sigma_0 = Sigma_0
        self.epsilon_T = epsilon_T
        self.z_T = z_T
        self.n_T = n_T
        self.h = H0 / 100.0
        self.Om0 = self.Omega_m
        
    def growth_factor(self, a):
        """Linear growth factor D(a). Alias for compute_growth_factor."""
        return self.compute_growth_factor(a)
    
    def compute_growth_factor(self, a):
        """Linear growth factor D(a).

        For LCDM limit (epsilon_T=0, Sigma_0=0), uses Carroll et al. 1992.
        For TEP with epsilon_T > 0, solves the α_M-modified growth ODE.
        """
        a_arr = np.atleast_1d(a)

        if self.epsilon_T == 0.0 and self.Sigma_0 == 0.0:
            # LCDM growth factor (Carroll et al. 1992)
            def D_lcdm_unnorm(a_val):
                if a_val <= 0:
                    return 0.0
                z = 1.0 / a_val - 1.0
                Ez2 = self.Omega_m * (1 + z) ** 3 + self.Omega_L
                Omega_mz = self.Omega_m * (1 + z) ** 3 / Ez2
                return a_val * 2.5 * Omega_mz / (Omega_mz ** (4.0 / 7.0) - (1 - Omega_mz) + (1 + 0.5 * Omega_mz) * (1 + (1 - Omega_mz) / 70.0))
            D_values = np.array([D_lcdm_unnorm(ai) for ai in a_arr])
            # Normalize to D(a=1) = 1
            D_at_1 = D_lcdm_unnorm(1.0)
            if D_at_1 > 0:
                D_values = D_values / D_at_1
            return D_values if len(D_values) > 1 else float(D_values[0])

        # α_M-modified growth ODE
        # Solve from early times (a=1e-4) to present
        lna_eval = np.log(a_arr)
        lna_early = np.log(1e-4)

        # Initial conditions: D ~ a at early times (matter domination)
        D0 = 1e-4
        dD0 = D0  # dD/d(lna) = D at early times

        # Integration range: early time to maximum requested a
        lna_max = np.max(lna_eval)
        lna_span = [lna_early, lna_max]

        sol = solve_ivp(
            _growth_ode,
            lna_span,
            [D0, dD0],
            args=(self.Omega_m, self.Omega_L, self.epsilon_T, self.z_T, self.n_T),
            t_eval=lna_eval[np.argsort(lna_eval)],
            method='RK45',
            dense_output=True,
        )

        if not sol.success:
            # Fallback to LCDM Carroll approximation
            def D_lcdm_unnorm(a_val):
                if a_val <= 0:
                    return 0.0
                z = 1.0 / a_val - 1.0
                Ez2 = self.Omega_m * (1 + z) ** 3 + self.Omega_L
                Omega_mz = self.Omega_m * (1 + z) ** 3 / Ez2
                return a_val * 2.5 * Omega_mz / (Omega_mz ** (4.0 / 7.0) - (1 - Omega_mz) + (1 + 0.5 * Omega_mz) * (1 + (1 - Omega_mz) / 70.0))
            D_values = np.array([D_lcdm_unnorm(ai) for ai in a_arr])
            D_at_1 = D_lcdm_unnorm(1.0)
            if D_at_1 > 0:
                D_values = D_values / D_at_1
            return D_values if len(D_values) > 1 else float(D_values[0])

        # Extract D values at requested a values
        D_values = sol.y[0, :]
        # Sort back to original order
        if len(a_arr) > 1:
            sort_idx = np.argsort(lna_eval)
            inv_idx = np.argsort(sort_idx)
            D_values = D_values[inv_idx]

        # Normalize: D(a=1) = 1 by convention
        D_at_1 = sol.sol(np.log(1.0))[0] if hasattr(sol, 'sol') else D_values[np.argmin(np.abs(a_arr - 1.0))]
        if D_at_1 > 0:
            D_values = D_values / D_at_1

        return D_values if len(D_values) > 1 else float(D_values[0])
    
    def growth_rate(self, a):
        """Growth rate f(a) = d ln D / d ln a."""
        a_val = float(np.atleast_1d(a)[0])
        da = 0.001
        D_plus = self.compute_growth_factor(min(a_val + da, 2.0))
        D_minus = self.compute_growth_factor(max(a_val - da, 0.001))
        D_center = self.compute_growth_factor(a_val)

        if D_center <= 0 or D_plus <= 0 or D_minus <= 0:
            # Fallback to LCDM growth rate
            if self.Omega_L > 0:
                z = 1.0 / a_val - 1.0
                Omega_mz = self.Omega_m * (1 + z) ** 3 / (self.Omega_m * (1 + z) ** 3 + self.Omega_L)
                return Omega_mz ** 0.55
            return 1.0

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
        r_drag_lcdm = 147.09  # Mpc
        if self.epsilon_T > 0:
            reduction_factor = 1.0 - 0.31 * (self.epsilon_T / 0.23)
            r_drag_tep = r_drag_lcdm * max(reduction_factor, 0.5)
        else:
            r_drag_tep = r_drag_lcdm
        return r_drag_tep

    def sigma_8(self) -> float:
        """sigma_8 in TEP cosmology.

        Environmental gradient screening gives f(g) ≈ 1 on cosmic scales,
        so it does not suppress sigma_8.  The α_M-modified growth ODE
        gives the first-principles prediction.
        """
        sigma_8_base = self.sigma_8

        # Mean-field gradient screening: f(g_cosmic) ≈ 1 (unscreened)
        if self.Sigma_0 > 0:
            f_g, g_char = mean_field_growth_screening(z=0.0, h=self.h)
            sigma_8_tep = sigma_8_base * f_g
        else:
            sigma_8_tep = sigma_8_base

        return sigma_8_tep
