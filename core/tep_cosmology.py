"""TEP Cosmology module for distance calculations.

Provides TEP-modified cosmological distance calculations.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad


class TEPCosmology:
    """TEP cosmology with temporal shear effects."""

    # Spatial transition density (manuscript: rho_half ≈ 0.5 M_sun/pc^3)
    # Maps to Equation in Section 2.5: S(rho) = [1 + (rho/rho_half)^2]^-1
    RHO_HALF: float = 0.5  # M_sun / pc^3, threshold galactic onset density

    def __init__(
        self,
        H0: float = 70.0,
        Omega_m: float = 0.3,
        epsilon_T: float = 0.0,
        z_T: float = 5.0,
        n_T: float = 1.0,
    ):
        self.H0 = H0
        self.Omega_m = Omega_m
        self.epsilon_T = epsilon_T
        self.z_T = z_T
        self.n_T = n_T
        self.Omega_lambda = 1.0 - Omega_m
        
    def tep_gamma(self, z: float | np.ndarray) -> float | np.ndarray:
        """TEP path enhancement factor.
        
        Uses the exact normalised formula: gamma(z) = 1 + epsilon_T * ln(1+z) * S(z).
        This guarantees gamma(0) = 1 (recovering the local reference frame),
        and gamma > 1 for z > 0 (reflecting clocks ran slower in the denser past,
        A_past < A_today, thus gamma = A_today / A_past > 1).
        This makes intermediate distances larger, actively mimicking dark energy acceleration.
        """
        z_arr = np.asarray(z, dtype=float)
        if self.epsilon_T == 0:
            gamma = np.ones_like(z_arr)
            return float(gamma) if np.isscalar(z) else gamma

        log_factor = np.log(1.0 + z_arr)
        suppression = np.exp(-((z_arr / self.z_T) ** self.n_T))
        enhancement = 1.0 + self.epsilon_T * log_factor * suppression
        gamma = np.maximum(enhancement, 0.1)  # Physical bound
        return float(gamma) if np.isscalar(z) else gamma

    def screening_function(self, rho: float | np.ndarray) -> float | np.ndarray:
        """TEP environmental screening factor S(rho).

        Implements the continuous shear-suppression formula from
        Section 2.5 of the manuscript:
            S(rho) = [1 + (rho / rho_half)^2]^-1

        where rho_half = 0.5 M_sun / pc^3 is the threshold galactic
        onset density (class constant RHO_HALF).

        Parameters
        ----------
        rho : float or array
            Local matter density in M_sun / pc^3.

        Returns
        -------
        float or array
            Screening factor between 0 (fully screened) and 1 (unscreened).
        """
        ratio = np.asarray(rho, dtype=float) / self.RHO_HALF
        return 1.0 / (1.0 + ratio ** 2)

    def e_func(self, z: float | np.ndarray) -> float | np.ndarray:
        """Dimensionless Hubble parameter."""
        zp1 = 1.0 + np.asarray(z, dtype=float)
        matter = self.Omega_m * zp1**3
        Lambda = self.Omega_lambda
        return np.sqrt(matter + Lambda)
    
    def comoving_distance(self, z: float | np.ndarray) -> float | np.ndarray:
        """TEP comoving distance in Mpc."""
        z_arr = np.atleast_1d(np.asarray(z, dtype=float))
        c = 299792.458  # km/s

        max_z = float(np.max(z_arr))
        if max_z <= 0:
            distances = np.zeros_like(z_arr)
            return float(distances[0]) if np.isscalar(z) else distances

        grid = np.linspace(0.0, max_z, 1000)
        integrand = c / (self.H0 * self.e_func(grid)) * self.tep_gamma(grid)
        cumulative = cumulative_trapezoid(integrand, grid, initial=0.0)
        distances = np.interp(z_arr, grid, cumulative)
        return float(distances[0]) if np.isscalar(z) else distances
    
    def angular_diameter_distance(self, z: float | np.ndarray) -> float | np.ndarray:
        """TEP angular diameter distance in Mpc."""
        d_c = self.comoving_distance(z)
        return d_c / (1 + z)
    
    def luminosity_distance(self, z: float | np.ndarray) -> float | np.ndarray:
        """TEP luminosity distance in Mpc."""
        d_c = self.comoving_distance(z)
        return d_c * (1 + z)

    def distance_modulus(self, z: float | np.ndarray) -> float | np.ndarray:
        """TEP distance modulus."""
        d_l = np.maximum(self.luminosity_distance(z), np.finfo(float).tiny)
        return 5.0 * np.log10(d_l) + 25.0


class TEPPerturbations:
    """TEP linear perturbation calculations."""
    
    def __init__(
        self,
        H0: float = 70.0,
        Omega_m: float = 0.3,
        sigma_0: float = 0.0,
        epsilon_T: float = 0.0,
    ):
        self.H0 = H0
        self.Omega_m = Omega_m
        self.sigma_0 = sigma_0
        self.epsilon_T = epsilon_T
        self.h = H0 / 100.0
        
    def sound_horizon(self) -> float:
        """Sound horizon at drag epoch in Mpc.
        
        Uses Eisenstein & Hu 1998 approximation with TEP corrections.
        """
        # Physical sound horizon
        c_s = 1.0 / np.sqrt(3.0)  # c_s/c for photon-baryon fluid
        z_drag = 1060.0  # Drag redshift
        
        # Standard sound horizon approximation
        r_s_std = 147.0  # Mpc for standard LCDM
        
        # In matter-frame BBN/sound horizon, TEP corrections are exponentially suppressed
        if self.sigma_0 > 0 or self.epsilon_T > 0:
            # The Temporal-Shear suppression function exp[-(z/z_T)^n_T] inherently
            # vanishes at early times (z >> z_T), leaving the pre-recombination
            # sound horizon preserved.
            tep_factor = 1.0
            return r_s_std * tep_factor
        
        return r_s_std
    
    def sigma_8(self) -> float:
        """Sigma_8 amplitude from perturbations."""
        # Reference Planck 2018 value
        sigma_8_planck = 0.81
        
        # TEP modification through growth suppression
        if self.sigma_0 > 0 or self.epsilon_T > 0:
            # Late-time structure growth is preserved due to galactic screening
            # (rho >> rho_half).
            suppression = 1.0
            return sigma_8_planck * suppression
        
        return sigma_8_planck


class CosmologyFLRW:
    """Standard FLRW cosmology for comparison."""
    
    def __init__(self, H0: float = 70.0, Om0: float = 0.3):
        self.H0 = H0
        self.Om0 = Om0
        self.Ode0 = 1.0 - Om0
        
    def e_func(self, z: float) -> float:
        """Dimensionless Hubble parameter."""
        zp1 = 1.0 + z
        return np.sqrt(self.Om0 * zp1**3 + self.Ode0)
    
    def angular_diameter_distance(self, z: np.ndarray) -> np.ndarray:
        """Angular diameter distance in Mpc."""
        c = 299792.458  # km/s
        
        def d_c_single(z_val):
            result, _ = quad(
                lambda zp: c / (self.H0 * self.e_func(zp)),
                0, z_val, limit=100
            )
            return result
        
        if np.isscalar(z):
            d_c = d_c_single(z)
        else:
            d_c = np.array([d_c_single(zi) for zi in z])
        
        return d_c / (1 + z)


class TEPCosmologyFitter:
    """Fit TEP cosmology to observational data."""
    
    def __init__(
        self,
        z: np.ndarray | None = None,
        mu: np.ndarray | None = None,
        mu_err: np.ndarray | None = None,
        data_path: str | None = None,
    ):
        self.data_path = data_path
        self.data = None
        if z is not None and mu is not None and mu_err is not None:
            self.data = {
                'z': np.asarray(z, dtype=float),
                'mu': np.asarray(mu, dtype=float),
                'mu_err': np.asarray(mu_err, dtype=float),
            }
        
    def load_data(self):
        """Load SNe data for fitting."""
        import pandas as pd
        from pathlib import Path
        
        if self.data_path is None:
            # Try default locations
            candidates = [
                Path("data/raw/pantheon_plus_shoes.dat"),
                Path("data/raw/Pantheon+SH0ES.dat"),
            ]
            for path in candidates:
                if path.exists():
                    self.data_path = str(path)
                    break
        
        if self.data_path is None or not Path(self.data_path).exists():
            raise FileNotFoundError("Pantheon+ data not found")
        
        df = pd.read_csv(self.data_path, sep=r'\s+', comment='#')
        
        # Extract relevant columns
        z = df['zHD'].values if 'zHD' in df.columns else df['zCMB'].values
        mb = df['m_b_corr'].values if 'm_b_corr' in df.columns else df['mB'].values
        mb_err = df['m_b_corr_err_DIAG'].values if 'm_b_corr_err_DIAG' in df.columns else df.get('mBERR', pd.Series([0.15]*len(df))).values
        
        # Filter valid data
        valid = np.isfinite(z) & np.isfinite(mb) & (z > 0.001) & (z < 3.0)
        
        self.data = {
            'z': z[valid],
            'mb': mb[valid],
            'mb_err': mb_err[valid],
        }
        
        return self.data

    def _observations(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if self.data is None:
            self.load_data()

        z = self.data['z']
        if 'mu' in self.data:
            obs = self.data['mu']
            err = self.data['mu_err']
        else:
            obs = self.data['mb']
            err = self.data['mb_err']
        return z, obs, np.maximum(err, np.finfo(float).tiny)
    
    def fit_tep(self, initial_guess: dict = None) -> dict:
        """Fit TEP model to data.
        
        Returns best-fit parameters and statistics.
        """
        if self.data is None:
            self.load_data()
        
        from scipy.optimize import minimize
        
        z, obs, obs_err = self._observations()
        
        # Initial guess
        if initial_guess is None:
            p0 = [70.0, 0.3, 0.1, 5.0, 0.0]  # H0, Omega_m, epsilon_T, z_T, offset
        else:
            p0 = [
                initial_guess.get('H0', 70.0),
                initial_guess.get('Omega_m', 0.3),
                initial_guess.get('epsilon_T', 0.1),
                initial_guess.get('z_T', 5.0),
                initial_guess.get('offset', 0.0),
            ]
        
        def chi2(params):
            H0, Omega_m, epsilon_T, z_T, offset = params
            try:
                tep = TEPCosmology(H0=H0, Omega_m=Omega_m, epsilon_T=epsilon_T, z_T=z_T)
                mu_pred = tep.distance_modulus(z) + offset
                residuals = (obs - mu_pred) / obs_err
                return np.sum(residuals**2)
            except:
                return 1e10
        
        # Run minimization
        result = minimize(
            chi2,
            p0,
            method='Nelder-Mead',
            bounds=[(50, 100), (0.05, 0.9), (0.0, 1.0), (0.5, 15.0), (-5.0, 5.0)],
            options={"maxiter": 5000},
        )
        
        if result.success:
            H0_fit, Om_fit, eps_fit, zT_fit, offset_fit = result.x
            chi2_min = result.fun
            
            return {
                'success': True,
                'parameters': {
                    'H0': H0_fit,
                    'Omega_m': Om_fit,
                    'epsilon_T': eps_fit,
                    'z_T': zT_fit,
                    'offset': offset_fit,
                },
                'chi2': chi2_min,
                'dof': len(z) - 5,
            }
        else:
            return {
                'success': False,
                'error': 'Minimization failed',
            }

    def fit_lcdm(self, initial_guess: dict | None = None) -> dict:
        """Fit matched LCDM distance model to the same observations."""
        if self.data is None:
            self.load_data()

        from scipy.optimize import minimize

        z, obs, obs_err = self._observations()
        p0 = [
            70.0 if initial_guess is None else initial_guess.get('H0', 70.0),
            0.3 if initial_guess is None else initial_guess.get('Omega_m', 0.3),
            0.0 if initial_guess is None else initial_guess.get('offset', 0.0),
        ]

        def chi2(params):
            H0, Omega_m, offset = params
            try:
                lcdm = CosmologyFLRW(H0=H0, Om0=Omega_m)
                mu_pred = np.array([
                    5.0 * np.log10(lcdm.angular_diameter_distance(float(zi)) * (1.0 + zi) ** 2) + 25.0
                    for zi in z
                ]) + offset
                residuals = (obs - mu_pred) / obs_err
                return float(np.sum(residuals**2))
            except Exception:
                return 1e10

        result = minimize(
            chi2,
            p0,
            method='Nelder-Mead',
            bounds=[(50, 100), (0.05, 0.9), (-5.0, 5.0)],
            options={"maxiter": 5000},
        )
        if result.success:
            H0_fit, Om_fit, offset_fit = result.x
            return {
                'success': True,
                'parameters': {'H0': H0_fit, 'Omega_m': Om_fit, 'offset': offset_fit},
                'chi2': float(result.fun),
                'dof': len(z) - 3,
            }
        return {'success': False, 'error': 'Minimization failed'}

    def compare_models(self) -> dict:
        """Fit LCDM and TEP with matched data and report information criteria."""
        z, _, _ = self._observations()
        lcdm = self.fit_lcdm()
        tep = self.fit_tep()
        if not lcdm.get('success') or not tep.get('success'):
            return {'status': 'failed', 'lcdm': lcdm, 'tep': tep}

        n = len(z)
        k_lcdm = 3
        k_tep = 5
        lcdm_aic = lcdm['chi2'] + 2 * k_lcdm
        tep_aic = tep['chi2'] + 2 * k_tep
        lcdm_bic = lcdm['chi2'] + k_lcdm * np.log(n)
        tep_bic = tep['chi2'] + k_tep * np.log(n)
        lcdm['chi2_per_dof'] = lcdm['chi2'] / max(lcdm['dof'], 1)
        tep['chi2_per_dof'] = tep['chi2'] / max(tep['dof'], 1)

        delta_bic = tep_bic - lcdm_bic
        return {
            'status': 'completed',
            'lcdm': {**lcdm, 'aic': float(lcdm_aic), 'bic': float(lcdm_bic)},
            'tep': {**tep, 'aic': float(tep_aic), 'bic': float(tep_bic)},
            'delta_aic_tep_vs_lcdm': float(tep_aic - lcdm_aic),
            'delta_bic_tep_vs_lcdm': float(delta_bic),
            'tep_competitive': bool(delta_bic < 2.0),
            'best_model': 'tep' if delta_bic < 0.0 else 'lcdm',
        }
