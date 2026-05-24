"""Full TEP Cosmology implementation for competitive model comparison.

This implements the Level 3 TEP-COSMO framework where cosmic expansion is replaced
by temporal transport along photon paths.

Response Coefficient Framework (aligned with TEP manuscripts):
- KAPPA_GAL: Canonical clock-sector coupling (9.6e5 mag from Kingston upon Hull)
- KAPPA_MSP: Pulsar response coefficient (~10^5-10^6 from Caracas)
- ALPHA_SAT: Wide-binary saturation amplitude (0.366 from Kilifi)
- epsilon_T: Phenomenological temporal shear strength parameter

The response coefficients connect different observational probes and ensure
cross-probe consistency as demonstrated across the TEP manuscript series.
"""

import numpy as np
from typing import Optional, Tuple, Dict
from scipy import integrate


# Canonical response coefficients from TEP manuscripts
KAPPA_GAL = 9.6e5  # mag (Kingston upon Hull: Cepheid clock bias, 11-TEP-H0-v0.6)
KAPPA_GAL_UNCERTAINTY = 4.0e5  # mag (uncertainty from Kingston upon Hull fit)
KAPPA_MSP = 1.0e6  # Pulsar response coefficient (Caracas, order-of-magnitude estimate ~10^5-10^6, 10-TEP-COS-v0.6)
ALPHA_SAT = 0.366  # Wide-binary saturation amplitude (Kilifi, 13-TEP-WB-v0.3)


class TEPCosmology:
    """Full TEP cosmological model with temporal transport.
    
    The key equation is:
        ln(1 + z_obs) = ∫_γ H_T dℓ
    
    where H_T is the temporal transport density (not expansion rate).
    
    For a homogeneous but non-integrable temporal field:
        H_T(z) = H_0 * f_T(z)
    
    where f_T(z) captures the temporal shear effects.
    """
    
    def __init__(
        self,
        H0: float = 70.0,
        Omega_m: float = 0.3,
        Omega_L: float = 0.7,
        epsilon_T: float = 0.0,  # TEP amplitude (coupling strength)
        z_T: float = 1.0,        # Characteristic redshift for TEP effects
        n: float = 1.0,          # Redshift dependence exponent
        kappa_gal: Optional[float] = None,  # Cepheid response coefficient (mag)
        kappa_msp: Optional[float] = None,  # Pulsar response coefficient
        alpha_sat: Optional[float] = None,  # Wide-binary saturation amplitude
    ):
        """Initialize TEP cosmology.
        
        Args:
            H0: Hubble parameter at z=0 (km/s/Mpc)
            Omega_m: Matter density parameter
            Omega_L: Dark energy density parameter (can be reinterpreted in TEP)
            epsilon_T: TEP coupling amplitude (not a fraction!)
            z_T: Characteristic redshift scale for TEP transition
            n: Power-law index for redshift dependence
            kappa_gal: Cepheid response coefficient (mag). If None, uses canonical KAPPA_GAL
            kappa_msp: Pulsar response coefficient. If None, uses canonical KAPPA_MSP
            alpha_sat: Wide-binary saturation amplitude. If None, uses canonical ALPHA_SAT
        """
        self.H0 = H0
        self.Omega_m = Omega_m
        self.Omega_L = Omega_L
        self.epsilon_T = epsilon_T
        self.z_T = z_T
        self.n = n
        self.c = 299792.458  # km/s
        
        # Response coefficients (from TEP manuscripts)
        # See: 11-TEP-H0-v0.6 (KAPPA_GAL), 10-TEP-COS-v0.6 (KAPPA_MSP), 13-TEP-WB-v0.3 (ALPHA_SAT)
        self.kappa_gal = kappa_gal if kappa_gal is not None else KAPPA_GAL
        self.kappa_msp = kappa_msp if kappa_msp is not None else KAPPA_MSP
        self.alpha_sat = alpha_sat if alpha_sat is not None else ALPHA_SAT
        
    def _H_T(self, z: float) -> float:
        """Temporal transport density at redshift z.
        
        This is the TEP replacement for the Hubble parameter.
        For epsilon_T = 0, reduces to standard H(z).
        """
        # Standard H(z) for reference with radiation component
        zp1 = 1 + z
        Omega_r = 9.0e-5  # Radiation density parameter (photons + neutrinos)
        Hz_standard = self.H0 * np.sqrt(
            self.Omega_m * zp1**3 + 
            Omega_r * zp1**4 + 
            self.Omega_L
        )
        
        # TEP modification: additional temporal transport term
        # This accumulates along the photon path
        if self.epsilon_T > 0:
            # TEP correction factor
            f_T_z = self._f_T(z)
            # Total temporal transport density
            H_T = Hz_standard * (1 + self.epsilon_T * f_T_z)
        else:
            H_T = Hz_standard
            
        return H_T
    
    def _f_T(self, z: float) -> float:
        """TEP transition function.
        
        Returns value between 0 and 1, representing the strength of
        temporal shear effects at redshift z.
        
        For z << z_T: ~0 (weak TEP, near-LCDM)
        For z >> z_T: ~1 (strong TEP, temporal transport dominates)
        """
        return 1.0 - np.exp(-((z / self.z_T) ** self.n))
    
    def luminosity_distance(self, z: float) -> float:
        """Compute TEP luminosity distance.
        
        In TEP, the luminosity distance is modified because:
        1. Photon energy loss from temporal redshift
        2. Arrival rate dilation from temporal transport
        3. Geometric area factor
        
        The key TEP prediction is that D_L ≠ (1+z)^2 D_A in general.
        """
        # Compute comoving distance via integration
        # D_C = c ∫_0^z dz' / H_T(z')
        
        def integrand(zp):
            return self.c / self._H_T(zp)
        
        D_C, _ = integrate.quad(integrand, 0, z, limit=100)
        
        # TEP luminosity distance
        # In full TEP, this includes temporal transport corrections
        # For the mixed model (epsilon_T interpolates between LCDM and TEP):
        
        if self.epsilon_T > 0:
            # TEP-modified distance
            # The temporal transport affects both the comoving distance
            # AND the (1+z) factors
            
            # Effective redshift from temporal transport
            z_TEP = self._compute_tep_redshift(z)
            
            # TEP luminosity distance
            # D_L = (1 + z_TEP) * D_C when TEP effects are included
            D_L = (1 + z_TEP) * D_C
        else:
            # LCDM limit
            D_L = (1 + z) * D_C
            
        return D_L
    
    def _compute_tep_redshift(self, z_geom: float) -> float:
        """Compute the TEP effective redshift.
        
        In TEP, the observed redshift is:
            1 + z_obs = (1 + z_geom) * (1 + z_TEP)
            
        where z_geom is from standard expansion and z_TEP is from
        temporal transport along the path.
        
        For simplicity in this implementation:
            z_TEP ≈ epsilon_T * f_T(z_geom) * min(z_geom, z_T * 3)
            
        The min() prevents blow-up at very high redshift (z >> z_T),
        which is needed for acoustic consistency at CMB redshifts.
        """
        if self.epsilon_T == 0:
            return z_geom
            
        # TEP temporal redshift contribution
        # Cap z contribution to prevent blow-up at z >> z_T
        z_effective = min(z_geom, self.z_T * 3)
        z_TEP = self.epsilon_T * self._f_T(z_geom) * z_effective
        
        # Total effective redshift
        z_total = (1 + z_geom) * (1 + z_TEP) - 1
        
        return z_total
    
    def distance_modulus(self, z: float) -> float:
        """Compute distance modulus mu = 5*log10(D_L/10pc)."""
        D_L = self.luminosity_distance(z)  # Mpc
        mu = 5.0 * np.log10(D_L) + 25.0
        return mu
    
    def angular_diameter_distance(self, z: float) -> float:
        """Compute angular diameter distance.
        
        In TEP, D_A ≠ D_L / (1+z)^2 in general.
        This is a key testable prediction.
        """
        D_L = self.luminosity_distance(z)
        z_TEP = self._compute_tep_redshift(z)
        
        # TEP angular diameter distance
        # For the mixed model:
        D_A = D_L / (1 + z_TEP)**2
        
        return D_A
    
    def distance_duality(self, z: float) -> float:
        """Compute distance duality violation η(z) = D_L / [(1+z)^2 D_A].
        
        LCDM: η(z) = 1 exactly
        TEP: η(z) ≠ 1 in general (testable prediction)
        """
        D_L = self.luminosity_distance(z)
        D_A = self.angular_diameter_distance(z)
        
        if D_A <= 0:
            return np.nan
            
        eta = D_L / (D_A * (1 + z)**2)
        return eta
    
    def sn_stretch_factor(self, z: float) -> float:
        """Compute SN light curve stretch factor from temporal transport.
        
        In TEP, light curves are stretched by (1 + z_TEP):
            Δt_obs = (1 + z_TEP) * Δt_em
            
        This replaces the standard (1+z) time dilation.
        """
        if self.epsilon_T == 0:
            return 1 + z
        
        z_TEP = self._compute_tep_redshift(z)
        return 1 + z_TEP

    def temporal_enhancement_factor(self, z: float, Phi_ref: float = 0.0) -> float:
        """Compute temporal enhancement factor Γ_t from Kos manuscript.
        
        Γ_t = exp[K * (Φ - Φ_ref) / c^2 * sqrt(1+z)]
        
        This affects stellar age inference: t_eff = Γ_t * t_cosmic
        and mass-to-light ratio: M/L ∝ t_eff^n
        
        Args:
            z: Redshift
            Phi_ref: Reference gravitational potential (default 0)
            
        Returns:
            Γ_t: Temporal enhancement factor
        """
        if self.epsilon_T == 0:
            return 1.0
        
        # Estimate gravitational potential from matter density
        # Φ ≈ -GM/R ∝ -Ω_m * (1+z)
        # This is a simplified approximation
        Phi = -1e5 * self.Omega_m * (1 + z)  # Approximate in m^2/s^2
        
        # K is related to kappa_gal and epsilon_T
        # K ~ kappa_gal * epsilon_T / (1e6 mag) for dimensional consistency
        K = self.kappa_gal * self.epsilon_T / 1e6
        
        # Compute Γ_t
        Gamma_t = np.exp(K * (Phi - Phi_ref) / (self.c * 1000)**2 * np.sqrt(1 + z))
        
        return Gamma_t

    def cepheid_clock_correction(self, sigma: float, sigma_ref: float = 75.25) -> float:
        """Compute Cepheid clock correction from Kingston upon Hull manuscript.
        
        Δμ = κ_Cep * S(ρ) * (σ^2 - σ_ref^2) / c^2
        Reference: 11-TEP-H0-v0.6, Eq. 47-48

        This corrects for environment-dependent Cepheid period contraction.

        Args:
            sigma: Host galaxy velocity dispersion (km/s)
            sigma_ref: Reference velocity dispersion (default 75.25 km/s from SH0ES anchor weighting)
            
        Returns:
            Δμ: Distance modulus correction (mag)
        """
        if self.epsilon_T == 0:
            return 0.0
        
        # Screening function S(ρ) - simplified as environmental factor
        # In full TEP, this depends on local density and Temporal Shear
        S = 1.0  # Assume unscreened for now
        
        # Compute correction
        Delta_mu = self.kappa_gal * S * (sigma**2 - sigma_ref**2) / self.c**2
        
        return Delta_mu

    def pulsar_spin_down_correction(self, rho: float, rho_ref: float = 1e-18) -> float:
        """Compute pulsar spin-down correction from Caracas manuscript.
        
        This models the suppressed density scaling observed in globular cluster pulsars.
        
        Args:
            rho: Local density (g/cm^3)
            rho_ref: Reference density (default 1e-18 g/cm^3)
            
        Returns:
            Correction factor for spin-down rate
        """
        if self.epsilon_T == 0:
            return 1.0
        
        # Screening factor based on density
        # In full TEP, this uses Temporal Topology saturation density ρ_T
        # Reference: 0-TEP-v0.8-Jakarta, 10-TEP-COS-v0.6 (ρ_T ≈ 20 g/cm³ from terrestrial calibration)
        rho_T = 20.0  # g/cm^3 from Jakarta/Kilifi terrestrial calibration
        
        # Screening function
        S = 1.0 / (1 + (rho / rho_T))
        
        # Spin-down enhancement
        enhancement = 1.0 + self.kappa_msp * S * (rho - rho_ref) / rho_ref * 1e-6
        
        return enhancement

    def wide_binary_velocity_enhancement(self, separation: float) -> float:
        """Compute wide binary velocity enhancement from Kilifi manuscript.
        
        ṽ(s) = 1 + α_sat * (1 - exp(-s/R_s))
        
        Args:
            separation: Projected separation (AU)
            
        Returns:
            Velocity enhancement factor
        """
        if self.epsilon_T == 0:
            return 1.0
        
        # Screening radius R_s depends on mass and environment
        # For typical wide binary (M ~ 1.2 M_sun), R_s ~ 2646 AU
        # Reference: 13-TEP-WB-v0.3, R_s = 2,646 ± 182 AU from fit to 341,315 systems
        R_s = 2646.0  # AU from Kilifi fit result
        
        # Velocity enhancement
        enhancement = 1.0 + self.alpha_sat * (1 - np.exp(-separation / R_s))
        
        return enhancement

    def get_response_coefficients(self) -> Dict[str, float]:
        """Get current response coefficient values.
        
        Returns:
            Dictionary of response coefficients
        """
        return {
            "kappa_gal": self.kappa_gal,
            "kappa_msp": self.kappa_msp,
            "alpha_sat": self.alpha_sat,
            "epsilon_T": self.epsilon_T,
        }

    def check_cross_probe_consistency(self) -> Dict[str, bool]:
        """Check consistency between different probe response coefficients.
        
        According to the manuscripts, the response coefficients should be
        consistent across different probes:
        - κ_Cep ~ κ_MSP ~ 10^5-10^6
        - α_sat ~ 0.36 (dimensionless)
        
        Returns:
            Dictionary of consistency checks
        """
        checks = {}
        
        # Check that kappa_gal and kappa_msp are within an order of magnitude
        kappa_ratio = self.kappa_gal / self.kappa_msp
        checks["kappa_consistency"] = (0.1 < kappa_ratio < 10.0)
        
        # Check that alpha_sat is in expected range
        checks["alpha_sat_range"] = (0.3 < self.alpha_sat < 0.5)
        
        # Check that epsilon_T is small (consistent with CMB constraints)
        checks["epsilon_T_range"] = (self.epsilon_T < 0.01)
        
        return checks


class TEPCosmologyFitter:
    """Fit TEP cosmology to observational data."""
    
    def __init__(self, data_z: np.ndarray, data_mu: np.ndarray, data_mu_err: np.ndarray):
        """Initialize with Pantheon+ data.
        
        Args:
            data_z: Redshift values
            data_mu: Distance modulus observations
            data_mu_err: Distance modulus uncertainties
        """
        self.z = data_z
        self.mu = data_mu
        self.mu_err = data_mu_err
        
    def fit(self, initial_guess: Optional[dict] = None) -> dict:
        """Fit TEP model parameters to data.
        
        Returns best-fit parameters and fit statistics.
        """
        from scipy.optimize import minimize
        
        if initial_guess is None:
            initial_guess = {
                'H0': 70.0,
                'Omega_m': 0.3,
                'epsilon_T': 0.1,
                'z_T': 0.5,
            }
        
        def chi2(params):
            """Compute chi-squared for given parameters."""
            H0, Omega_m, epsilon_T, z_T = params
            
            # Create TEP model
            tep = TEPCosmology(
                H0=H0,
                Omega_m=Omega_m,
                Omega_L=1.0 - Omega_m,
                epsilon_T=epsilon_T,
                z_T=z_T,
            )
            
            # Compute model predictions
            mu_model = np.array([tep.distance_modulus(z) for z in self.z])
            
            # Chi-squared
            chi2_val = np.sum(((self.mu - mu_model) / self.mu_err) ** 2)
            
            return chi2_val
        
        # Initial parameter vector
        x0 = [
            initial_guess['H0'],
            initial_guess['Omega_m'],
            initial_guess['epsilon_T'],
            initial_guess['z_T'],
        ]
        
        # Parameter bounds
        bounds = [
            (50, 100),     # H0
            (0.1, 0.9),    # Omega_m
            (0.0, 2.0),    # epsilon_T
            (0.1, 5.0),    # z_T
        ]
        
        # Minimize chi-squared
        result = minimize(chi2, x0, bounds=bounds, method='L-BFGS-B')
        
        # Compute best-fit statistics
        H0_best, Omega_m_best, epsilon_T_best, z_T_best = result.x
        chi2_min = result.fun
        ndof = len(self.z) - 4  # 4 parameters
        
        return {
            'H0': H0_best,
            'Omega_m': Omega_m_best,
            'epsilon_T': epsilon_T_best,
            'z_T': z_T_best,
            'chi2': chi2_min,
            'ndof': ndof,
            'chi2_per_dof': chi2_min / ndof if ndof > 0 else chi2_min,
            'success': result.success,
        }
    
    def compare_models(self) -> dict:
        """Compare TEP, LCDM, and static metric models."""
        from core.cosmology import CosmologyFLRW
        from core.static_metric import StaticCosmology
        
        results = {}
        
        # 1. LCDM fit
        lcdm = CosmologyFLRW(H0=70, Om0=0.3, Ode0=0.7, Ok0=0.0)
        mu_lcdm = np.array([5.0 * np.log10(float(lcdm.luminosity_distance(np.array([z]))[0])) + 25.0 
                           for z in self.z])
        chi2_lcdm = np.sum(((self.mu - mu_lcdm) / self.mu_err) ** 2)
        
        results['lcdm'] = {
            'chi2': chi2_lcdm,
            'ndof': len(self.z) - 2,  # H0, Omega_m
            'n_params': 2,
        }
        
        # 2. TEP fit
        tep_fit = self.fit()
        tep = TEPCosmology(
            H0=tep_fit['H0'],
            Omega_m=tep_fit['Omega_m'],
            epsilon_T=tep_fit['epsilon_T'],
            z_T=tep_fit['z_T'],
        )
        mu_tep = np.array([tep.distance_modulus(z) for z in self.z])
        chi2_tep = np.sum(((self.mu - mu_tep) / self.mu_err) ** 2)
        
        results['tep'] = {
            'chi2': chi2_tep,
            'ndof': len(self.z) - 4,  # H0, Omega_m, epsilon_T, z_T
            'n_params': 4,
            'parameters': tep_fit,
        }
        
        # 3. Static metric
        static = StaticCosmology(H0=70, T0=2.725, L_c=3000, alpha=1.0, Omega_m=0.3, Omega_L=0.7)
        mu_static = np.array([5.0 * np.log10(static.luminosity_distance(z)) + 25.0 
                             for z in self.z])
        chi2_static = np.sum(((self.mu - mu_static) / self.mu_err) ** 2)
        
        results['static'] = {
            'chi2': chi2_static,
            'ndof': len(self.z) - 1,  # H0 only
            'n_params': 1,
        }
        
        # Model comparison statistics
        for model in ['lcdm', 'tep', 'static']:
            r = results[model]
            r['chi2_per_dof'] = r['chi2'] / r['ndof'] if r['ndof'] > 0 else r['chi2']
            r['aic'] = r['chi2'] + 2 * r['n_params']
            r['bic'] = r['chi2'] + r['n_params'] * np.log(len(self.z))
        
        # Determine best model
        bic_values = {m: results[m]['bic'] for m in ['lcdm', 'tep', 'static']}
        best_model = min(bic_values, key=bic_values.get)
        
        results['best_model'] = best_model
        results['delta_bic_tep_vs_lcdm'] = results['tep']['bic'] - results['lcdm']['bic']
        results['tep_competitive'] = results['delta_bic_tep_vs_lcdm'] < 10  # Within 10 BIC points
        
        return results


if __name__ == "__main__":
    # Test TEP cosmology
    print("Testing TEP Cosmology...")
    
    tep = TEPCosmology(H0=70, Omega_m=0.3, epsilon_T=0.5, z_T=0.5)
    
    for z in [0.1, 0.5, 1.0, 2.0]:
        D_L = tep.luminosity_distance(z)
        D_A = tep.angular_diameter_distance(z)
        eta = tep.distance_duality(z)
        mu = tep.distance_modulus(z)
        
        print(f"z={z:.1f}: D_L={D_L:.2f} Mpc, D_A={D_A:.2f} Mpc, eta={eta:.4f}, mu={mu:.2f}")
    
    print("\nTEP Cosmology implementation complete!")
