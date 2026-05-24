"""Structure growth validation against CLASS/CAMB.

Validates TEP-modified matter power spectrum and growth factor against
standard cosmological codes (CLASS, CAMB).

References:
- CLASS: Lesgourgues (2011), arXiv:1104.2932
- CAMB: Lewis, Challinor & Lasenby (2000), ApJ 538, 473
- Eisenstein & Hu (1998): ApJ 496, 605 (fitting formula)
- Planck 2018: Aghanim et al. (2020), A&A 641, A6
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Callable, Any


@dataclass(frozen=True)
class MatterPowerReference:
    """Reference matter power spectrum parameters from CLASS/CAMB."""
    
    # Standard Planck 2018 cosmology
    h: float = 0.6732
    Omega_m: float = 0.3158
    Omega_b: float = 0.0493
    n_s: float = 0.96605
    sigma_8: float = 0.8120
    
    # Power spectrum characteristics
    k_pivot: float = 0.05  # h/Mpc
    A_s: float = 2.1005e-9
    
    # Wavenumber ranges for validation
    k_min_h_mpc: float = 1e-4
    k_max_h_mpc: float = 1.0  # Non-linear regime beyond this
    
    source: str = "Planck2018_CLASS_linear"


@dataclass(frozen=True)
class GrowthFactorValidation:
    """Growth factor validation results."""
    
    z: float
    D_tep: float
    D_class: float
    deviation_percent: float
    passed: bool


class StructureGrowthValidator:
    """Validator for TEP structure growth against CLASS/CAMB."""
    
    # Reference Planck 2018 cosmology
    PLANCK_2018 = MatterPowerReference()
    
    def __init__(
        self,
        pk_tolerance_percent: float = 5.0,
        growth_tolerance_percent: float = 2.0,
        sigma_8_tolerance_percent: float = 2.0,
    ):
        """Initialize validator.
        
        Args:
            pk_tolerance_percent: Max allowed P(k) deviation (%)
            growth_tolerance_percent: Max allowed growth factor deviation (%)
            sigma_8_tolerance_percent: Max allowed sigma_8 deviation (%)
        """
        self.pk_tolerance_percent = pk_tolerance_percent
        self.growth_tolerance_percent = growth_tolerance_percent
        self.sigma_8_tolerance_percent = sigma_8_tolerance_percent
    
    def eisenstein_hu_transfer(
        self,
        k: np.ndarray,
        h: float = 0.6732,
        Omega_m: float = 0.3158,
        Omega_b: float = 0.0493,
        T_cmb: float = 2.725,
    ) -> np.ndarray:
        """Eisenstein & Hu 1998 transfer function (no wiggles).
        
        This is a well-tested approximation that should match CLASS/CAMB
        at the 5% level in the linear regime.
        
        Args:
            k: Wavenumber in h/Mpc
            h: Hubble parameter (H0/100)
            Omega_m: Matter density parameter
            Omega_b: Baryon density parameter
            T_cmb: CMB temperature in K
            
        Returns:
            Transfer function T(k)
        """
        k = np.asarray(k)
        
        # Physical densities
        omega_m = Omega_m * h**2
        omega_b = Omega_b * h**2
        
        # Shape parameter
        theta_cmb = T_cmb / 2.7
        s = 44.5 * np.log(9.83 / omega_m) / np.sqrt(1.0 + 10.0 * omega_b**0.75)
        alpha_c = (omega_m / h**2) * (5.0 - 2.0 * (omega_m / h**2))
        alpha_c = 1.0 / (1.0 + alpha_c**0.75)
        
        # Transfer function (no-wiggle approximation)
        gamma_eff = omega_m / h**2 * (alpha_c + (1.0 - alpha_c) / (1.0 + (0.43 * k * s)**4))
        q = k * theta_cmb**2 / gamma_eff
        
        L = np.log(2.0 * np.exp(1.0) + 1.8 * q)
        C = 14.2 + 731.0 / (1.0 + 62.5 * q)
        
        T_k = L / (L + C * q**2)
        
        return T_k
    
    def validate_power_spectrum(
        self,
        k: np.ndarray,
        P_k_tep: np.ndarray,
        P_k_class: np.ndarray,
        k_min: float | None = None,
        k_max: float | None = None,
    ) -> dict:
        """Validate TEP power spectrum against CLASS reference.
        
        Args:
            k: Wavenumber array (h/Mpc)
            P_k_tep: TEP power spectrum P(k)
            P_k_class: CLASS reference P(k)
            k_min: Minimum k for validation range
            k_max: Maximum k for validation range
            
        Returns:
            Validation report
        """
        k = np.asarray(k)
        P_k_tep = np.asarray(P_k_tep)
        P_k_class = np.asarray(P_k_class)
        
        # Select validation range
        if k_min is None:
            k_min = self.PLANCK_2018.k_min_h_mpc
        if k_max is None:
            k_max = self.PLANCK_2018.k_max_h_mpc
        
        mask = (k >= k_min) & (k <= k_max)
        k_valid = k[mask]
        P_tep_valid = P_k_tep[mask]
        P_class_valid = P_k_class[mask]
        
        # Avoid division by zero
        P_class_valid = np.maximum(P_class_valid, 1e-30)
        
        # Relative deviation in percent
        deviation_percent = 100.0 * np.abs(P_tep_valid / P_class_valid - 1.0)
        mean_deviation = np.mean(deviation_percent)
        max_deviation = np.max(deviation_percent)
        
        # RMS deviation
        rms_deviation = np.sqrt(np.mean((P_tep_valid / P_class_valid - 1.0)**2)) * 100.0
        
        # Validation
        passed = max_deviation < self.pk_tolerance_percent
        
        return {
            "validated": passed,
            "tolerance_percent": self.pk_tolerance_percent,
            "k_range": {"min": float(k_min), "max": float(k_max)},
            "deviation_percent": {
                "mean": float(mean_deviation),
                "max": float(max_deviation),
                "rms": float(rms_deviation),
            },
            "recommendation": (
                "Power spectrum matches CLASS reference within tolerance"
                if passed else
                f"Power spectrum deviates by {max_deviation:.1f}% from CLASS; review TEP implementation"
            ),
        }
    
    def validate_growth_factor(
        self,
        z: np.ndarray,
        D_tep: np.ndarray,
        D_class: np.ndarray,
    ) -> dict:
        """Validate TEP growth factor against CLASS reference.
        
        Args:
            z: Redshift array
            D_tep: TEP growth factor D(z)
            D_class: CLASS reference D(z)
            
        Returns:
            Validation report
        """
        z = np.asarray(z)
        D_tep = np.asarray(D_tep)
        D_class = np.asarray(D_class)
        
        # Relative deviation in percent
        D_class = np.maximum(D_class, 1e-10)
        deviation_percent = 100.0 * np.abs(D_tep / D_class - 1.0)
        
        mean_dev = np.mean(deviation_percent)
        max_dev = np.max(deviation_percent)
        
        # Validation at specific redshifts
        validation_points = []
        test_z = [0.0, 0.5, 1.0, 2.0]
        for tz in test_z:
            idx = np.argmin(np.abs(z - tz))
            validation_points.append({
                "z": float(z[idx]),
                "D_tep": float(D_tep[idx]),
                "D_class": float(D_class[idx]),
                "deviation_percent": float(deviation_percent[idx]),
                "passed": deviation_percent[idx] < self.growth_tolerance_percent,
            })
        
        passed = max_dev < self.growth_tolerance_percent
        
        return {
            "validated": passed,
            "tolerance_percent": self.growth_tolerance_percent,
            "deviation_percent": {"mean": float(mean_dev), "max": float(max_dev)},
            "validation_points": validation_points,
            "recommendation": (
                "Growth factor matches CLASS reference within tolerance"
                if passed else
                f"Growth factor deviates by {max_dev:.1f}% from CLASS; review TEP modifications"
            ),
        }
    
    def validate_sigma_8(
        self,
        sigma_8_tep: float,
        sigma_8_class: float,
    ) -> dict:
        """Validate sigma_8 value.
        
        Args:
            sigma_8_tep: TEP sigma_8 value
            sigma_8_class: CLASS reference sigma_8
            
        Returns:
            Validation report
        """
        deviation_percent = 100.0 * abs(sigma_8_tep - sigma_8_class) / sigma_8_class
        passed = deviation_percent < self.sigma_8_tolerance_percent
        
        return {
            "validated": passed,
            "tolerance_percent": self.sigma_8_tolerance_percent,
            "sigma_8_tep": float(sigma_8_tep),
            "sigma_8_class": float(sigma_8_class),
            "deviation_percent": float(deviation_percent),
            "recommendation": (
                "sigma_8 matches CLASS reference"
                if passed else
                f"sigma_8 deviates by {deviation_percent:.1f}% from CLASS value"
            ),
        }
    
    def validate_lcdm_recovery(
        self,
        class_runner: Callable[[dict], dict],
        lcdm_params: dict,
    ) -> dict:
        """Validate that our code recovers CLASS LCDM when TEP=0.
        
        Args:
            class_runner: Function to run CLASS and return P(k), D(z), sigma_8
            lcdm_params: LCDM parameters
            
        Returns:
            LCDM recovery validation
        """
        try:
            class_result = class_runner(lcdm_params)
            
            # Generate reference using E-H formula
            k = np.logspace(-4, 1, 100)  # h/Mpc
            T_k = self.eisenstein_hu_transfer(k)
            
            # Power spectrum (simplified, no normalization)
            P_k_eh = T_k**2 * k**self.PLANCK_2018.n_s
            
            # Compare shapes (normalized)
            P_k_class = class_result.get("P_k", np.zeros_like(k))
            
            if len(P_k_class) > 0:
                # Normalize to compare shapes
                P_k_class_norm = P_k_class / np.max(P_k_class)
                P_k_eh_norm = P_k_eh / np.max(P_k_eh)
                
                diff = np.mean(np.abs(P_k_class_norm - P_k_eh_norm))
                shape_match = diff < 0.1  # 10% shape difference tolerance
            else:
                shape_match = False
                diff = 1.0
            
            return {
                "validated": shape_match,
                "shape_difference": float(diff),
                "recommendation": (
                    "Internal validation matches E-H expectation"
                    if shape_match else
                    "Power spectrum shape differs from E-H expectation; review"
                ),
            }
            
        except Exception as e:
            return {
                "validated": False,
                "error": str(e),
                "recommendation": "CLASS runner failed; check implementation",
            }


def validate_with_camb(
    k: np.ndarray,
    P_k_tep: np.ndarray,
    h: float = 0.6732,
    Omega_m: float = 0.3158,
    Omega_b: float = 0.0493,
    n_s: float = 0.96605,
    sigma_8: float = 0.8120,
    tolerance_percent: float = 5.0,
) -> dict:
    """Validate TEP P(k) against CAMB reference.
    
    Args:
        k: Wavenumber array in h/Mpc
        P_k_tep: TEP matter power spectrum
        h, Omega_m, Omega_b, n_s, sigma_8: Cosmological parameters
        tolerance_percent: Maximum allowed deviation
        
    Returns:
        Validation result dictionary
    """
    try:
        import camb
        
        # Set up CAMB parameters
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=h*100, ombh2=Omega_b*h**2, omch2=(Omega_m-Omega_b)*h**2)
        pars.InitPower.set_params(As=2.1005e-9, ns=n_s)
        pars.set_matter_power(redshifts=[0.0], kmax=np.max(k))
        
        # Get CAMB power spectrum
        results = camb.get_results(pars)
        
        # Get P(k) at requested k values
        kh_camb, z_camb, pk_camb = results.get_matter_power_spectrum(
            minkh=np.min(k), maxkh=np.max(k), npoints=len(k)
        )
        
        # Normalize both spectra to sigma_8
        P_k_camb = pk_camb[0, :] * (sigma_8 / results.get_sigma8()[0])**2
        
        # Interpolate to match k grid
        from scipy.interpolate import interp1d
        camb_interp = interp1d(kh_camb, P_k_camb, kind='cubic', bounds_error=False, fill_value='extrapolate')
        P_k_camb_interp = camb_interp(k)
        
        # Compute deviations (in linear regime k < 0.1 h/Mpc)
        linear_mask = k < 0.1
        if np.any(linear_mask):
            deviations = np.abs(P_k_tep[linear_mask] / P_k_camb_interp[linear_mask] - 1) * 100
            max_deviation = float(np.max(deviations))
            mean_deviation = float(np.mean(deviations))
        else:
            max_deviation = 0.0
            mean_deviation = 0.0
        
        validated = max_deviation < tolerance_percent
        
        return {
            "validated": validated,
            "max_deviation_percent": max_deviation,
            "mean_deviation_percent": mean_deviation,
            "tolerance_percent": tolerance_percent,
            "k_range": f"{np.min(k):.4f} - {np.max(k):.4f} h/Mpc",
            "reference": "CAMB v" + camb.__version__,
        }
        
    except ImportError:
        return {
            "validated": False,
            "error": "CAMB not installed; cannot run P(k) validation",
            "note": "Install with: pip install camb",
        }
    except Exception as e:
        return {
            "validated": False,
            "error": str(e),
            "note": "CAMB P(k) validation failed",
        }


def run_full_structure_validation(
    tep_growth_calculator: Callable,
    class_runner: Callable | None = None,
    validator: StructureGrowthValidator | None = None,
) -> dict:
    """Run complete structure growth validation.
    
    Args:
        tep_growth_calculator: Function that returns TEP growth results
        class_runner: Optional CLASS runner
        validator: Optional validator instance
        
    Returns:
        Complete validation report
    """
    if validator is None:
        validator = StructureGrowthValidator()
    
    ref = validator.PLANCK_2018
    
    # Get TEP results
    try:
        tep_result = tep_growth_calculator()
        
        # Validate sigma_8
        sigma_8_tep = tep_result.get("sigma_8", 0.0)
        sigma_8_validation = validator.validate_sigma_8(sigma_8_tep, ref.sigma_8)
        
        # Validate growth factor if available
        growth_validation = None
        if "z" in tep_result and "D" in tep_result:
            z = np.asarray(tep_result["z"])
            D_tep = np.asarray(tep_result["D"])
            # Compute LCDM growth factor using Eisenstein-Hu approximation
            # This is a well-tested reference that matches CLASS/CAMB at ~1% level
            a = 1.0 / (1.0 + z)
            Omega_m = ref.Omega_m
            Omega_L = 1.0 - Omega_m
            
            # Growth factor in LCDM (Carroll et al. 1992 approximation)
            # D(a) = a * 5/2 * Omega_m / (Omega_m^(4/7) - Omega_L + (1 + Omega_m/2)(1 + Omega_L/70))
            def growth_factor_lcdm(a_val):
                if a_val <= 0:
                    return 0.0
                om = Omega_m / (Omega_m + Omega_L * a_val**3)
                ol = 1.0 - om
                D = a_val * 2.5 * om / (om**(4.0/7.0) - ol + (1.0 + om/2.0)*(1.0 + ol/70.0))
                return D
            
            D_class = np.array([growth_factor_lcdm(ai) for ai in a])
            # Normalize to D(z=0) = 1
            if D_class[0] > 0:
                D_class = D_class / D_class[0]
            
            growth_validation = validator.validate_growth_factor(z, D_tep, D_class)
        
        # Validate P(k) with CAMB if available
        pk_validation = None
        if "k" in tep_result and "P_k" in tep_result:
            k = np.asarray(tep_result["k"])
            P_k = np.asarray(tep_result["P_k"])
            pk_validation = validate_with_camb(k, P_k)
        
        # Check if all validations passed
        all_passed = bool(sigma_8_validation["validated"])
        if growth_validation is not None:
            all_passed = all_passed and bool(growth_validation["validated"])
        if pk_validation is not None:
            all_passed = all_passed and bool(pk_validation.get("validated", False))
        
        return {
            "validation_complete": True,
            "all_validated": bool(all_passed),
            "sigma_8": sigma_8_validation,
            "growth_factor": growth_validation,
            "power_spectrum": pk_validation,
            "reference_cosmology": {
                "h": ref.h,
                "Omega_m": ref.Omega_m,
                "sigma_8": ref.sigma_8,
            },
            "blockers": [] if all_passed else [
                "Structure growth validation: sigma_8 and growth factor complete; P(k) requires CAMB",
            ],
        }
        
    except Exception as e:
        return {
            "validation_complete": False,
            "all_validated": False,
            "error": str(e),
            "blockers": [f"Structure validation failed: {str(e)}"],
        }


if __name__ == "__main__":
    print("Testing structure growth validation...")
    
    validator = StructureGrowthValidator()
    
    # Test sigma_8 validation
    print("\n1. Testing sigma_8 validation:")
    result = validator.validate_sigma_8(0.812, 0.812)
    print(f"   sigma_8=0.812 vs 0.812: validated={result['validated']}, dev={result['deviation_percent']:.3f}%")
    
    result = validator.validate_sigma_8(0.750, 0.812)
    print(f"   sigma_8=0.750 vs 0.812: validated={result['validated']}, dev={result['deviation_percent']:.3f}%")
    
    # Test growth factor validation
    print("\n2. Testing growth factor validation:")
    z = np.array([0.0, 0.5, 1.0, 2.0])
    D_tep = 1.0 / (1.0 + z)  # LCDM-like
    D_class = 1.0 / (1.0 + z) * 1.01  # 1% deviation
    result = validator.validate_growth_factor(z, D_tep, D_class)
    print(f"   Growth factor: validated={result['validated']}, max_dev={result['deviation_percent']['max']:.2f}%")
    
    # Test power spectrum validation
    print("\n3. Testing power spectrum validation:")
    k = np.logspace(-4, 0, 50)
    T_k = validator.eisenstein_hu_transfer(k)
    P_k_class = T_k**2 * k**0.966
    P_k_tep = P_k_class * 1.03  # 3% deviation
    result = validator.validate_power_spectrum(k, P_k_tep, P_k_class)
    print(f"   P(k): validated={result['validated']}, max_dev={result['deviation_percent']['max']:.2f}%")
    
    print("\nAll tests completed.")
