"""CMB TEP CLASS parameter validation.

Validates TEP-modified CLASS predictions against standard CLASS/LCDM reference.
Provides validation gates for CMB spectra (TT, TE, EE) and derived parameters.

References:
- CLASS: Lesgourgues (2011), arXiv:1104.2932
- Planck 2018: Aghanim et al. (2020), A&A 641, A6
- Planck 2018 parameters: Table 1, TT,TE,EE+lowE+lensing
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Callable, Any


@dataclass(frozen=True)
class CLASSReferenceParameters:
    """Standard CLASS/Planck 2018 reference parameters."""
    
    # Cosmological parameters from Planck 2018 TT,TE,EE+lowE+lensing
    h: float = 0.6732
    omega_b: float = 0.022383
    omega_cdm: float = 0.12011
    A_s: float = 2.1005e-9
    n_s: float = 0.96605
    tau_reio: float = 0.0543
    N_ur: float = 2.0328  # 3.046 massless neutrinos, but 1 massive
    # Omega_k = 0 (flat)
    
    # Derived parameters
    H0: float = 67.32  # km/s/Mpc
    Omega_m: float = 0.3158
    Omega_Lambda: float = 0.6842
    sigma_8: float = 0.8120
    
    # CMB acoustic scale
    theta_s: float = 1.04109  # 100 * theta_star
    r_drag: float = 147.09  # Mpc
    
    source: str = "Planck2018_TTTEEE_lowE_lensing"


@dataclass(frozen=True)
class CMBSpectrumValidation:
    """CMB spectrum validation results."""
    
    spectrum_type: str  # "TT", "TE", "EE"
    l_max: int
    chi2_per_dof: float
    max_deviation_sigma: float
    passed: bool
    

class CLASSValidator:
    """Validator for TEP CLASS predictions against reference."""
    
    # Planck 2018 TT,TE,EE+lowE+lensing best-fit power spectrum
    # From publicly available CLASS/Planck 2018 runs
    PLANCK_2018 = CLASSReferenceParameters()
    
    def __init__(self, tolerance_sigma: float = 5.0, chi2_tolerance: float = 1.5):
        """Initialize validator.
        
        Args:
            tolerance_sigma: Max allowed deviation in single bin (sigma)
            chi2_tolerance: Max allowed chi2/dof for spectrum comparison
        """
        self.tolerance_sigma = tolerance_sigma
        self.chi2_tolerance = chi2_tolerance
    
    def validate_spectrum(
        self,
        l: np.ndarray,
        Dl_tep: np.ndarray,
        Dl_lcdm: np.ndarray,
        spectrum_type: str = "TT",
        cosmic_variance: np.ndarray | None = None,
    ) -> CMBSpectrumValidation:
        """Validate TEP spectrum against LCDM reference.
        
        Args:
            l: Multipole array
            Dl_tep: TEP D_ell values (uK^2)
            Dl_lcdm: LCDM D_ell values (uK^2)
            spectrum_type: "TT", "TE", or "EE"
            cosmic_variance: Optional cosmic variance estimate
            
        Returns:
            Validation result
        """
        l = np.asarray(l)
        Dl_tep = np.asarray(Dl_tep)
        Dl_lcdm = np.asarray(Dl_lcdm)
        
        # Use cosmic variance or 1% floor as uncertainty
        if cosmic_variance is not None:
            sigma = np.asarray(cosmic_variance)
        else:
            # Approximate cosmic variance for TT
            if spectrum_type == "TT":
                sigma = np.sqrt(2.0 / (2.0 * l + 1.0)) * Dl_lcdm
            else:
                sigma = 0.01 * Dl_lcdm  # 1% floor for TE/EE
        
        sigma = np.maximum(sigma, 1e-10 * Dl_lcdm)  # Avoid division by zero
        
        # Chi2 per degree of freedom
        residuals = (Dl_tep - Dl_lcdm) / sigma
        chi2 = np.sum(residuals**2)
        dof = len(l)
        chi2_per_dof = chi2 / dof if dof > 0 else 0.0
        
        # Maximum deviation
        max_dev = np.max(np.abs(residuals))
        
        # Validation criteria
        passed = (
            chi2_per_dof < self.chi2_tolerance and
            max_dev < self.tolerance_sigma
        )
        
        return CMBSpectrumValidation(
            spectrum_type=spectrum_type,
            l_max=int(np.max(l)),
            chi2_per_dof=float(chi2_per_dof),
            max_deviation_sigma=float(max_dev),
            passed=passed,
        )
    
    def validate_derived_parameters(
        self,
        H0: float,
        Omega_m: float,
        sigma_8: float,
        theta_s: float | None = None,
        r_drag: float | None = None,
    ) -> dict:
        """Validate derived cosmological parameters.
        
        Args:
            H0: Hubble constant (km/s/Mpc)
            Omega_m: Matter density parameter
            sigma_8: Matter fluctuation amplitude
            theta_s: Acoustic angular scale (100*theta_star, optional)
            r_drag: Sound horizon at drag epoch (Mpc, optional)
            
        Returns:
            Validation report
        """
        ref = self.PLANCK_2018
        
        # Uncertainties from Planck 2018 (approximate 1-sigma)
        H0_unc = 0.54  # km/s/Mpc
        Omega_m_unc = 0.007
        sigma_8_unc = 0.0073
        theta_s_unc = 0.00030
        r_drag_unc = 0.26
        
        results = {}
        
        def check_param(name: str, val: float, ref_val: float, unc: float) -> dict:
            dev = abs(val - ref_val) / unc if unc > 0 else 0.0
            return {
                "value": float(val),
                "reference": float(ref_val),
                "deviation_sigma": float(dev),
                "passed": dev < self.tolerance_sigma,
            }
        
        results["H0"] = check_param("H0", H0, ref.H0, H0_unc)
        results["Omega_m"] = check_param("Omega_m", Omega_m, ref.Omega_m, Omega_m_unc)
        results["sigma_8"] = check_param("sigma_8", sigma_8, ref.sigma_8, sigma_8_unc)
        
        if theta_s is not None:
            results["theta_s"] = check_param("theta_s", theta_s, ref.theta_s, theta_s_unc)
        if r_drag is not None:
            results["r_drag"] = check_param("r_drag", r_drag, ref.r_drag, r_drag_unc)
        
        all_passed = all(r["passed"] for r in results.values())
        
        return {
            "validated": all_passed,
            "tolerance_sigma": self.tolerance_sigma,
            "parameters": results,
            "reference_source": ref.source,
        }
    
    def validate_class_runner(
        self,
        class_runner: Callable[[dict], dict],
        tep_params: dict,
        lcdm_params: dict | None = None,
    ) -> dict:
        """Validate by running CLASS with TEP and LCDM parameters.
        
        Args:
            class_runner: Function that takes CLASS params dict and returns spectra
            tep_params: TEP parameter dictionary for CLASS
            lcdm_params: LCDM parameter dictionary (uses Planck 2018 defaults if None)
            
        Returns:
            Complete validation report
        """
        if lcdm_params is None:
            ref = self.PLANCK_2018
            lcdm_params = {
                "h": ref.h,
                "omega_b": ref.omega_b,
                "omega_cdm": ref.omega_cdm,
                "A_s": ref.A_s,
                "n_s": ref.n_s,
                "tau_reio": ref.tau_reio,
            }
        
        try:
            # Run CLASS for LCDM
            lcdm_result = class_runner(lcdm_params)
            
            # Run CLASS for TEP (with Sigma_0=0 should recover LCDM)
            tep_result = class_runner(tep_params)
            
            # Validate spectra
            validations = []
            for spectrum_type in ["TT", "TE", "EE"]:
                l_key = f"ell_{spectrum_type.lower()}"
                Dl_lcdm_key = f"cl_{spectrum_type.lower()}"
                Dl_tep_key = f"cl_{spectrum_type.lower()}"
                
                if l_key in lcdm_result and Dl_lcdm_key in lcdm_result:
                    val = self.validate_spectrum(
                        l=lcdm_result[l_key],
                        Dl_tep=tep_result[Dl_tep_key],
                        Dl_lcdm=lcdm_result[Dl_lcdm_key],
                        spectrum_type=spectrum_type,
                    )
                    validations.append(val)
            
            # Validate derived parameters
            derived = self.validate_derived_parameters(
                H0=tep_result.get("H0", 0.0),
                Omega_m=tep_result.get("Omega_m", 0.0),
                sigma_8=tep_result.get("sigma_8", 0.0),
            )
            
            all_spectra_passed = all(v.passed for v in validations)
            
            return {
                "validated": all_spectra_passed and derived["validated"],
                "spectra_validations": [
                    {
                        "type": v.spectrum_type,
                        "l_max": v.l_max,
                        "chi2_per_dof": v.chi2_per_dof,
                        "max_deviation_sigma": v.max_deviation_sigma,
                        "passed": v.passed,
                    }
                    for v in validations
                ],
                "derived_parameters": derived,
                "lcdm_params_used": lcdm_params,
                "tep_params_used": tep_params,
                "recommendation": (
                    "TEP CLASS predictions are consistent with Planck 2018 reference"
                    if (all_spectra_passed and derived["validated"]) else
                    "TEP CLASS predictions show significant deviations; review parameter mapping"
                ),
            }
            
        except Exception as e:
            return {
                "validated": False,
                "error": str(e),
                "recommendation": "CLASS runner failed; check implementation and CLASS installation",
            }


def validate_lcdm_limit(
    class_runner: Callable[[dict], dict],
    lcdm_params: dict,
    tep_lcdm_equivalent_params: dict,
    tolerance_sigma: float = 3.0,
) -> dict:
    """Validate that TEP with Sigma_0=0 recovers standard LCDM.
    
    This is the critical validation: TEP must reduce to LCDM when the shear
    parameter is zero.
    
    Args:
        class_runner: CLASS runner function
        lcdm_params: Standard LCDM parameters
        tep_lcdm_equivalent_params: TEP parameters with sigma_0=0
        tolerance_sigma: Acceptance threshold
        
    Returns:
        LCDM limit validation report
    """
    validator = CLASSValidator(tolerance_sigma=tolerance_sigma)
    
    result = validator.validate_class_runner(
        class_runner=class_runner,
        tep_params=tep_lcdm_equivalent_params,
        lcdm_params=lcdm_params,
    )
    
    # Add explicit LCDM limit test
    result["lcdm_limit_test"] = {
        "description": "TEP with sigma_0=0 must recover LCDM exactly",
        "passed": result["validated"],
        "tolerance_sigma": tolerance_sigma,
    }
    
    return result


if __name__ == "__main__":
    # Test validation with mock data
    validator = CLASSValidator()
    
    print("Testing CLASS parameter validation...")
    
    # Test derived parameter validation
    result = validator.validate_derived_parameters(
        H0=67.32,
        Omega_m=0.3158,
        sigma_8=0.8120,
    )
    print(f"Derived parameters validated: {result['validated']}")
    for param, details in result["parameters"].items():
        print(f"  {param}: dev={details['deviation_sigma']:.2f} sigma, passed={details['passed']}")
    
    # Test with deviant parameters
    print("\nTesting with deviant parameters:")
    result = validator.validate_derived_parameters(
        H0=75.0,  # Too high
        Omega_m=0.25,  # Too low
        sigma_8=0.70,  # Too low
    )
    print(f"Derived parameters validated: {result['validated']}")
    for param, details in result["parameters"].items():
        print(f"  {param}: dev={details['deviation_sigma']:.2f} sigma, passed={details['passed']}")
    
    # Test spectrum validation
    print("\nTesting spectrum validation...")
    l = np.arange(2, 2509)
    # Mock LCDM spectrum (simple approximation)
    Dl_lcdm = 1000.0 * np.exp(-((l - 220)**2) / (2 * 50**2)) + 1000.0
    # TEP spectrum with small deviation
    Dl_tep = Dl_lcdm * (1.0 + 0.001 * np.sin(l / 100.0))
    
    result = validator.validate_spectrum(l, Dl_tep, Dl_lcdm, "TT")
    print(f"TT spectrum: chi2/dof={result.chi2_per_dof:.3f}, max_dev={result.max_deviation_sigma:.2f} sigma, passed={result.passed}")
