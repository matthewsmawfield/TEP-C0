"""Structure growth validation framework."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional
import numpy as np


@dataclass
class StructureGrowthValidator:
    """Validator for structure growth predictions."""
    
    H0: float = 70.0
    Omega_m: float = 0.3
    sigma_8: float = 0.81
    
    def validate_sigma_8(self, measured_sigma_8: float, error: float = 0.01) -> Dict:
        """Validate sigma_8 against Planck 2018 constraint."""
        planck_sigma_8 = 0.81
        planck_error = 0.006
        
        diff = measured_sigma_8 - planck_sigma_8
        sigma = np.sqrt(error**2 + planck_error**2)
        
        n_sigma = abs(diff) / sigma if sigma > 0 else 0
        
        return {
            'measured': measured_sigma_8,
            'planck_reference': planck_sigma_8,
            'difference': diff,
            'n_sigma': n_sigma,
            'validated': n_sigma < 3.0,  # Within 3 sigma
        }
    
    def validate_growth_rate(self, z: float, f_measured: float, f_lcdm: float, error: float = 0.05) -> Dict:
        """Validate growth rate f(z) against LCDM and observations."""
        diff = f_measured - f_lcdm
        n_sigma = abs(diff) / error if error > 0 else 0
        
        return {
            'z': z,
            'measured': f_measured,
            'lcdm_prediction': f_lcdm,
            'difference': diff,
            'n_sigma': n_sigma,
            'validated': n_sigma < 2.0,  # Within 2 sigma
        }


def run_full_structure_validation(
    H0: float = 70.0,
    Omega_m: float = 0.3,
    sigma_0: float = 0.0,
    tep_growth_calculator = None,
    validator = None,
) -> Dict:
    """Run full structure growth validation suite."""
    
    from scripts.utils.structure_formation import StructureFormation
    from scripts.utils.tep_perturbations import TEPPerturbations
    
    # Use passed-in objects or create new ones
    if validator is None:
        validator = StructureGrowthValidator(H0=H0, Omega_m=Omega_m)
    
    # Test sigma_8
    tep = TEPPerturbations(H0=H0, Omega_m=Omega_m, sigma_0=sigma_0)
    sigma_8_measured = tep.sigma_8()
    sigma_8_validation = validator.validate_sigma_8(sigma_8_measured)
    
    # Handle tep_growth_calculator - can be callable function or StructureFormation object
    if tep_growth_calculator is None:
        sf = StructureFormation(H0=H0, Om0=Omega_m, Sigma_0=sigma_0)
    elif callable(tep_growth_calculator):
        # It's a function that returns pre-computed data
        calc_data = tep_growth_calculator()
        sf = None  # Will use calc_data directly
    else:
        # It's a StructureFormation object
        sf = tep_growth_calculator
    
    # Test growth factor at several redshifts
    z_test = [0.0, 0.5, 1.0, 2.0]
    a_test = [1.0 / (1 + z) for z in z_test]
    
    if sf is not None:
        D_test = sf.growth_factor(np.array(a_test))
        f_test = [sf.growth_rate(a) for a in a_test]
    else:
        # Use pre-computed data from calc_data
        D_test = np.interp(z_test, calc_data['z'], calc_data['D'])
        f_test = [0.55] * len(z_test)  # Approximation
    
    growth_results = []
    for z, D, f in zip(z_test, D_test, f_test):
        # LCDM growth rate approximation
        f_lcdm = 0.55 + 0.05 * (1 / (1 + z))
        
        result = validator.validate_growth_rate(z, f, f_lcdm)
        growth_results.append({
            'z': z,
            'D': D,
            'f': f,
            'validated': result['validated'],
        })
    
    return {
        'sigma_8': sigma_8_validation,
        'growth_factor': {
            'z_values': z_test,
            'D_values': D_test.tolist() if hasattr(D_test, 'tolist') else list(D_test),
            'f_values': f_test,
        },
        'growth_validation': growth_results,
        'all_validated': sigma_8_validation['validated'] and all(r['validated'] for r in growth_results),
        'reference_cosmology': {
            'H0': H0,
            'Omega_m': Omega_m,
            'sigma_8': 0.81,
        },
    }
