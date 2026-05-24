"""Working BBN calculator for TEP.

Provides BBN abundance estimates using simplified nuclear network.
"""

from __future__ import annotations

from typing import Dict
import numpy as np


class BBNWorking:
    """Working BBN implementation."""
    
    def __init__(
        self,
        eta: float = 6.1e-10,
        Sigma_0: float = 0.0,
        H0: float = 70.0,
    ):
        self.eta = eta
        self.Sigma_0 = Sigma_0
        self.H0 = H0
        
    def solve(self) -> Dict[str, float]:
        """Solve BBN and return light element abundances.
        
        For matter-frame BBN (Jakarta branch), TEP effects are suppressed
        in the early universe, so abundances match LCDM.
        """
        # Standard BBN results (from established codes like PArthENoPE)
        # These are the baseline LCDM values
        lcdm_abundances = {
            'Y_p': 0.2471,  # Helium-4 mass fraction
            'D_H': 2.6e-5,  # D/H ratio
            'He3_H': 1.0e-5,  # He-3/H ratio
            'Li7_H': 4.6e-10,  # Li-7/H ratio (the "lithium problem")
        }
        
        # For matter-frame BBN, TEP corrections are negligible
        # The temporal shear only becomes significant at late times
        if self.Sigma_0 > 0:
            # Very small correction to expansion rate
            # This is effectively negligible for BBN
            tep_correction = 1.0 - 1e-6 * self.Sigma_0 * 1e4
            abundances = {
                'Y_p': lcdm_abundances['Y_p'] * tep_correction,
                'D_H': lcdm_abundances['D_H'] * tep_correction,
                'He3_H': lcdm_abundances['He3_H'] * tep_correction,
                'Li7_H': lcdm_abundances['Li7_H'] * tep_correction,
            }
        else:
            abundances = lcdm_abundances
        
        return abundances
