"""BBN working model stub.

Provides a minimal analytic BBN solver as fallback when nuclear network
packages are unavailable.
"""

from __future__ import annotations
from typing import Any


class BBNWorking:
    """Minimal BBN working model for light-element abundances."""

    def __init__(self, eta: float = 6.1e-10, Sigma_0: float = 0.0, H0: float = 70.0):
        self.eta = eta
        self.Sigma_0 = Sigma_0
        self.H0 = H0

    def solve(self) -> dict[str, float]:
        """Return standard BBN abundances (no TEP correction in stub).

        Values from standard BBN (Coc et al. 2014 / Planck 2018):
        """
        return {
            "Y_p": 0.24709,      # primordial helium mass fraction
            "D_H": 2.6e-5,       # deuterium to hydrogen ratio
            "He3_H": 1.0e-5,     # helium-3 to hydrogen ratio
            "Li7_H": 4.6e-10,    # lithium-7 to hydrogen ratio
            "n_p_ratio": 1.0e-10,
        }
