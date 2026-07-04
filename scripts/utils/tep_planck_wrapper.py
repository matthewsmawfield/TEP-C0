"""Wrapped Planck likelihoods for TEP-C0 Cobaya runs.

Catches NaN/inf from the Planck clipy likelihood when TEP-CLASS produces
unphysical spectra, returning -inf so Cobaya rejects the proposal cleanly.
"""

from __future__ import annotations

import numpy as np
from typing import Dict, Any

from cobaya.likelihood import Likelihood


class PlanckTTTEEE_Wrapper(Likelihood):
    """Wrapper for Planck 2018 high-l TTTEEE that guards against TEP NaN/inf."""

    def initialize(self):
        from pathlib import Path
        from cobaya.likelihoods.planck_2018_highl_plik import TTTEEE
        # Auto-detect clik file if not provided
        clik_file = getattr(self, 'clik_file', None)
        if clik_file is None:
            search_paths = [
                Path(self.path) / "baseline" / "plc_3.0" / "hi_l" / "plik" / "plik_rd12_HM_v22b_TTTEEE.clik",
                Path(self.path) / "hi_l" / "plik" / "plik_rd12_HM_v22b_TTTEEE.clik",
            ]
            for p in search_paths:
                if p.exists():
                    clik_file = str(p.relative_to(self.path))
                    break
        self._inner = TTTEEE({"path": self.path, "clik_file": clik_file})
        self._inner.initialize()
        self._nan_count = 0

    def get_requirements(self):
        return self._inner.get_requirements()

    def logp(self, **params_values) -> float:
        try:
            lkl = self._inner.logp(**params_values)
        except Exception as e:
            self._nan_count += 1
            return -np.inf
        if lkl is None or not np.isfinite(lkl):
            self._nan_count += 1
            return -np.inf
        return float(lkl)


class PlanckLowT_Wrapper(Likelihood):
    """Wrapper for Planck 2018 low-l TT that guards against TEP NaN/inf."""

    def initialize(self):
        from cobaya.likelihoods.planck_2018_lowl import TT
        self._inner = TT({"path": self.path})
        self._inner.initialize()
        self._nan_count = 0

    def get_requirements(self):
        return self._inner.get_requirements()

    def logp(self, **params_values) -> float:
        try:
            lkl = self._inner.logp(**params_values)
        except Exception:
            self._nan_count += 1
            return -np.inf
        if lkl is None or not np.isfinite(lkl):
            self._nan_count += 1
            return -np.inf
        return float(lkl)


class PlanckLowE_Wrapper(Likelihood):
    """Wrapper for Planck 2018 low-l EE that guards against TEP NaN/inf."""

    def initialize(self):
        from cobaya.likelihoods.planck_2018_lowl import EE
        self._inner = EE({"path": self.path})
        self._inner.initialize()
        self._nan_count = 0

    def get_requirements(self):
        return self._inner.get_requirements()

    def logp(self, **params_values) -> float:
        try:
            lkl = self._inner.logp(**params_values)
        except Exception:
            self._nan_count += 1
            return -np.inf
        if lkl is None or not np.isfinite(lkl):
            self._nan_count += 1
            return -np.inf
        return float(lkl)
