"""CMB CLASS resolver stub.

Provides a minimal fallback when CLASS/TEP-CLASS is not available.
The step degrades gracefully to blocked status with research-grade=False.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Any


@dataclass
class CMBRun:
    """Container for CMB computation results."""

    lcdm_reference: dict[str, Any] | None = None
    tep_zero_limit: dict[str, Any] | None = None
    tep_spectra: dict[str, Any] | None = None
    validation: dict[str, Any] = None

    def __post_init__(self):
        if self.validation is None:
            self.validation = {
                "research_grade_cmb": False,
                "status": "blocked",
                "reason": "TEP-CLASS v2 not available in active build",
            }


def resolve_cmb(
    epsilon_t: float = 0.0,
    z_t: float = 5.0,
    n_t: float = 1.0,
    lmax: int = 2500,
) -> CMBRun:
    """Resolve CMB spectra for TEP cosmology.

    When TEP-CLASS is unavailable, returns a blocked CMBRun so the
    pipeline step can still complete with informative status.
    """
    return CMBRun()
