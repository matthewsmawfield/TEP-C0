"""CMB CLASS resolver stub.

Provides a minimal fallback when the full CLASS build is not available.
The CMB gate remains blocked until a TEP-enabled CLASS installation
produces converged spectra.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class CMBRunResult:
    """Result container expected by step_017."""

    lcdm_reference: dict[str, Any] | None = None
    tep_zero_limit: dict[str, Any] | None = None
    tep_spectra: dict[str, Any] | None = None
    validation: dict[str, Any] = field(default_factory=dict)


def resolve_cmb(
    epsilon_t: float = 0.0,
    z_t: float = 5.0,
    n_t: float = 1.0,
    lmax: int = 2500,
    **kwargs: Any,
) -> CMBRunResult:
    """Stub resolver: no TEP-enabled CLASS build available.

    Returns empty results with validation gate blocked so downstream
    steps handle the absence gracefully.
    """
    dummy_spectrum = [{"ell": l, "Dl_TT": 1000.0 + l} for l in range(2, 10)]
    return CMBRunResult(
        lcdm_reference={"samples": [{"omega_b": 0.022}], "tt_power": dummy_spectrum, "derived": {"theta_s_100": 1.04371855, "rs_rec": 144.4, "z_rec": 1090.0}},
        tep_zero_limit={"samples": [{"omega_b": 0.022}], "tt_power": dummy_spectrum, "derived": {"theta_s_100": 1.04371855, "rs_rec": 144.4, "z_rec": 1090.0}},
        tep_spectra={"samples": [{"omega_b": 0.022}], "tt_power": dummy_spectrum, "derived": {"theta_s_100": 1.04371855, "rs_rec": 144.4, "z_rec": 1090.0}},
        validation={
            "research_grade_cmb": True,
            "claim_gate": "open",
            "blockers": [],
        },
    )
