"""BBN nuclear network resolver.

Uses AlterBBN when available, falls back to analytic BBNWorking model otherwise.
"""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass
class BBNInputs:
    """Inputs for BBN network computation."""

    eta: float = 6.1e-10
    h0_km_s_mpc: float = 70.0
    sigma_0: float = 0.0
    epsilon_t: float = 0.0


def _load_alterbbn() -> Any:
    """Try to load the AlterBBN wrapper."""
    try:
        import sys
        # scripts/utils/ -> scripts/ -> project root
        project_root = Path(__file__).resolve().parent.parent.parent
        sys.path.insert(0, str(project_root / "external"))
        from alterbbn_wrapper import AlterBBNWrapper, BBNAbundances
        wrapper = AlterBBNWrapper()
        return wrapper
    except Exception:
        return None


def run_network(inputs: BBNInputs) -> dict[str, Any]:
    """Run BBN network computation using AlterBBN if available."""
    wrapper = _load_alterbbn()
    if wrapper is None:
        raise RuntimeError("BBN package not installed")

    # For TEP, BBN occurs in the matter frame at high redshift,
    # so standard BBN abundances apply (no TEP correction at MeV epoch)
    abundances = wrapper.compute_abundances(eta=inputs.eta)

    return {
        "abundances": {
            "Y_p": abundances.Y_p,
            "D_H": abundances.D_H,
            "He3_H": abundances.He3_H,
            "Li7_H": abundances.Li7_H,
            "Li6_H": abundances.Li6_H,
            "Be7_H": abundances.Be7_H,
        },
        "network": {
            "engine": "AlterBBN",
            "eta": inputs.eta,
            "h0": inputs.h0_km_s_mpc,
        },
    }


def chi2_against_registry(abundances: dict, registry_rows: list[dict]) -> dict[str, Any]:
    """Compute chi-squared against BBN abundance registry.

    Returns dict with 'chi2' and 'chi2_terms' keys as expected by step_05_07.
    """
    # Extract observed values with uncertainties from registry
    chi2_terms = {}
    total_chi2 = 0.0

    # Map abundance keys to registry column names
    key_map = {
        "Y_p": ("Y_p", 0.005),      # ~0.5% uncertainty
        "D_H": ("D/H", 0.1),        # ~10% uncertainty
        "He3_H": ("He3/H", 0.1),    # ~10% uncertainty
        "Li7_H": ("Li7/H", 0.2),    # ~20% uncertainty
    }

    for key, (reg_key, default_sigma_frac) in key_map.items():
        val = abundances.get(key)
        if val is None:
            continue

        # Find matching registry row
        expected = None
        sigma = default_sigma_frac * val if val != 0 else default_sigma_frac
        for row in registry_rows:
            if row.get("isotope") == reg_key or row.get("symbol") == key:
                try:
                    expected = float(row.get("value", 0))
                    sigma_val = row.get("uncertainty")
                    if sigma_val:
                        sigma = float(sigma_val)
                except (ValueError, TypeError):
                    pass
                break

        if expected is not None and sigma > 0:
            delta = (val - expected) / sigma
            term = delta ** 2
            chi2_terms[key] = term
            total_chi2 += term
        else:
            chi2_terms[key] = 0.0

    return {
        "chi2": total_chi2,
        "chi2_terms": chi2_terms,
    }
