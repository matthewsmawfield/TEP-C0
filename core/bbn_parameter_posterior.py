"""BBN parameter posterior MCMC stub.

Provides a minimal fallback for BBN parameter inference.
"""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import csv


@dataclass
class BBNPosteriorConfig:
    """Configuration for BBN parameter MCMC."""

    sigma_0: float = 0.0
    epsilon_t: float = 0.0
    h0_km_s_mpc: float = 70.0
    n_walkers: int = 16
    n_steps: int = 500
    burn_in: int = 100


def load_bbn_registry() -> list[dict]:
    """Load BBN abundance registry from processed data."""
    from pathlib import Path
    path = Path("data/processed/tep_c0_bbn_abundance_registry.csv")
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def run_bbn_mcmc(
    registry_rows: list[dict],
    config: BBNPosteriorConfig,
    seed: int = 42,
) -> dict[str, Any]:
    """Run BBN parameter posterior MCMC.

    Fallback: returns empty posterior with informative status.
    """
    return {
        "status": "skipped",
        "reason": "BBN package not installed; no posterior available",
        "parameters": {
            "sigma_0": config.sigma_0,
            "epsilon_t": config.epsilon_t,
            "h0": config.h0_km_s_mpc,
        },
    }
