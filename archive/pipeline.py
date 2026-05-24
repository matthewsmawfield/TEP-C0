#!/usr/bin/env python3
"""Compatibility entrypoint for the canonical TEP-C0 pipeline."""

from __future__ import annotations

import argparse
import sys

from run_pipeline import run_pipeline


def main() -> int:
    parser = argparse.ArgumentParser(description="TEP-C0 Pipeline")
    parser.add_argument("--steps", nargs="+", help="Specific step module names to run")
    parser.add_argument("--core", action="store_true", help="Run core analysis only")
    args = parser.parse_args()

    if args.core:
        steps = [
            "step_022_three_model_comparison",
            "step_018_tep_mcmc_inference",
            "step_030_model_comparison_statistics",
            "step_023_sn_time_dilation_test",
            "step_024_tolman_surface_brightness",
            "step_025_distance_duality_test",
            "step_026_redshift_drift_forecast",
            "step_027_environment_residuals",
        ]
        results = run_pipeline(steps)
    else:
        results = run_pipeline(args.steps)

    return 1 if any("error" in result for result in results.values()) else 0


if __name__ == "__main__":
    sys.exit(main())
