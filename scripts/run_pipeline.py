#!/usr/bin/env python3
"""TEP-C0 canonical pipeline runner.

The runner is intentionally strict: failed dependencies are recorded as failed
results, not silently ignored.
"""

from __future__ import annotations

import argparse
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, str(Path(__file__).parent / "steps"))

from c0_common import step_json_path

PIPELINE_STEPS = [
    ("step_000_data_download", "Download Pantheon+ and auxiliary data", []),
    ("step_000_data_ingestion", "Ingest and validate data", ["step_000_data_download"]),
    ("step_001_transport_kernel", "Define TEP transport kernel", []),
    ("step_022_three_model_comparison", "Fit M0/LCDM, M1/Mixed, M2/TEP", ["step_000_data_ingestion", "step_001_transport_kernel"]),
    ("step_018_tep_mcmc_inference", "Independent posterior check for TEP parameters", ["step_022_three_model_comparison"]),
    ("step_023_sn_time_dilation_test", "SN time dilation test", ["step_022_three_model_comparison"]),
    ("step_024_tolman_surface_brightness", "Tolman surface brightness", ["step_022_three_model_comparison"]),
    ("step_025_distance_duality_test", "Distance-duality relation", ["step_022_three_model_comparison"]),
    ("step_026_redshift_drift_forecast", "Redshift drift predictions", ["step_022_three_model_comparison"]),
    ("step_027_environment_residuals", "Environment correlation", ["step_022_three_model_comparison"]),
    ("step_010_cmb_blackbody_preservation", "CMB blackbody check", []),
    ("step_011_bao_acoustic_projection", "BAO constraints", ["step_022_three_model_comparison"]),
    ("step_017_tep_boltzmann_solver", "CMB spectra diagnostic", ["step_022_three_model_comparison"]),
    ("step_014_cmb_acoustic_projection", "CMB acoustic peaks", ["step_017_tep_boltzmann_solver"]),
    ("step_028_cmb_full_spectra", "CMB full-spectra gate", ["step_017_tep_boltzmann_solver"]),
    ("step_012_bbn_preservation_registry", "BBN registry check", ["step_022_three_model_comparison"]),
    ("step_029_bbn_preservation", "BBN abundance diagnostic", ["step_022_three_model_comparison"]),
    ("step_015_structure_growth_solver", "Structure growth", ["step_022_three_model_comparison"]),
    ("step_016_global_likelihood_synthesis", "Global likelihood", ["step_022_three_model_comparison", "step_018_tep_mcmc_inference"]),
    ("step_030_model_comparison_statistics", "Model-comparison statistics gate", ["step_022_three_model_comparison"]),
    ("step_031_sensitivity_analysis", "Sensitivity analysis", ["step_022_three_model_comparison"]),
    ("step_013_explanatory_evidence_matrix", "Evidence matrix", ["step_022_three_model_comparison", "step_018_tep_mcmc_inference"]),
    ("step_007_level3_gate_registry", "Evidence gate registry", ["step_013_explanatory_evidence_matrix"]),
    ("step_008_claim_consistency_audit", "Claim audit", ["step_022_three_model_comparison", "step_013_explanatory_evidence_matrix"]),
    ("step_009_evidence_gate_summary", "Final summary", ["step_007_level3_gate_registry", "step_008_claim_consistency_audit"]),
    ("step_020_tep_blind_injection_recovery", "Blind injection recovery", ["step_022_three_model_comparison"]),
    ("step_021_expansion_falsifier", "Expansion discriminator", ["step_025_distance_duality_test", "step_024_tolman_surface_brightness"]),
    ("step_032_full_physics_implementation", "Full physics implementation audit", ["step_022_three_model_comparison"]),
    ("step_033_cobaya_tep_inference", "Cobaya TEP-CLASS inference", ["step_022_three_model_comparison"]),
]


def run_step(step_module: str) -> dict:
    try:
        module = __import__(step_module)
        if not hasattr(module, "run"):
            return {"step": step_module, "error": "No run() function"}
        return module.run()
    except Exception as exc:
        return {"step": step_module, "error": str(exc), "traceback": traceback.format_exc()}


def selected_steps(specific_steps: list[str] | None) -> list[tuple[str, str, list[str]]]:
    if specific_steps is None:
        return PIPELINE_STEPS
    requested = set(specific_steps)
    step_deps = {name: deps for name, _, deps in PIPELINE_STEPS}
    known = set(step_deps)
    unknown = sorted(requested - known)
    if unknown:
        raise SystemExit(f"Unknown step(s): {', '.join(unknown)}")

    expanded = set(requested)
    changed = True
    while changed:
        changed = False
        for name in list(expanded):
            for dep in step_deps[name]:
                if dep not in expanded:
                    expanded.add(dep)
                    changed = True
    return [step for step in PIPELINE_STEPS if step[0] in expanded]


def run_pipeline(specific_steps: list[str] | None = None, resume: bool = False) -> dict:
    print("=" * 70)
    print("TEP-C0 COSMOLOGICAL PIPELINE (STRICT EMPIRICAL MODE)")
    print("=" * 70)
    print()

    steps_to_run = selected_steps(specific_steps)
    results: dict[str, dict] = {}
    completed: set[str] = set()
    failed: list[str] = []
    skipped: list[str] = []

    for step_module, description, deps in steps_to_run:
        if resume and step_json_path(step_module).exists():
            print(f"[SKIP] {step_module}: existing result found")
            results[step_module] = {
                "step": step_module,
                "status": "skipped_existing_result",
                "result_path": str(step_json_path(step_module)),
            }
            completed.add(step_module)
            continue

        missing = [dep for dep in deps if dep not in completed]
        if missing:
            print(f"[SKIP] {step_module}: Missing dependencies {missing}")
            results[step_module] = {
                "step": step_module,
                "error": "missing dependencies",
                "missing_dependencies": missing,
            }
            skipped.append(step_module)
            continue

        print(f"[RUN] {step_module}: {description}")
        start = time.time()
        result = run_step(step_module)
        elapsed = time.time() - start
        results[step_module] = result

        if "error" in result:
            print(f"  [FAIL] {elapsed:.1f}s: {result['error']}")
            failed.append(step_module)
        else:
            print(f"  [OK] {elapsed:.1f}s")
            completed.add(step_module)

    print()
    print("=" * 70)
    print(f"PIPELINE COMPLETE: {len(completed)}/{len(steps_to_run)} steps succeeded")
    if failed:
        print(f"FAILED: {', '.join(failed)}")
    if skipped:
        print(f"SKIPPED: {', '.join(skipped)}")
    print("=" * 70)
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description="TEP-C0 Pipeline")
    parser.add_argument("--steps", nargs="+", help="Specific step module names to run")
    parser.add_argument("--core", action="store_true", help="Run the core statistical steps")
    parser.add_argument("--resume", action="store_true", help="Skip steps whose result JSON already exists")
    args = parser.parse_args()

    if args.core:
        steps = [
            "step_022_three_model_comparison",
            "step_018_tep_mcmc_inference",
            "step_030_model_comparison_statistics",
        ]
    else:
        steps = args.steps

    results = run_pipeline(steps, resume=args.resume)
    return 1 if any("error" in result for result in results.values()) else 0


if __name__ == "__main__":
    sys.exit(main())
