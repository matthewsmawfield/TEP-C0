#!/usr/bin/env python3
"""TEP-C0 canonical pipeline runner.

The runner is intentionally strict: failed dependencies are recorded as failed
results, not silently ignored.

Manuscript-to-code mapping for external auditors:
    rho_half (manuscript, Section 2.5)  ->  core.cosmology.TEPCosmology.RHO_HALF
    Screening formula S(rho)            ->  TEPCosmology.screening_function(rho)
    Value: 0.5 M_sun / pc^3
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
    ("step_01_01_data_download", "Data download", []),
    ("step_01_02_data_ingestion", "Data ingestion", ["step_01_01_data_download"]),
    ("step_01_03_download_ddr", "Download ddr", []),
    ("step_01_04_download_sb", "Download sb", []),
    ("step_01_05_download_sz", "Download sz", []),
    ("step_01_06_download_sgl", "Download sgl", []),
    ("step_01_07_download_desi", "Download desi", []),
    ("step_01_08_compile_sb", "Compile sb", ["step_01_04_download_sb"]),
    ("step_02_01_transport_kernel", "Transport kernel", []),
    ("step_02_02_theory_derivation", "Theory derivation", []),
    ("step_02_03_physics_implementation", "Physics implementation", []),
    ("step_03_01_three_model_comparison", "Three model comparison", ["step_01_02_data_ingestion", "step_02_01_transport_kernel"]),
    ("step_03_02_independent_mcmc", "Independent mcmc", ["step_03_01_three_model_comparison"]),
    ("step_03_04_cobaya_mcmc", "Cobaya mcmc", ["step_03_01_three_model_comparison"]),
    ("step_03_05_analyze_cobaya", "Analyze cobaya", ["step_03_04_cobaya_mcmc"]),
    ("step_03_06_cobaya_verbose", "Cobaya verbose", ["step_03_01_three_model_comparison"]),
    ("step_03_07_likelihood_synthesis", "Likelihood synthesis", ["step_03_01_three_model_comparison", "step_03_02_independent_mcmc", "step_03_04_cobaya_mcmc"]),
    ("step_04_01_sn_time_dilation", "Sn time dilation", ["step_03_01_three_model_comparison"]),
    ("step_04_02_sn_tolman", "Sn tolman", ["step_03_01_three_model_comparison"]),
    ("step_04_03_tolman_sb", "Tolman sb", ["step_03_01_three_model_comparison"]),
    ("step_04_04_distance_duality", "Distance duality", ["step_03_01_three_model_comparison", "step_01_03_download_ddr"]),
    ("step_04_05_ddr_threeway", "Ddr threeway", ["step_04_04_distance_duality"]),
    ("step_04_06_screening_fit", "Screening fit", ["step_03_01_three_model_comparison"]),
    ("step_04_07_highz_ddr", "Highz ddr", ["step_03_01_three_model_comparison", "step_01_07_download_desi"]),
    ("step_05_01_cmb_blackbody", "Cmb blackbody", []),
    ("step_05_03_cmb_boltzmann", "Cmb boltzmann", ["step_03_01_three_model_comparison"]),
    ("step_05_04_cmb_spectra", "Cmb spectra", ["step_05_03_cmb_boltzmann"]),
    ("step_05_05_cmb_consistency", "Cmb consistency", ["step_05_03_cmb_boltzmann"]),
    ("step_05_06_bbn_registry", "Bbn registry", ["step_03_01_three_model_comparison"]),
    ("step_05_07_bbn_preservation", "Bbn preservation", ["step_03_01_three_model_comparison"]),
    ("step_05_08_cmb_acoustic", "Cmb acoustic", ["step_05_03_cmb_boltzmann"]),
    ("step_05_09_jordan_frame_proof", "Jordan frame proof", []),
    ("step_06_01_bao_projection", "Bao projection", ["step_03_01_three_model_comparison"]),
    ("step_06_02_bao_likelihood", "Bao likelihood", ["step_03_01_three_model_comparison"]),
    ("step_06_03_growth_solver", "Growth solver", ["step_03_01_three_model_comparison"]),
    ("step_06_04_growth_validation", "Growth validation", ["step_06_03_growth_solver"]),
    ("step_06_05_growth_rsd", "Growth rsd", ["step_03_01_three_model_comparison"]),
    ("step_07_01_mixed_forecast", "Mixed forecast", ["step_03_01_three_model_comparison"]),
    ("step_07_02_redshift_drift", "Redshift drift", ["step_03_01_three_model_comparison"]),
    ("step_07_03_jwst_test", "Jwst test", ["step_03_01_three_model_comparison"]),
    ("step_07_04_gw_sirens", "Gw sirens", []),
    ("step_07_05_weak_lensing_plan", "Weak lensing plan", []),
    ("step_07_06_weak_lensing", "Weak lensing", []),
    ("step_07_07_blind_injection", "Blind injection", ["step_03_01_three_model_comparison"]),
    ("step_08_01_expansion_falsifier", "Expansion falsifier", ["step_04_04_distance_duality", "step_04_03_tolman_sb"]),
    ("step_08_02_comparison_stats", "Comparison stats", ["step_03_01_three_model_comparison"]),
    ("step_08_03_sensitivity_analysis", "Sensitivity analysis", ["step_03_01_three_model_comparison"]),
    ("step_08_04_evidence_matrix", "Evidence matrix", ["step_03_01_three_model_comparison", "step_03_02_independent_mcmc", "step_03_04_cobaya_mcmc"]),
    ("step_08_05_gate_registry", "Gate registry", ["step_08_04_evidence_matrix"]),
    ("step_08_06_claim_audit", "Claim audit", ["step_03_01_three_model_comparison", "step_08_04_evidence_matrix"]),
    ("step_08_07_final_summary", "Final summary", ["step_08_05_gate_registry", "step_08_06_claim_audit"]),
    ("step_08_08_diagnostic_plots", "Diagnostic plots", ["step_03_01_three_model_comparison", "step_04_04_distance_duality"]),
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

        if result is None:
            result = {"step": step_module, "status": "failed", "error": "Step returned None (missing return value)"}
        results[step_module] = result

        if isinstance(result, dict) and "error" in result:
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
            "step_03_01_three_model_comparison",
            "step_03_02_independent_mcmc",
            "step_08_02_comparison_stats",
        ]
    else:
        steps = args.steps

    results = run_pipeline(steps, resume=args.resume)
    return 1 if any("error" in result for result in results.values()) else 0


if __name__ == "__main__":
    sys.exit(main())
