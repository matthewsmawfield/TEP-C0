#!/usr/bin/env python3
"""Step 030: model-comparison statistics gate.

This step computes information criteria from the step 022 likelihoods and then
checks whether publication-grade Bayesian evidence is actually present. It does
not promote AIC/BIC into a claim of proof.
"""

from __future__ import annotations

from math import log
from pathlib import Path

from c0_common import (
    TEPLogger,
    ensure_dirs,
    print_status,
    read_json,
    rounded,
    set_step_logger,
    step_csv_path,
    step_json_path,
    write_csv,
    write_json,
)

STEP_ID = "step_08_02_comparison_stats"
SOURCE_STEP = "step_03_01_three_model_comparison"


def _parameter_count(model_payload: dict) -> int:
    return len(model_payload.get("parameters_mle", {}))


def run() -> dict:
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    source_path = step_json_path(SOURCE_STEP)
    if not source_path.exists():
        raise FileNotFoundError(f"{SOURCE_STEP} results not found: {source_path}")

    source = read_json(source_path)
    n_data = source.get("data", {}).get("n_supernovae", source.get("n_supernovae"))
    if not n_data:
        raise RuntimeError(f"{SOURCE_STEP} did not report data.n_supernovae")

    rows = []
    evidence_available = True
    for model_id, model_payload in source.get("models", {}).items():
        log_likelihood = model_payload.get("log_likelihood_mle", model_payload.get("log_likelihood"))
        if log_likelihood is None:
            raise RuntimeError(f"{model_id} is missing log_likelihood")
        k = _parameter_count(model_payload)
        aic = 2 * k - 2 * log_likelihood
        bic = k * log(n_data) - 2 * log_likelihood
        log_evidence = model_payload.get("log_evidence")
        log_evidence_err = model_payload.get("log_evidence_error")
        evidence_available = evidence_available and log_evidence is not None
        rows.append({
            "model": model_id,
            "n_parameters": k,
            "log_likelihood": rounded(log_likelihood, 6),
            "aic": rounded(aic, 6),
            "bic": rounded(bic, 6),
            "log_evidence": "" if log_evidence is None else rounded(log_evidence, 6),
            "log_evidence_error": "" if log_evidence_err is None else rounded(log_evidence_err, 6),
        })

    best_aic = min(row["aic"] for row in rows)
    best_bic = min(row["bic"] for row in rows)
    for row in rows:
        row["delta_aic"] = rounded(row["aic"] - best_aic, 6)
        row["delta_bic"] = rounded(row["bic"] - best_bic, 6)

    csv_path = step_csv_path(STEP_ID)
    write_csv(csv_path, rows)

    source_validation = source.get("validation", {})
    # research_grade_model_comparison flag. A model comparison is research grade
    # if it uses real data, maintains strict provenance, and achieves MCMC convergence,
    # regardless of whether the TEP model 'wins' or 'loses'.
    source_research_grade = (
        source_validation.get("research_grade_provenance") == True
        and source_validation.get("mcmc_converged") == True
        and evidence_available
    )
    blockers = [] if source_research_grade else source_validation.get("blockers", [
        "Upstream step_03_01 must achieve MCMC convergence and strict provenance.",
    ])

    payload = {
        "step": STEP_ID,
        "source_step": SOURCE_STEP,
        "model_statistics": rows,
        "validation": {
            "aic_bic_available": True,
            "nested_sampling_evidence_available": evidence_available,
            "research_grade_model_comparison": source_research_grade,
            "claim_gate": "open" if source_research_grade else "blocked",
            "blockers": blockers if evidence_available else [
                "Step 022 must report log_evidence and log_evidence_error for every model.",
                "AIC/BIC are useful diagnostics but are not a replacement for the stated nested-sampling evidence.",
            ],
        },
        "artifacts": {"csv": str(csv_path)},
    }
    write_json(step_json_path(STEP_ID), payload)
    print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
    return payload


if __name__ == "__main__":
    print(run())
