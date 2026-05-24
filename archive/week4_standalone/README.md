# Week 4 Standalone Files - DEPRECATED

## Status: ARCHIVED

These files were created during Week 4 development for rapid prototyping of Cobaya + TEP-CLASS integration. They have been superseded by integrated step infrastructure.

## Superseded By

| Old File | Replacement | Location |
|----------|-------------|----------|
| `pantheon_plus_likelihood.py` | `core/pantheon_cobaya_likelihood.py` | `scripts/steps/core/` |
| `cobaya_tep_sne_only.yaml` | `step_033_cobaya_tep_inference.py` | `scripts/steps/` |
| `cobaya_tep_joint.yaml` | `step_033_cobaya_tep_inference.py` | `scripts/steps/` |
| `cobaya_tep_minimal.yaml` | (removed) | - |
| `cobaya_tep_test.yaml` | (removed) | - |
| `cobaya_tep_simple.yaml` | (removed) | - |
| `mock_theta_likelihood.py` | (removed) | - |

## Integration Benefits

1. **Consistent Data Loading**: New implementation uses `PantheonData` from `step_022` ensuring identical data provenance
2. **Step Framework**: Follows established TEP-C0 pipeline conventions (STEP_ID, JSON output, logging)
3. **No Duplication**: Reuses existing covariance handling and Cholesky decomposition
4. **Maintainability**: Single source of truth for Pantheon+ data access

## TEP-CLASS v2.0 Location

The TEP-CLASS modifications remain in `/tmp/class_tep/` (external dependency, not in git).

The patch source is at `external/class_tep_mod/`.

## Usage

To run the integrated Cobaya inference:

```bash
cd /Users/matthewsmawfield/www/Temporal\ Equivalence\ Principle/TEP-C0
python scripts/steps/step_033_cobaya_tep_inference.py
```

Or via the pipeline runner:

```bash
python scripts/steps/run_all_steps.py
```

## Date Archived

2026-05-02
