# TEP-C0 Repository Organization Audit

**Date**: 2026-05-02
**Auditor**: Cascade

## Summary

The repository has accumulated development artifacts that need organization for professional presentation.

## Issues Found

### 1. Root Directory Clutter
**Location**: `/`

**Problem**: Development notes and temporary files at root level.

| File | Size | Issue |
|------|------|-------|
| `feedback` | 14.7 KB | Development notes, not documentation |
| `plan` | 33.7 KB | Development planning, superseded by ROADMAP |
| `WEEK2_ACOUSTIC_ANALYSIS.md` | - | Should be in docs/development/ |
| `WEEK3_STATUS.md` | - | Should be in docs/development/ |
| `WEEK4_PROGRESS.md` | - | Should be in docs/development/ |
| `WEEK4_JOINT_SUCCESS.md` | - | Should be in docs/development/ |

### 2. External Directory Cleanup
**Location**: `external/`

| Item | Issue |
|------|-------|
| `class_backup_20260502_063702/` | Temporary backup, should be removed |

### 3. Chain Files (Development Artifacts)
**Location**: `chains/`

| Pattern | Count | Issue |
|---------|-------|-------|
| `tep_*.txt` | 2 | Test chains from development |
| `tep_*.yaml` | 4 | Cobaya configs from testing |
| `tep_*.progress` | 2 | Progress files |

**Solution**: Move to `archive/cobaya_testing/`

### 4. Documentation Organization
**Location**: Root level MD files

Current state: 13 markdown files at root
Recommended: Only README, CHANGELOG, LICENSE at root

### 5. Scripts Directory
**Status**: ✅ Well organized
- `scripts/steps/` follows naming convention
- `scripts/steps/core/` has supporting modules
- `c0_common.py` provides shared utilities

### 6. Results Directory
**Status**: ✅ Well organized
- `results/outputs/` - JSON and CSV outputs
- `results/figures/` - Plots
- `results/tables/` - LaTeX tables

## Cleanup Plan

### Phase 1: Archive Development Files
```
feedback → archive/development_notes/feedback.md
plan → archive/development_notes/plan.md
WEEK*.md → docs/development/progress/
```

### Phase 2: Remove Temporary Files
```
rm -rf external/class_backup_20260502_063702/
mv chains/tep_* archive/cobaya_testing/
```

### Phase 3: Consolidate Documentation
```
Root MD files → Organize by type:
- docs/development/ - Progress reports
- docs/status/ - Gate assessments
- docs/roadmap/ - Planning documents
```

## Directory Structure (Target)

```
TEP-C0/
├── README.md                    # Main entry point
├── LICENSE                      # License
├── CHANGELOG.md                 # Version history
├── VERSION.json               # Current version
│
├── archive/                     # Historical artifacts
│   ├── cobaya_testing/          # Week 4 test chains
│   ├── development_notes/       # feedback, plan
│   └── week4_standalone/        # Already organized
│
├── docs/                        # Documentation
│   ├── development/             # Progress reports
│   │   ├── WEEK2_ACOUSTIC_ANALYSIS.md
│   │   ├── WEEK3_STATUS.md
│   │   ├── WEEK4_JOINT_SUCCESS.md
│   │   └── WEEK4_PROGRESS.md
│   ├── status/                    # Current status
│   │   ├── GATE_STATUS_REALISTIC.txt
│   │   ├── PARTIAL_GATE_ASSESSMENT.md
│   │   └── RESEARCH_GRADE_STATUS.md
│   └── roadmap/                   # Planning
│       └── ROADMAP_FULL_IMPLEMENTATION.md
│
├── external/                    # External dependencies
│   ├── class/                     # CLASS (git submodule)
│   └── class_tep_mod/           # TEP patch source
│
├── scripts/                     # Pipeline code
│   └── steps/                     # Analysis steps
│       ├── core/                  # Shared modules
│       └── step_*.py            # Individual steps
│
├── results/                     # Pipeline outputs
│   ├── outputs/                   # JSON/CSV data
│   ├── figures/                   # Plots
│   └── tables/                    # LaTeX tables
│
├── data/                        # Input data
│   └── raw/                       # Pantheon+, etc.
│
└── tests/                       # Test suite
    ├── integration/
    └── unit/
```

## Post-Cleanup Root Files

Only these files should remain at root:
- `README.md`
- `LICENSE`
- `CHANGELOG.md`
- `VERSION.json`
- `.gitignore`
- `CITATION.cff`
- `requirements*.txt`

All other MD/TXT files move to appropriate `docs/` subdirectories.
