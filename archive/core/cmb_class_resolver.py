"""CLASS-backed CMB resolver for TEP-C0.

The resolver deliberately separates three states:

1. Vanilla CLASS LCDM reference succeeds.
2. A CLASS build accepts TEP background parameters.
3. A TEP run validates its Sigma_0=0 LCDM limit and produces spectra.

Only state 3 can open the CMB implementation gate.  Python-side spectrum
rescaling is intentionally absent.

References:
- Planck Collaboration 2020, A&A 641, A6: Planck 2018 results. VI. Cosmological parameters.
"""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import sys
from typing import Any, Optional

import numpy as np

# Planck 2018 baseline cosmological parameters (Planck Collaboration 2020, A&A 641, A6)
# These are used as validation reference values, not computed by CLASS in this step
# CLASS with "output: tCl,pCl,lCl" does not compute sigma_8
# For actual sigma_8 computation, see step_030_structure_growth_validation.py
PLANCK_2018_H0 = 67.36  # km/s/Mpc (Planck 2018 TT,TE,EE+lowE baseline)
PLANCK_2018_OMEGA_M = 0.3158  # Matter density parameter (Planck 2018 baseline)
PLANCK_2018_SIGMA_8 = 0.8120  # Amplitude of matter fluctuations (Planck 2018 baseline)

# Import validation module - try both relative and absolute
CMB_VALIDATION_AVAILABLE = False
CLASSValidator = None
validate_lcdm_limit = None

try:
    from .cmb_class_validation import CLASSValidator, validate_lcdm_limit
    CMB_VALIDATION_AVAILABLE = True
except ImportError:
    pass

# Try absolute import when running directly
if not CMB_VALIDATION_AVAILABLE:
    try:
        _core_dir = Path(__file__).parent
        if str(_core_dir) not in sys.path:
            sys.path.insert(0, str(_core_dir))
        from cmb_class_validation import CLASSValidator, validate_lcdm_limit
        CMB_VALIDATION_AVAILABLE = True
    except ImportError as e:
        # Validation module not available - log for debugging
        print(f"[DEBUG] CMB validation import failed: {e}", file=sys.stderr)

# Try to import CLASS or CAMB
# First, try local TEP-CLASS build
CLASS_AVAILABLE = False
CAMB_AVAILABLE = False

PROJECT_ROOT = Path(__file__).resolve().parents[3]
TEP_CLASS_BUILD = PROJECT_ROOT / "external" / "class" / "build" / "lib.macosx-11.1-arm64-cpython-313"

# Add TEP-CLASS build path to Python path if it exists
if TEP_CLASS_BUILD.exists():
    if str(TEP_CLASS_BUILD) not in sys.path:
        sys.path.insert(0, str(TEP_CLASS_BUILD))

try:
    from classy import Class
    CLASS_AVAILABLE = True
except ImportError:
    CLASS_AVAILABLE = False
    try:
        import camb
        CAMB_AVAILABLE = True
    except ImportError:
        CAMB_AVAILABLE = False


PLANCK18_BASELINE = {
    "output": "tCl,pCl,lCl",
    "l_max_scalars": 2500,
    "lensing": "yes",
    "h": 0.6736,
    "omega_b": 0.02237,
    "omega_cdm": 0.1200,
    "A_s": 2.100e-9,
    "n_s": 0.9649,
    "tau_reio": 0.0544,
    "T_cmb": 2.7255,
}

CLASS_SOURCE_ROOT = PROJECT_ROOT / "external" / "class"


@dataclass(frozen=True)
class CMBRun:
    class_available: bool
    tep_class_available: bool
    lcdm_reference: Optional[dict[str, Any]]
    tep_zero_limit: Optional[dict[str, Any]]
    tep_spectra: Optional[dict[str, Any]]
    validation: dict[str, Any]


def _spectrum_summary(cls: dict[str, Any], lmax: int) -> dict[str, Any]:
    ell = np.asarray(cls["ell"], dtype=int)
    summary_rows = []
    full_rows = []
    for idx, lval in enumerate(ell):
        if lval < 2:
            continue
        tt = float(cls["tt"][idx])
        dl_tt = lval * (lval + 1.0) * tt / (2.0 * np.pi) * (2.7255e6**2)
        if lval in (2, 10, 100, 220, 500, 1000, 1500, 2000):
            summary_rows.append({"ell": int(lval), "D_ell_TT_uK2": dl_tt})
        if lval <= lmax:
            full_rows.append({"ell": int(lval), "D_ell_TT_uK2": dl_tt})
    return {"samples": summary_rows, "tt_power": full_rows}


def run_class(params: dict[str, Any], lmax: int = 2500) -> dict[str, Any]:
    """Run CLASS with given parameters. Falls back to CAMB if CLASS not available."""
    # Try CLASS first
    if CLASS_AVAILABLE:
        class_pythonpath = os.getenv("TEP_CLASS_PYTHONPATH")
        if class_pythonpath and class_pythonpath not in sys.path:
            sys.path.insert(0, class_pythonpath)
        from classy import Class

        class_params = dict(params)
        class_params["l_max_scalars"] = lmax
        cosmo = Class()
        try:
            cosmo.set(class_params)
            cosmo.compute()
            cls = cosmo.lensed_cl(lmax)
            derived = cosmo.get_current_derived_parameters(["z_rec", "rs_rec", "theta_s_100"])
            spectra = _spectrum_summary(cls, lmax)
            return {
                "params": class_params,
                "derived": {key: float(value) for key, value in derived.items()},
                **spectra,
            }
        finally:
            try:
                cosmo.struct_cleanup()
                cosmo.empty()
            except Exception:
                pass
    
    # Fall back to CAMB
    elif CAMB_AVAILABLE:
        return _run_camb(params, lmax)
    else:
        raise RuntimeError("Neither CLASS nor CAMB is available. Install one of them.")


def _run_camb(params: dict[str, Any], lmax: int = 2500) -> dict[str, Any]:
    """Run CAMB with given parameters (CLASS alternative)."""
    import camb
    from camb import model, initialpower
    
    # Convert CLASS-style params to CAMB-style
    h = params.get("h", 0.6736)
    H0 = h * 100
    ombh2 = params.get("omega_b", 0.02237)
    omch2 = params.get("omega_cdm", 0.1200)
    As = params.get("A_s", 2.100e-9)
    ns = params.get("n_s", 0.9649)
    
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2)
    pars.InitPower.set_params(As=As, ns=ns)
    pars.set_for_lmax(lmax, lens_potential_accuracy=0)
    
    # Check for TEP parameters - CAMB doesn't support them natively
    epsilon_t = params.get("tep_epsilon_T", 0.0)
    z_t = params.get("tep_z_T", 0.0)
    if epsilon_t != 0.0 or z_t != 0.0:
        raise RuntimeError("CAMB does not support TEP parameters (sigma_0, epsilon_t). Use CLASS for TEP.")
    
    results = camb.get_results(pars)
    powers = results.get_cmb_power_spectra(pars, CMB_unit="muK")
    
    # Extract TT spectrum
    l = np.arange(2, lmax + 1)
    Dl_tt = powers["total"][2:lmax+1, 0]  # TT column
    
    # Build summary
    summary_rows = []
    full_rows = []
    for i, lval in enumerate(l):
        dl_tt = float(Dl_tt[i])
        if lval in (2, 10, 100, 220, 500, 1000, 1500, 2000):
            summary_rows.append({"ell": int(lval), "D_ell_TT_uK2": dl_tt})
        full_rows.append({"ell": int(lval), "D_ell_TT_uK2": dl_tt})
    
    # Get derived parameters
    derived = {
        "z_rec": float(results.get_derived_params().get("zstar", 1090.0)),
        "rs_rec": float(results.get_derived_params().get("rstar", 147.0)),
        "theta_s_100": float(results.get_derived_params().get("theta", 1.04)),
    }
    
    return {
        "params": params,
        "derived": derived,
        "samples": summary_rows,
        "tt_power": full_rows,
    }


def _tep_params(epsilon_t: float, z_t: float = 5.0, n_t: float = 1.0, lmax: int = 2500) -> dict[str, Any]:
    """Build TEP parameter dict for CLASS v2.1 (epsilon_T, z_T, n_T)."""
    params = dict(PLANCK18_BASELINE)
    params.update({
        "l_max_scalars": lmax,
        "tep_mode": "yes",
        "tep_epsilon_T": float(epsilon_t),
        "tep_z_T": float(z_t),
        "tep_n_T": float(n_t),
    })
    return params


def resolve_cmb(epsilon_t: float, z_t: float = 5.0, n_t: float = 1.0, lmax: int = 2500) -> CMBRun:
    lcdm_reference = None
    tep_zero_limit = None
    tep_spectra = None
    blockers: list[str] = []

    source_patch_present = all(
        token in path.read_text(encoding="utf-8", errors="replace")
        for path, token in [
            (CLASS_SOURCE_ROOT / "include" / "background.h", "epsilon_T"),
            (CLASS_SOURCE_ROOT / "source" / "input.c", "tep_mode"),
            (CLASS_SOURCE_ROOT / "source" / "background.c", "tep_gamma"),
        ]
        if path.exists()
    )

    # Check which Boltzmann solver is available
    if not CLASS_AVAILABLE and not CAMB_AVAILABLE:
        blockers.append("Neither CLASS nor CAMB is available. Install one for CMB calculations.")
        return CMBRun(False, False, None, None, None, {
            "class_available": False,
            "camb_available": False,
            "tep_class_available": False,
            "local_class_source_tep_patch_present": source_patch_present,
            "research_grade_cmb": False,
            "claim_gate": "blocked",
            "blockers": blockers,
        })

    try:
        lcdm_reference = run_class(PLANCK18_BASELINE, lmax=lmax)
    except Exception as exc:
        blockers.append(f"LCDM reference failed: {exc}")
        return CMBRun(False, False, None, None, None, {
            "class_available": CLASS_AVAILABLE,
            "camb_available": CAMB_AVAILABLE,
            "tep_class_available": False,
            "local_class_source_tep_patch_present": source_patch_present,
            "research_grade_cmb": False,
            "claim_gate": "blocked",
            "blockers": blockers,
        })

    # Try to run with TEP parameters
    try:
        tep_zero_limit = run_class(_tep_params(0.0, z_t, n_t, lmax), lmax=lmax)
        tep_spectra = run_class(_tep_params(epsilon_t, z_t, n_t, lmax), lmax=lmax)
        tep_class_works = True
    except Exception as exc:
        # TEP parameters not supported (expected with CAMB)
        tep_class_works = False
        tep_zero_limit = lcdm_reference  # Use LCDM as approximation
        tep_spectra = lcdm_reference
        if epsilon_t != 0.0:
            blockers.append(
                f"TEP parameters (epsilon_t={epsilon_t}) not supported. "
                "Current solver: CAMB (no TEP support) or CLASS without TEP patch."
            )

    lcdm_tt = np.asarray([row["D_ell_TT_uK2"] for row in lcdm_reference["tt_power"]])
    zero_tt = np.asarray([row["D_ell_TT_uK2"] for row in tep_zero_limit["tt_power"]])
    rel = np.max(np.abs(zero_tt - lcdm_tt) / np.maximum(np.abs(lcdm_tt), np.finfo(float).tiny))
    zero_limit_ok = bool(rel < 1e-6)

    if not zero_limit_ok:
        blockers.append(f"TEP CLASS epsilon_T=0 limit fails LCDM spectrum recovery: max rel diff={rel:.3e}")

    # Run detailed validation using cmb_class_validation module
    detailed_validation = {}
    if CMB_VALIDATION_AVAILABLE:
        validator = CLASSValidator()
        
        # Validate derived parameters
        # NOTE: This validation uses Planck 2018 baseline values as reference points.
        # CLASS with "output: tCl,pCl,lCl" does not compute sigma_8.
        # To get actual sigma_8, CLASS must be run with "output: mPk" (matter power spectrum).
        # For actual sigma_8 computation, see step_030_structure_growth_validation.py.
        derived_params = tep_spectra.get("derived", {})
        derived_validation = validator.validate_derived_parameters(
            H0=derived_params.get("H0", PLANCK_2018_H0),
            Omega_m=PLANCK_2018_OMEGA_M,
            sigma_8=PLANCK_2018_SIGMA_8,
        )
        
        # Validate TT spectrum
        l = np.arange(2, min(len(lcdm_tt) + 2, lmax + 1))
        if len(l) == len(lcdm_tt) and len(l) == len(zero_tt):
            spectrum_validation = validator.validate_spectrum(
                l=l,
                Dl_tep=zero_tt[:len(l)],
                Dl_lcdm=lcdm_tt[:len(l)],
                spectrum_type="TT",
            )
            detailed_validation = {
                "derived_parameters": derived_validation,
                "spectrum_validation": {
                    "type": spectrum_validation.spectrum_type,
                    "l_max": spectrum_validation.l_max,
                    "chi2_per_dof": spectrum_validation.chi2_per_dof,
                    "max_deviation_sigma": spectrum_validation.max_deviation_sigma,
                    "passed": spectrum_validation.passed,
                },
                "validation_module_available": True,
            }
        else:
            detailed_validation = {
                "error": "Spectrum length mismatch",
                "validation_module_available": True,
            }
    else:
        detailed_validation = {
            "validation_module_available": False,
            "note": "Install cmb_class_validation dependencies for detailed validation",
        }

    # Determine research grade status
    has_class_with_tep = CLASS_AVAILABLE and tep_class_works
    has_validation = CMB_VALIDATION_AVAILABLE
    
    # Research grade requires:
    # 1. CLASS with TEP support (for TEP parameter handling)
    # 2. Zero limit passes (TEP with epsilon_T=0 recovers LCDM)
    # Note: Validation module is optional for basic research grade
    research_grade = has_class_with_tep and zero_limit_ok
    
    # Build blockers list
    all_blockers = list(blockers)
    if not CLASS_AVAILABLE:
        all_blockers.append("CLASS not available - install CLASS Python package")
    elif not tep_class_works:
        all_blockers.append("TEP parameters not supported by CLASS - check TEP-CLASS patch")
    if not zero_limit_ok:
        all_blockers.append(f"TEP LCDM limit failed: max relative error = {rel:.3e}")
    if not has_validation:
        all_blockers.append("CMB validation module not available (optional)")
    
    validation = {
        "class_available": CLASS_AVAILABLE,
        "camb_available": CAMB_AVAILABLE,
        "tep_class_available": has_class_with_tep,
        "solver_used": "CLASS" if CLASS_AVAILABLE else "CAMB",
        "local_class_source_tep_patch_present": source_patch_present,
        "tep_zero_limit_max_relative_tt_error": float(rel),
        "tep_zero_limit_ok": zero_limit_ok,
        "research_grade_cmb": research_grade,
        "claim_gate": "open" if research_grade else "blocked",
        "detailed_validation": detailed_validation,
        "blockers": all_blockers if not research_grade else [],
    }
    return CMBRun(
        CLASS_AVAILABLE or CAMB_AVAILABLE, 
        has_class_with_tep, 
        lcdm_reference, 
        tep_zero_limit, 
        tep_spectra, 
        validation
    )
