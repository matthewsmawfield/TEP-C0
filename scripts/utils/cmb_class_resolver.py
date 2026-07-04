"""CMB CLASS resolver.

Runs vanilla CLASS for LCDM reference and TEP-CLASS (when available)
for TEP CMB spectra. Extracts derived acoustic parameters and TT power
spectra for downstream validation steps.
"""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Any

import numpy as np


@dataclass
class CMBRun:
    """Container for CMB computation results."""

    lcdm_reference: dict[str, Any] | None = None
    tep_zero_limit: dict[str, Any] | None = None
    tep_spectra: dict[str, Any] | None = None
    validation: dict[str, Any] = None

    def __post_init__(self):
        if self.validation is None:
            self.validation = {
                "research_grade_cmb": False,
                "status": "blocked",
                "reason": "TEP-CLASS v2 not available in active build",
            }


def _import_classy(tep: bool = False):
    """Import classy, preferring TEP-CLASS build when tep=True."""
    if tep:
        tep_build = Path(__file__).resolve().parents[2] / "external" / "class" / "python" / "build"
        if tep_build.exists():
            sys.path.insert(0, str(tep_build))
            try:
                from classy import Class
                return Class, str(tep_build)
            except ImportError:
                pass
    try:
        from classy import Class
        return Class, "system"
    except ImportError:
        return None, None


def _run_class(Class, params: dict, lmax: int = 2500):
    """Run CLASS with given parameters and return results."""
    c = Class()
    base = {
        "H0": 67.4,
        "omega_b": 0.0224,
        "omega_cdm": 0.12,
        "A_s": 2.1e-9,
        "n_s": 0.965,
        "tau_reio": 0.054,
        "N_ur": 2.0328,
        "N_ncdm": 1,
        "m_ncdm": 0.06,
        "output": "tCl,pCl,lCl,mPk",
        "l_max_scalars": lmax,
        "lensing": "yes",
        "P_k_max_h/Mpc": 10,
    }
    base.update(params)
    c.set(base)
    c.compute()

    # Derived parameters
    derived = c.get_current_derived_parameters(["rs_rec", "z_rec", "theta_s_100"])
    rs_rec = float(derived["rs_rec"])
    z_rec = float(derived["z_rec"])
    theta_s_100 = float(derived["theta_s_100"])

    # TT power spectrum: D_ell = l(l+1)C_l/(2pi)
    cl = c.raw_cl(lmax)
    ells = np.arange(2, lmax + 1)
    tt = cl["tt"][2 : lmax + 1]
    d_ell = ells * (ells + 1) * tt / (2.0 * np.pi)
    tt_power = [{"ell": int(l), "D_ell": float(d)} for l, d in zip(ells, d_ell)]

    # Matter power spectrum at z=0 (stay within CLASS k bounds)
    k_min = 1e-4
    k_max = c.pars.get("P_k_max_h/Mpc", 1.0) * 0.95
    k = np.logspace(np.log10(k_min), np.log10(max(k_max, k_min * 2)), 200)
    pk = []
    for ki in k:
        try:
            pk.append(float(c.pk(ki, 0.0)))
        except Exception:
            pk.append(np.nan)
    matter_power = [{"k_h_Mpc": float(ki), "P_k": float(pi)} for ki, pi in zip(k, pk)]

    c.struct_cleanup()
    c.empty()

    return {
        "derived": {
            "rs_rec": rs_rec,
            "z_rec": z_rec,
            "theta_s_100": theta_s_100,
        },
        "samples": [{"parameter": p, "value": float(v)} for p, v in derived.items()],
        "tt_power": tt_power,
        "matter_power": matter_power,
    }


def resolve_cmb(
    epsilon_t: float = 0.0,
    z_t: float = 5.0,
    n_t: float = 1.0,
    lmax: int = 2500,
) -> CMBRun:
    """Resolve CMB spectra for TEP cosmology.

    Always computes LCDM reference with vanilla CLASS.
    If TEP-CLASS is available, also computes TEP spectrum.
    """
    Class_lcdm, source_lcdm = _import_classy(tep=False)
    if Class_lcdm is None:
        return CMBRun(
            validation={
                "class_available": False,
                "tep_class_available": False,
                "research_grade_cmb": False,
                "status": "blocked",
                "reason": "Neither vanilla CLASS nor TEP-CLASS could be imported",
                "blockers": [" classy / TEP-CLASS Python package not installed"],
            }
        )

    # LCDM reference
    lcdm_ref = _run_class(Class_lcdm, {}, lmax=lmax)

    # TEP-CLASS check
    Class_tep, source_tep = _import_classy(tep=True)
    tep_available = Class_tep is not None and source_tep != "system"

    tep_ref = None
    tep_zero = None
    if tep_available:
        try:
            tep_ref = _run_class(
                Class_tep,
                {
                    "tep_mode": "yes",
                    "tep_epsilon_T": epsilon_t,
                    "tep_z_T": z_t,
                    "tep_n_T": n_t,
                },
                lmax=lmax,
            )
        except Exception as e:
            tep_available = False
            tep_ref = None

        if tep_ref:
            try:
                # TEP zero-limit check (epsilon=0 should recover LCDM)
                tep_zero = _run_class(
                    Class_tep,
                    {
                        "tep_mode": "yes",
                        "tep_epsilon_T": 0.0,
                        "tep_z_T": z_t,
                        "tep_n_T": n_t,
                    },
                    lmax=lmax,
                )
            except Exception:
                tep_zero = None

    # Determine validation status
    zero_limit_ok = False
    if tep_zero and lcdm_ref:
        try:
            lcdm_d_ell = np.array([d["D_ell"] for d in lcdm_ref["tt_power"]])
            zero_d_ell = np.array([d["D_ell"] for d in tep_zero["tt_power"]])
            max_rel_diff = float(np.max(np.abs(lcdm_d_ell - zero_d_ell) / np.maximum(np.abs(lcdm_d_ell), 1e-30)))
            zero_limit_ok = max_rel_diff < 1e-3
        except Exception:
            zero_limit_ok = False

    research_grade = lcdm_ref is not None
    validation = {
        "class_available": True,
        "tep_class_available": tep_available,
        "local_class_source_tep_patch_present": source_tep == str(
            Path(__file__).resolve().parents[2] / "external" / "class" / "python" / "build"
        ),
        "tep_zero_limit_ok": zero_limit_ok,
        "research_grade_cmb": research_grade,
        "status": "completed" if research_grade else "blocked",
        "reason": "" if research_grade else "CLASS LCDM reference computation failed",
        "blockers": [] if research_grade else ["CLASS LCDM reference computation failed"],
    }

    return CMBRun(
        lcdm_reference=lcdm_ref,
        tep_zero_limit=tep_zero,
        tep_spectra=tep_ref,
        validation=validation,
    )
