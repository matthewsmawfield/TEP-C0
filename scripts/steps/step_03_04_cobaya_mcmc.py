#!/usr/bin/env python3
"""Step 033: Cobaya-based TEP inference with CLASS backend.

This step integrates the TEP-CLASS v2.0 modifications (ε_T, z_T, n_T parameters)
with Cobaya for joint SNe + CMB parameter estimation. It provides an alternative
to the emcee-based inference in step_018 with proper Boltzmann-solver backing.

Tier: RESEARCH GRADE

Requirements:
  - TEP-CLASS v2.0 built at /tmp/class_tep with TEP parameters enabled
  - Pantheon+ data in data/raw/
  - Cobaya installed with classy theory support

References:
  - TEP-CLASS patch: external/class_tep_mod/
  - Cobaya docs: https://cobaya.readthedocs.io/
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from c0_common import TEPLogger, ensure_dirs, print_status, step_json_path, write_json, set_step_logger

# Import Pantheon data handling from step_022
try:
    from step_03_01_three_model_comparison import PantheonData, RESEARCH_GRADE_RHAT_MAX
    HAS_STEP022 = True
except ImportError:
    HAS_STEP022 = False

STEP_ID = "step_03_04_cobaya_mcmc"

# TEP-CLASS path - built in external/class
# Path resolution: script -> steps -> scripts -> TEP-C0
PROJECT_ROOT = Path(__file__).resolve().parents[2]  # TEP-C0 directory
TEP_CLASS_PATH = PROJECT_ROOT / "external" / "class"
TEP_CLASS_BUILD = TEP_CLASS_PATH / "build" / "lib.macosx-11.1-arm64-cpython-313"


def check_tep_class_available() -> Tuple[bool, str]:
    """Check if TEP-CLASS v2.0 is available."""
    if not TEP_CLASS_PATH.exists():
        return False, f"TEP-CLASS not found at {TEP_CLASS_PATH}"
    
    if not TEP_CLASS_BUILD.exists():
        return False, f"TEP-CLASS build not found at {TEP_CLASS_BUILD}"
    
    # Try importing classy from build directory
    original_path = sys.path.copy()
    try:
        sys.path.insert(0, str(TEP_CLASS_BUILD))
        from classy import Class
        return True, "TEP-CLASS v2.0 available"
    except ImportError as e:
        return False, f"Cannot import classy from {TEP_CLASS_BUILD}: {e}"
    finally:
        sys.path = original_path


def check_cobaya_available() -> Tuple[bool, str]:
    """Check if Cobaya is installed."""
    try:
        import cobaya
        return True, f"Cobaya {cobaya.__version__} available"
    except ImportError:
        return False, "Cobaya not installed"


def create_cobaya_config(
    output_prefix: str,
    use_planck: bool = True,
    max_samples: int = 5000,
) -> Dict[str, Any]:
    """Create Cobaya configuration for TEP cosmology."""
    
    # Base theory configuration
    theory = {
        "classy": {
            "path": str(TEP_CLASS_PATH),
            "extra_args": {
                "N_ur": 2.0328,
                "N_ncdm": 1,
                "m_ncdm": 0.06,
                "output": "tCl,pCl,lCl,mPk",
                "P_k_max_h/Mpc": 10,
                "l_max_scalars": 2500,
                "lensing": "yes",
                "tep_mode": "yes",
                "non_linear": "halofit",
            }
        }
    }
    
    # Add Planck likelihoods if requested
    likelihood: Dict[str, Any] = {}
    if use_planck:
        try:
            # Full Planck 2018 TTTEEE + lowl
            likelihood["planck_2018_highl_plik.TTTEEE"] = {}
            likelihood["planck_2018_lowl.TT"] = {}
            likelihood["planck_2018_lowl.EE"] = {}
        except Exception:
            print_status("Planck likelihoods setup error", "WARNING")
    
    # Always add Pantheon+ (via custom likelihood)
    # This will be loaded from core/ directory
    likelihood["core.pantheon_cobaya_likelihood.PantheonCobaya"] = {}
    
    # Parameters
    params = {
        # TEP parameters
        "tep_epsilon_T": {
            "prior": {"min": -0.05, "max": 0.05},
            "ref": 0.001,
            "proposal": 0.00005,
            "latex": r"\epsilon_T",
        },
        "tep_z_T": {
            "prior": {"min": 0.5, "max": 15.0},
            "ref": 3.0,
            "proposal": 0.01,
            "latex": r"z_T",
        },
        "tep_n_T": {
            "value": 1.0,
            "latex": r"n_T",
        },
        # LCDM parameters
        "H0": {
            "prior": {"min": 60.0, "max": 80.0},
            "ref": 67.5,
            "proposal": 0.01,
            "latex": r"H_0",
        },
        "omega_b": {
            "prior": {"min": 0.018, "max": 0.025},
            "ref": 0.0224,
            "proposal": 0.00001,
            "latex": r"\Omega_b h^2",
        },
        "omega_cdm": {
            "prior": {"min": 0.10, "max": 0.14},
            "ref": 0.120,
            "proposal": 0.0001,
            "latex": r"\Omega_{cdm} h^2",
        },
        "tau_reio": {
            "prior": {"min": 0.01, "max": 0.10},
            "ref": 0.054,
            "proposal": 0.001,
            "latex": r"\tau_{reio}",
        },
        "logA": {
            "prior": {"min": 2.5, "max": 3.5},
            "ref": 3.044,
            "proposal": 0.001,
            "latex": r"\ln(10^{10}A_s)",
            "drop": True,
        },
        "A_s": {
            "value": "lambda logA: 1e-10 * np.exp(logA)",
            "latex": r"A_s",
        },
        "n_s": {
            "prior": {"min": 0.90, "max": 1.05},
            "ref": 0.966,
            "proposal": 0.0001,
            "latex": r"n_s",
        },
    }
    
    # Sampler configuration for production run
    sampler = {
        "mcmc": {
            "max_tries": 100000,
            "burn_in": 0,
            "Rminus1_stop": 0.02,  # 0.02 is standard publication threshold
            "Rminus1_cl_stop": 0.2,
            "covmat": "auto",
            "learn_every": "40d",
            "proposal_scale": 2.4,
            "max_samples": max_samples, # Safety limit to prevent unbounded memory growth
        }
    }
    
    config = {
        "output": output_prefix,
        "resume": True,
        "force": False,
        "theory": theory,
        "likelihood": likelihood,
        "params": params,
        "sampler": sampler,
        "debug": False,
        "verbose": 2,
    }
    
    return config


def sanity_check_tep_active() -> bool:
    """Verify that TEP parameters are actively modifying the CLASS output."""
    print_status("Running TEP sanity check...", "PROCESS")
    original_path = sys.path.copy()
    try:
        sys.path.insert(0, str(TEP_CLASS_BUILD))
        from classy import Class
        
        cosmo1 = Class()
        cosmo1.set({'tep_mode': 'yes', 'tep_epsilon_T': 0.0, 'tep_z_T': 3.0, 'output': 'tCl,lCl,mPk'})
        cosmo1.compute()
        d1 = cosmo1.luminosity_distance(1.0)
        
        cosmo2 = Class()
        cosmo2.set({'tep_mode': 'yes', 'tep_epsilon_T': 0.2, 'tep_z_T': 3.0, 'output': 'tCl,lCl,mPk'})
        cosmo2.compute()
        d2 = cosmo2.luminosity_distance(1.0)
        
        if np.isclose(d1, d2, rtol=1e-5):
            print_status("TEP sanity check failed: epsilon_T has no effect on distances!", "ERROR")
            return False
            
        print_status(f"TEP sanity check passed: dL(eps=0)={d1:.2f}, dL(eps=0.2)={d2:.2f}", "SUCCESS")
        return True
    except Exception as e:
        print_status(f"TEP sanity check failed with exception: {e}", "ERROR")
        return False
    finally:
        sys.path = original_path



def run_cobaya_mcmc(
    config: Dict[str, Any],
) -> Dict[str, Any]:
    """Run Cobaya MCMC with given configuration."""
    try:
        from mpi4py import MPI
        rank = MPI.COMM_WORLD.Get_rank()
    except ImportError:
        rank = 0

    try:
        from cobaya.run import run
        import contextlib
        import sys
        
        log_path = Path(f"logs/{STEP_ID}.log")
    except ImportError:
        return {"success": False, "error": "Cobaya not available"}
    
    # Set environment for TEP-CLASS
    os.environ["TEP_CLASS_PYTHONPATH"] = str(TEP_CLASS_BUILD)
    os.environ["PYTHONPATH"] = str(TEP_CLASS_BUILD) + os.pathsep + os.environ.get("PYTHONPATH", "")
    
    result = {"success": False, "samples": None, "error": None}
    
    try:
        if rank == 0:
            # Run Cobaya and redirect its stdout/stderr to the log file on rank 0
            with open(log_path, 'a') as f:
                with contextlib.redirect_stdout(f), contextlib.redirect_stderr(f):
                    updated_info = run(config, resume=True)
        else:
            # Discard stdout/stderr on other ranks
            with open(os.devnull, 'w') as f:
                with contextlib.redirect_stdout(f), contextlib.redirect_stderr(f):
                    updated_info = run(config, resume=True)
                    
        if rank == 0:
            result["success"] = True
            result["updated_info"] = updated_info
            # Try to extract samples
            try:
                from cobaya.output import load_samples
                samples = load_samples(config["output"])
                result["samples"] = samples
            except Exception as e:
                result["samples_error"] = str(e)
    
    except Exception as e:
        if rank == 0:
            result["error"] = str(e)
    
    return result


def run() -> dict:
    """Run Cobaya-based TEP inference."""
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()
    
    # Check prerequisites
    tep_class_ok, tep_class_msg = check_tep_class_available()
    cobaya_ok, cobaya_msg = check_cobaya_available()
    
    print_status(f"TEP-CLASS: {tep_class_msg}", "SUCCESS" if tep_class_ok else "WARNING")
    print_status(f"Cobaya: {cobaya_msg}", "SUCCESS" if cobaya_ok else "WARNING")
    
    if not tep_class_ok or not cobaya_ok:
        print_status("Skipping Cobaya inference - dependencies not available", "INFO")
        payload = {
            "step": STEP_ID,
            "status": "skipped",
            "reason": "Missing dependencies",
            "tep_class_available": tep_class_ok,
            "cobaya_available": cobaya_ok,
        }
        write_json(step_json_path(STEP_ID), payload)
        return payload
    
    # Check Pantheon data availability
    pantheon_available = False
    if HAS_STEP022:
        try:
            data = PantheonData()
            data.load()
            pantheon_available = True
            print_status(f"Pantheon+ data: {len(data.z)} SNe loaded", "SUCCESS")
        except Exception as e:
            print_status(f"Pantheon+ data error: {e}", "WARNING")
    
    # Create output directory
    output_dir = Path("results/outputs")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = str(output_dir / "tep_cobaya_sne")
    
    # Create configuration
    max_samples = int(os.getenv("TEP_COBAYA_SAMPLES", "500000"))
    
    # Run sanity check
    if not sanity_check_tep_active():
        print_status("Aborting inference due to failed sanity check.", "ERROR")
        payload = {"step": STEP_ID, "status": "failed", "error": "TEP sanity check failed"}
        write_json(step_json_path(STEP_ID), payload)
        return payload

    config = create_cobaya_config(
        output_prefix=output_prefix,
        use_planck=True,  # Enable Planck for joint inference
        max_samples=max_samples,
    )
    
    try:
        from mpi4py import MPI
        rank = MPI.COMM_WORLD.Get_rank()
    except ImportError:
        rank = 0

    if rank == 0:
        print_status(f"Running Cobaya MCMC (max_samples={max_samples})...", "PROCESS")
    
    # Run MCMC without timeout constraint
    mcmc_result = run_cobaya_mcmc(config)
    
    # Process results
    if mcmc_result.get("success"):
        print_status("Cobaya MCMC completed successfully", "SUCCESS")
        status = "completed"
    elif mcmc_result.get("partial"):
        print_status("Cobaya MCMC partial (timeout)", "WARNING")
        status = "partial"
    else:
        print_status(f"Cobaya MCMC failed: {mcmc_result.get('error')}", "ERROR")
        status = "failed"
    
    # Prepare output
    payload = {
        "step": STEP_ID,
        "description": "Cobaya-based TEP inference with TEP-CLASS v2.0",
        "status": status,
        "tep_class_available": tep_class_ok,
        "cobaya_available": cobaya_ok,
        "pantheon_data_available": pantheon_available,
        "config": {
            "output_prefix": output_prefix,
            "max_samples": max_samples,
        },
        "mcmc_result": {
            "success": mcmc_result.get("success", False),
            "error": mcmc_result.get("error"),
        },
        "validation": {
            "can_run_joint_analysis": tep_class_ok and cobaya_ok and pantheon_available,
            "tep_class_path": str(TEP_CLASS_PATH),
        },
    }
    
    if rank == 0:
        write_json(step_json_path(STEP_ID), payload)
        print_status(f"Step {STEP_ID} completed with status: {status}", 
                     "SUCCESS" if status == "completed" else "WARNING")
    
    return payload


if __name__ == "__main__":
    run()
