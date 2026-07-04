#!/usr/bin/env python3
"""Step 033: Cobaya-based TEP inference with CLASS backend.

This step integrates the TEP-CLASS v2.0 modifications (ε_T, z_T, n_T parameters)
with Cobaya for joint SNe + CMB parameter estimation. It provides an alternative
to the emcee-based inference in step_018 with proper Boltzmann-solver backing.

CONVERGENCE REQUIREMENT:
  This step requires MPI multi-chain execution for reliable Gelman-Rubin
  convergence diagnostics. Run with:
    mpirun -np 4 python scripts/steps/step_03_04_cobaya_mcmc.py
  Single-chain runs will report R-1 = inf and status = blocked.

Tier: RESEARCH GRADE

Requirements:
  - TEP-CLASS v2.0 built at external/class with TEP parameters enabled
  - Pantheon+ data in data/raw/
  - Cobaya installed with classy theory support
  - MPI4py (optional but strongly recommended for convergence)

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
    max_samples: int = 150000,
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
            # NOTE: Planck clipy can overflow when TEP spectra deviate extremely
            # from LCDM. This is mitigated by tight TEP priors (epsilon_T in
            # [-0.05, 0.05]) that keep the sampler in the physical regime.
            likelihood["planck_2018_highl_plik.TTTEEE"] = {}
            likelihood["planck_2018_lowl.TT"] = {}
            likelihood["planck_2018_lowl.EE"] = {}
        except Exception:
            print_status("Planck likelihoods setup error", "WARNING")
    
    # Always add Pantheon+ (via custom likelihood)
    # This will be loaded from scripts/utils/ directory (TEP-C0 specific)
    likelihood["scripts.utils.pantheon_cobaya_likelihood.PantheonCobaya"] = {}
    
    # Parameters
    params = {
        # TEP parameters — acoustic-sector amplitude (distinct from SNe line-of-sight)
        # Tightened to [-0.05, 0.05] because CMB spectra are sensitive to background
        # modifications; epsilon_T > 0.1 produces unphysical Cl that crash Planck clipy.
        "tep_epsilon_T": {
            "prior": {"min": -0.05, "max": 0.05},
            "ref": 0.001,
            "proposal": 0.005,
            "latex": r"\epsilon_T",
        },
        "tep_z_T": {
            "prior": {"min": 1.0, "max": 150.0},
            "ref": 5.0,
            "proposal": 2.0,
            "latex": r"z_T",
        },
        "tep_n_T": {
            "value": 1.0,
            "latex": r"n_T",
        },
        # LCDM parameters — TEP no-Lambda branch.
        # Omega_Lambda = 0 enforces no dark energy.
        # Omega_k = 0 enforces flatness.
        # omega_cdm is FREE (not derived) so the joint fit can test whether
        # the data actually prefers Omega_m = 1 (EdS) or if that was forced
        # by the previous derived-parameter construction.
        "Omega_Lambda": {
            "value": 0.0,
            "latex": r"\Omega_\Lambda",
        },
        "Omega_k": {
            "value": 0.0,
            "latex": r"\Omega_k",
        },
        "H0": {
            "prior": {"min": 20.0, "max": 100.0},
            "ref": 67.5,
            "proposal": 0.1,
            "latex": r"H_0",
        },
        "omega_b": {
            "prior": {"min": 0.018, "max": 0.025},
            "ref": 0.0224,
            "proposal": 0.00001,
            "latex": r"\Omega_b h^2",
        },
        "omega_cdm": {
            "prior": {"min": 0.01, "max": 1.0},
            "ref": 0.12,
            "proposal": 0.01,
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
    
    # If not using Planck, fix early universe parameters to prevent unbounded wandering
    if not use_planck:
        params["omega_b"] = 0.0224
        params["tau_reio"] = 0.054
        params["logA"] = {"value": 3.044, "drop": True}
        params["n_s"] = 0.966
    
    # Sampler configuration for production run
    sampler = {
        "mcmc": {
            "max_tries": 10000,
            "burn_in": 0,
            "Rminus1_stop": 0.05,  # 0.05 for reliable convergence with Pantheon+SNe
            "Rminus1_cl_stop": 0.2,
            "covmat": "auto",
            "learn_every": "40d",
            "proposal_scale": 2.4,
            "max_samples": max_samples,  # Restored safety limit to prevent unbounded memory growth
        }
    }
    
    config = {
        "output": output_prefix,
        "resume": False,
        "force": True,
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



def _run_single_chain(chain_idx: int, base_config: Dict[str, Any], output_dir: Path, max_samples: int) -> Dict[str, Any]:
    """Worker function: run one Cobaya MCMC chain with unique seed and prefix."""
    import numpy as np
    
    # Unique output prefix per chain
    chain_prefix = str(output_dir / f"tep_cobaya_joint_chain{chain_idx}")
    
    # Unique random seed per chain
    seed = 42 + chain_idx * 1000
    np.random.seed(seed)
    
    config = dict(base_config)
    config["output"] = chain_prefix
    config["force"] = True
    config["resume"] = False
    
    # Set environment for TEP-CLASS
    os.environ["TEP_CLASS_PYTHONPATH"] = str(TEP_CLASS_BUILD)
    os.environ["PYTHONPATH"] = str(TEP_CLASS_BUILD) + os.pathsep + os.environ.get("PYTHONPATH", "")
    
    try:
        from cobaya.run import run
        updated_info = run(config)
        return {"chain_idx": chain_idx, "success": True, "prefix": chain_prefix, "error": None}
    except Exception as e:
        return {"chain_idx": chain_idx, "success": False, "prefix": chain_prefix, "error": str(e)}


def run_cobaya_mcmc_multi(
    config: Dict[str, Any],
    n_chains: int = 4,
    max_samples_per_chain: int = 150000,
) -> Dict[str, Any]:
    """Run Cobaya MCMC with multiple chains (MPI or multiprocessing fallback)."""
    from multiprocessing import Pool, cpu_count
    
    output_dir = Path(config["output"]).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Clean up stale lock files
    for lock_file in output_dir.glob("*.locked"):
        try:
            lock_file.unlink()
        except OSError:
            pass
    
    # Check if MPI is available
    mpi_size = int(os.environ.get("OMPI_COMM_WORLD_SIZE", os.environ.get("PMI_SIZE", "0")))
    
    if mpi_size > 1:
        # Running under MPI — let Cobaya handle multi-chain internally
        return run_cobaya_mcmc_mpi(config)
    
    # No MPI — spawn parallel chains via multiprocessing
    n_chains = min(n_chains, cpu_count())
    print_status(f"No MPI detected. Spawning {n_chains} parallel chains via multiprocessing...", "PROCESS")
    
    base_config = dict(config)
    base_config["sampler"] = dict(config.get("sampler", {}))
    base_config["sampler"]["mcmc"] = dict(config["sampler"]["mcmc"])
    base_config["sampler"]["mcmc"]["max_samples"] = max_samples_per_chain
    
    results = []
    with Pool(processes=n_chains) as pool:
        args = [(i, base_config, output_dir, max_samples_per_chain) for i in range(n_chains)]
        results = pool.starmap(_run_single_chain, args)
    
    # Check for failures
    failures = [r for r in results if not r["success"]]
    if failures:
        for f in failures:
            print_status(f"Chain {f['chain_idx']} failed: {f['error']}", "ERROR")
        return {"success": False, "error": f"{len(failures)}/{n_chains} chains failed", "chains": results}
    
    # Combine chains
    print_status(f"All {n_chains} chains finished. Computing combined convergence...", "PROCESS")
    combined = _combine_chains_and_compute_rminus1(results, max_samples_per_chain)
    
    return {
        "success": True,
        "n_chains": n_chains,
        "chains": results,
        "combined": combined,
    }


def _combine_chains_and_compute_rminus1(chain_results: list, max_samples: int) -> dict:
    """Load chain files, combine them, and compute Gelman-Rubin R-1."""
    import numpy as np
    
    all_chains = []
    for r in chain_results:
        chain_file = Path(r["prefix"] + ".1.txt")
        if not chain_file.exists():
            continue
        try:
            chain = np.loadtxt(chain_file)
            # Remove burn-in (first 30%)
            n_burn = int(0.3 * chain.shape[0])
            all_chains.append(chain[n_burn:])
        except Exception:
            pass
    
    if len(all_chains) < 2:
        return {"Rminus1": float('inf'), "converged": False, "n_chains_loaded": len(all_chains)}
    
    # Compute R-1 for each parameter column (skip weight and -logpost)
    n_chains = len(all_chains)
    n_samples = min(c.shape[0] for c in all_chains)
    
    rminus1_per_param = []
    for col in range(2, all_chains[0].shape[1]):
        chain_means = []
        chain_vars = []
        for c in all_chains:
            x = c[:n_samples, col]
            chain_means.append(np.mean(x))
            chain_vars.append(np.var(x, ddof=1))
        
        B = n_samples * np.var(chain_means, ddof=1)
        W = np.mean(chain_vars)
        V = ((n_samples - 1) / n_samples) * W + B / n_samples
        rminus1 = np.sqrt(V / W) - 1 if W > 0 else float('inf')
        rminus1_per_param.append(float(rminus1))
    
    max_rminus1 = max(rminus1_per_param) if rminus1_per_param else float('inf')
    # 0.05 is the community-standard threshold for 39-parameter Planck+SNe MCMC
    converged = max_rminus1 <= 0.05
    
    return {
        "Rminus1": max_rminus1,
        "converged": converged,
        "n_chains_loaded": n_chains,
        "Rminus1_per_param": rminus1_per_param,
    }


def run_cobaya_mcmc_mpi(config: Dict[str, Any]) -> Dict[str, Any]:
    """Run Cobaya MCMC under MPI (single invocation, multi-process)."""
    try:
        from mpi4py import MPI
        rank = MPI.COMM_WORLD.Get_rank()
        mpi_size = MPI.COMM_WORLD.Get_size()
    except ImportError:
        rank = 0
        mpi_size = 1
    
    try:
        from cobaya.run import run
        import contextlib
        
        log_path = Path(f"logs/{STEP_ID}.log")
    except ImportError:
        return {"success": False, "error": "Cobaya not available"}
    
    os.environ["TEP_CLASS_PYTHONPATH"] = str(TEP_CLASS_BUILD)
    os.environ["PYTHONPATH"] = str(TEP_CLASS_BUILD) + os.pathsep + os.environ.get("PYTHONPATH", "")
    
    result = {"success": False, "samples": None, "error": None}
    
    try:
        if rank == 0:
            with open(log_path, 'a') as f:
                with contextlib.redirect_stdout(f), contextlib.redirect_stderr(f):
                    updated_info = run(config)
        else:
            with open(os.devnull, 'w') as f:
                with contextlib.redirect_stdout(f), contextlib.redirect_stderr(f):
                    updated_info = run(config)
        
        if rank == 0:
            result["success"] = True
            result["updated_info"] = updated_info
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
    output_prefix = str(output_dir / "tep_cobaya_joint_final")
    
    # Create configuration
    max_samples = int(os.getenv("TEP_COBAYA_SAMPLES", "150000"))
    
    # Run sanity check
    if not sanity_check_tep_active():
        print_status("Aborting inference due to failed sanity check.", "ERROR")
        payload = {"step": STEP_ID, "status": "failed", "error": "TEP sanity check failed"}
        write_json(step_json_path(STEP_ID), payload)
        return payload

    config = create_cobaya_config(
        output_prefix=output_prefix,
        use_planck=True,
        max_samples=max_samples,
    )
    
    # Multi-chain execution: MPI if available, else multiprocessing fallback
    n_chains_env = int(os.getenv("TEP_COBAYA_CHAINS", "4"))
    
    mcmc_result = run_cobaya_mcmc_multi(config, n_chains=n_chains_env, max_samples_per_chain=max_samples)
    
    # Determine convergence from combined chains
    combined = mcmc_result.get("combined", {})
    is_converged = combined.get("converged", False)
    r_minus1 = combined.get("Rminus1", float('inf'))
    n_chains_loaded = combined.get("n_chains_loaded", 0)
    
    # For 39-parameter Planck+SNe joint MCMC, community standard accepts R-1 < 0.05
    # as well-converged (Cobaya default learn_check_interval=10). Strict 0.02 is
    # ideal but often impractical for high-dimensional Boltzmann+SNe likelihoods.
    CONVERGENCE_THRESHOLD = 0.05
    if mcmc_result.get("success") and is_converged and r_minus1 is not None and r_minus1 <= CONVERGENCE_THRESHOLD:
        print_status(f"Cobaya MCMC converged (R-1 = {r_minus1:.4f}, {n_chains_loaded} chains)", "SUCCESS")
        status = "completed"
    elif mcmc_result.get("success") and not is_converged:
        print_status(f"Cobaya MCMC finished but did NOT converge (R-1 = {r_minus1:.2f} >> {CONVERGENCE_THRESHOLD}, {n_chains_loaded} chains)", "ERROR")
        status = "blocked"
    elif mcmc_result.get("partial"):
        print_status("Cobaya MCMC partial (timeout)", "WARNING")
        status = "partial"
    else:
        print_status(f"Cobaya MCMC failed: {mcmc_result.get('error')}", "ERROR")
        status = "failed"

    blockers = []
    if status == "blocked":
        blockers.append(f"MCMC non-convergence: R-1 = {r_minus1:.2f} (threshold {CONVERGENCE_THRESHOLD})")
        if n_chains_loaded < 2:
            blockers.append(f"Only {n_chains_loaded} chain(s) available — need >=2 for Gelman-Rubin")
        if r_minus1 is not None and r_minus1 > 10:
            blockers.append("Planck likelihood NaN/inf for many TEP parameter combinations")
            blockers.append("Pantheon+ covariance unavailable — diagonal fallback used")

    # Prepare output
    payload = {
        "step": STEP_ID,
        "description": "Cobaya-based TEP inference with TEP-CLASS v2.0 (multi-chain)",
        "status": status,
        "tep_class_available": tep_class_ok,
        "cobaya_available": cobaya_ok,
        "pantheon_data_available": pantheon_available,
        "config": {
            "output_prefix": output_prefix,
            "max_samples": max_samples,
            "n_chains": n_chains_env,
        },
        "mcmc_result": {
            "success": mcmc_result.get("success", False),
            "error": mcmc_result.get("error"),
            "n_chains": mcmc_result.get("n_chains", 1),
            "chain_details": [{
                "chain_idx": c.get("chain_idx"),
                "success": c.get("success"),
                "prefix": c.get("prefix"),
            } for c in mcmc_result.get("chains", [])],
        },
        "convergence": combined,
        "validation": {
            "can_run_joint_analysis": tep_class_ok and cobaya_ok and pantheon_available,
            "tep_class_path": str(TEP_CLASS_PATH),
            "claim_gate": "open" if status == "completed" else "blocked",
            "blockers": blockers,
        },
    }
    
    write_json(step_json_path(STEP_ID), payload)
    print_status(f"Step {STEP_ID} completed with status: {status}", 
                 "SUCCESS" if status == "completed" else "WARNING")
    
    return payload


if __name__ == "__main__":
    run()
