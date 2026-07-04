#!/usr/bin/env python3
"""Step 033b: Analyze existing Cobaya TEP chains.

Loads the completed/partial chains from step_033 and produces summary statistics
and triangle plots.
"""

import sys
from pathlib import Path
import getdist
from getdist import plots, MCSamples
from getdist.mcsamples import loadMCSamples
import numpy as np

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "utils"))
from plot_style import apply_tep_style

# Apply TEP manuscript style (affects matplotlib backend used by getdist)
apply_tep_style()

from c0_common import TEPLogger, ensure_dirs, print_status, set_step_logger

STEP_ID = "step_03_05_analyze_cobaya"

def run():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()
    
    # Support both old single-chain and new multi-chain outputs
    combined_prefix = "results/outputs/tep_cobaya_joint_combined"
    
    if Path(f"{combined_prefix}.1.txt").exists():
        # Use pre-combined multi-chain file
        chain_prefix = combined_prefix
    else:
        # Fall back to single chain
        chain_prefix = "results/outputs/tep_cobaya_joint_final"
    
    if not Path(f"{chain_prefix}.1.txt").exists():
        print_status(f"Chain files not found at {chain_prefix}.1.txt", "WARNING")
        return {"step": STEP_ID, "status": "skipped", "reason": "No chain files from cobaya MCMC"}
        
    print_status(f"Loading chains from {chain_prefix}...", "PROCESS")
    
    try:
        # Load samples, removing 30% as burn-in.
        # getdist loadMCSamples automatically reads .input.yaml for parameter names and labels.
        samples = loadMCSamples(chain_prefix, settings={'ignore_rows': 0.3})
        
        print_status("Chains loaded successfully.", "SUCCESS")
        
        # Define the parameters we actually care about plotting
        # Use logA if available (new runs), fall back to A_s (old runs)
        plot_params = ['tep_epsilon_T', 'tep_z_T', 'H0', 'omega_b', 'omega_cdm', 'tau_reio', 'logA', 'n_s']
        
        # Check if they exist in the chains
        avail_params = samples.getParamNames().list()
        plot_params = [p for p in plot_params if p in avail_params]
        
        # Print constraints
        stats = samples.getMargeStats()
        print("\n--- Parameter Constraints (68% limits) ---")
        for param in plot_params:
            par_stat = stats.parWithName(param)
            if par_stat:
                print(f"{param}: {par_stat.mean:.4f} +/- {par_stat.err:.4f}")
        print("--------------------------------------------\n")
        
        # Plot triangle
        print_status("Generating triangle plot...", "PROCESS")
        g = plots.get_subplot_plotter(subplot_size=2.8)
        g.settings.figure_legend_frame = True
        g.settings.alpha_filled_add = 0.85
        g.settings.title_limit_fontsize = 14
        g.settings.axes_fontsize = 12
        g.settings.lab_fontsize = 13
        g.settings.legend_fontsize = 11
        
        try:
            g.triangle_plot([samples], plot_params, filled=True)
            
            out_path = Path(f"results/figures/{STEP_ID}_triangle.png")
            out_path.parent.mkdir(parents=True, exist_ok=True)
            g.export(str(out_path))
            print_status(f"Plot saved to {out_path}", "SUCCESS")
            print_status(
                "NOTE: This is a joint SNe+CMB background/acoustic MCMC diagnostic. "
                "tep_epsilon_T here is the homogeneous acoustic-sector amplitude, "
                "distinct from epsilon_shear_los fitted to SNe alone. "
                "H0 boundary behaviour is separately stress-tested.",
                "INFO",
            )
        except Exception as e:
            print_status(f"Could not generate triangle plot (chains may lack variance): {e}", "WARNING")
        
        # Save a basic stats dump
        stats_path = Path("results/outputs/tep_cobaya_joint_converged_stats.txt")
        constraint_list = []
        with open(stats_path, 'w') as f:
            f.write("--- Parameter Constraints ---\n")
            for param in plot_params:
                par_stat = stats.parWithName(param)
                if par_stat:
                    f.write(f"{param}: {par_stat.mean:.4g} +/- {par_stat.err:.4g}\n")
                    constraint_list.append({
                        "parameter": param,
                        "mean": float(par_stat.mean),
                        "std": float(par_stat.err),
                        "lower_68": float(par_stat.limits[0].lower),
                        "upper_68": float(par_stat.limits[0].upper),
                    })
        print_status(f"Summary stats saved to {stats_path}", "SUCCESS")
        
        # Read convergence from step_03_04 JSON (multi-chain combined R-1)
        import json
        step04_path = Path("results/step_03_04_cobaya_mcmc.json")
        convergence_info = {}
        if step04_path.exists():
            try:
                with open(step04_path) as f:
                    step04 = json.load(f)
                convergence_info = step04.get("convergence", {})
            except Exception:
                pass

        from c0_common import step_json_path, write_json
        
        # Determine overall status based on convergence
        is_converged = convergence_info.get("converged", False)
        r_minus1 = convergence_info.get("Rminus1", float('inf'))
        status = "completed" if is_converged else "blocked"
        
        payload = {
            "step": STEP_ID,
            "status": status,
            "constraints": constraint_list,
            "n_parameters": len(plot_params),
            "figure_path": str(out_path) if Path(out_path).exists() else None,
            "stats_path": str(stats_path),
            "chain_files": [f.name for f in Path("results/outputs").glob("tep_cobaya_joint_chain*")],
            "convergence": convergence_info,
            "note": (
                "Joint SNe+CMB background/acoustic MCMC diagnostic. "
                "tep_epsilon_T here is the homogeneous acoustic-sector amplitude, "
                "distinct from epsilon_shear_los fitted to SNe alone. "
                f"R-1 = {r_minus1:.3f} (threshold 0.05 for 39-parameter Planck+SNe). "
                + ("Converged." if is_converged else "Not yet converged — constraints are preliminary.")
            ),
        }
        write_json(step_json_path(STEP_ID), payload)
        return payload
        
    except Exception as e:
        print_status(f"Error analyzing chains: {e}", "ERROR")
        import traceback
        traceback.print_exc()
        return {"step": STEP_ID, "status": "failed", "error": str(e)}

if __name__ == "__main__":
    run()
