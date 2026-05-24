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
from c0_common import TEPLogger, ensure_dirs, print_status, set_step_logger

STEP_ID = "step_033b_analyze_cobaya"

def run():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()
    
    chain_prefix = "results/outputs/tep_cobaya_sne"
    
    if not Path(f"{chain_prefix}.1.txt").exists():
        print_status(f"Chain files not found at {chain_prefix}.1.txt", "ERROR")
        return
        
    print_status(f"Loading chains from {chain_prefix}...", "PROCESS")
    
    try:
        # Load samples, removing 30% as burn-in.
        # getdist loadMCSamples automatically reads .input.yaml for parameter names and labels.
        samples = loadMCSamples(chain_prefix, settings={'ignore_rows': 0.3})
        
        print_status("Chains loaded successfully.", "SUCCESS")
        
        # Define the parameters we actually care about plotting
        plot_params = ['tep_epsilon_T', 'tep_z_T', 'H0', 'omega_b', 'omega_cdm', 'tau_reio', 'A_s', 'n_s']
        
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
        g = plots.get_subplot_plotter()
        g.settings.figure_legend_frame = False
        g.settings.alpha_filled_add = 0.85
        g.settings.title_limit_fontsize = 14
        
        g.triangle_plot([samples], plot_params, filled=True, title_limit=1)
        
        out_path = Path("results/figures/step_033_cobaya_tep_triangle.png")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        g.export(str(out_path))
        print_status(f"Plot saved to {out_path}", "SUCCESS")
        
        # Save a basic stats dump
        stats_path = Path("results/outputs/tep_cobaya_sne_stats.txt")
        with open(stats_path, 'w') as f:
            f.write("--- Parameter Constraints ---\n")
            for param in plot_params:
                par_stat = stats.parWithName(param)
                if par_stat:
                    f.write(f"{param}: {par_stat.mean:.4f} +/- {par_stat.err:.4f}\n")
        print_status(f"Summary stats saved to {stats_path}", "SUCCESS")
        
    except Exception as e:
        print_status(f"Error analyzing chains: {e}", "ERROR")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    run()
