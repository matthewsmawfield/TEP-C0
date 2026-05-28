import sys
sys.path.insert(0, ".")
from getdist import loadMCSamples
import getdist.plots as gdplt
import matplotlib.pyplot as plt
import numpy as np

# Load samples from the 4 chains
# By default, loadMCSamples looks for root_name_*.txt
try:
    samples = loadMCSamples('results/outputs/tep_cobaya_sne', settings={'ignore_rows': 0.3})
    
    stats = samples.getMargeStats()
    
    print("=== MCMC Convergence ===")
    print(f"R-1 (Gelman-Rubin): {samples.getGelmanRubin()}")
    
    print("\n=== Parameter Constraints (68% limits) ===")
    params = ['tep_epsilon_T', 'tep_z_T', 'H0', 'omega_cdm', 'n_s']
    for p in params:
        try:
            lims = stats.parWithName(p).limits[0] # 68% limit
            mean = stats.parWithName(p).mean
            print(f"{p}: {mean:.5f} (Limits: {lims.lower:.5f} - {lims.upper:.5f})")
        except AttributeError:
            print(f"{p}: Not found or insufficient data")

    # Generate a triangle plot
    g = gdplt.getSubplotPlotter(width_inch=8)
    g.triangle_plot(samples, params, filled=True, title_limit=1)
    g.export('assets/tep_mcmc_recovered_triangle.png')
    print("\nTriangle plot exported to assets/tep_mcmc_recovered_triangle.png")
    
except Exception as e:
    print(f"Error loading samples: {e}")

