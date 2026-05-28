import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np
from scipy import stats
from core.tep_cosmology import TEPCosmology
from c0_common import read_json, RAW_DIR
import os
import pandas as pd

H0 = 70.0
epsilon_T = 0.425647
z_T = 1.0

# Load data
sn_path = os.path.join(RAW_DIR, 'supernovae', 'pantheon_plus', 'Pantheon+SH0ES.dat')
df = pd.read_csv(sn_path, sep='\s+')
z = df['zHD'].values
mb = df['mB'].values
host_mass = df['HOST_LOGMASS'].values
valid = (z > 0.01) & (host_mass > 0)
z = z[valid]
mb = mb[valid]
host_mass = host_mass[valid] > 10.0

cosmo_tep = TEPCosmology(H0=H0, Omega_m=0.3, epsilon_T=epsilon_T, z_T=z_T)
mu_tep = np.array([5.0 * np.log10(cosmo_tep.luminosity_distance(zi)) + 25.0 for zi in z])
R_tep = mb - mu_tep

slope, intercept, r_value, p_value, std_err = stats.linregress(host_mass, R_tep)
print(f"Corrected rho: {r_value:.3f}, p: {p_value:.2e}")
