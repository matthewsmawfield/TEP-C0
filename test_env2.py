import numpy as np
import pandas as pd
from scipy import stats
from core.tep_cosmology import TEPCosmology
from core.cosmology import CosmologyFLRW
import os

H0 = 70.0
epsilon_T = 0.425647
z_T = 1.0

sn_path = 'data/raw/supernovae/pantheon_plus/Pantheon+SH0ES.dat'
df = pd.read_csv(sn_path, sep='\s+')
z = df['zHD'].values
mb = df['mB'].values
host_mass = df['HOST_LOGMASS'].values
valid = (z > 0.01) & (host_mass > 0)
z = z[valid]
mb = mb[valid]
host_mass = host_mass[valid]

# Calculate LCDM residuals
cosmo_lcdm = CosmologyFLRW(H0=H0, Om0=0.3, Ok0=0.0)
mu_lcdm = cosmo_lcdm.distance_modulus(z)
R_lcdm = mb - mu_lcdm

# Calculate global TEP residuals (NO explicit host mass injection)
cosmo_tep = TEPCosmology(H0=H0, Omega_m=0.3, epsilon_T=epsilon_T, z_T=z_T)
mu_tep = np.array([5.0 * np.log10(cosmo_tep.luminosity_distance(zi)) + 25.0 for zi in z])
R_tep = mb - mu_tep

# Test correlation with host mass
r_lcdm, p_lcdm = stats.pearsonr(host_mass, R_lcdm)
r_tep, p_tep = stats.pearsonr(host_mass, R_tep)

print(f"LCDM rho(R, mass): {r_lcdm:.3f}, p: {p_lcdm:.2e}")
print(f"TEP rho(R, mass): {r_tep:.3f}, p: {p_tep:.2e}")
