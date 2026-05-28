import numpy as np
import pandas as pd
import glob

dfs = []
for chain_file in glob.glob('results/outputs/tep_cobaya_sne.*.txt'):
    with open(chain_file, 'r') as f:
        lines = f.readlines()
    if len(lines) < 2: continue
    header = lines[0].strip().lstrip('#').split()
    data = []
    for line in lines[1:]:
        data.append([float(x) for x in line.strip().split()])
    
    df = pd.DataFrame(data, columns=header)
    burnin = int(len(df) * 0.3)
    dfs.append(df.iloc[burnin:])

df = pd.concat(dfs, ignore_index=True)

print("Current MCMC Posteriors (Combined 4 Chains, Post-Burnin):")
params_of_interest = ['tep_epsilon_T', 'tep_z_T', 'H0', 'omega_cdm', 'A_s', 'n_s']

for param in params_of_interest:
    if param in df.columns:
        mean = df[param].mean()
        std = df[param].std()
        print(f"{param}: {mean:.5f} +/- {std:.5f}")

print("\nChain statistics:")
print(f"Total steps (post-burnin): {len(df)}")
