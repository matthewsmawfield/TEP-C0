#!/usr/bin/env python3
"""
TEP-C0 Step 09: Pantheon+ Cosmological Analysis (REFINED GAMMA)
===============================================================
REFINED gamma formulation matching CLASS patch:
gamma(z) = 1 + epsilon_T * log(1+z) * f_T(z)
This gives gamma > 1, H_TEP < H_LCDM, larger distances.

Data: Pantheon+SH0ES.dat (1701 SNe Ia)
"""

import sys, json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.integrate import quad

PROJECT_ROOT = Path(__file__).resolve().parents[2]

def parse_pantheon_plus(filepath):
    """Parse Pantheon+ SALT2 parameters."""
    data = []
    with open(filepath) as f:
        header = f.readline().strip().split()
        for line in f:
            if not line.strip():
                continue
            p = line.strip().split()
            if len(p) >= 21:
                z = float(p[2])
                mB = float(p[19])
                mB_err = float(p[20])
                c = float(p[15]) if p[15] != 'NA' else 0
                x1 = float(p[17]) if p[17] != 'NA' else 0
                data.append({'name': p[0], 'z': z, 'mB': mB, 'mB_err': mB_err,
                             'c': c, 'x1': x1})
    return data

def gamma_tep(z, eps_T, z_T, n_T):
    """REFINED gamma: matches CLASS implementation exactly."""
    if z_T <= 0 or n_T <= 0:
        return 1.0
    f_T = 1.0 / (1.0 + (z / z_T)**n_T)
    return 1.0 + eps_T * np.log(1.0 + z) * f_T

def h_lcdm(z, Om, OL, H0):
    return H0 * np.sqrt(Om * (1+z)**3 + OL)

def h_tep_lcdm(z, Om, OL, eps_T, z_T, n_T, H0):
    """TEP-C0: LCDM H divided by gamma (REFINED + sign)."""
    return h_lcdm(z, Om, OL, H0) / gamma_tep(z, eps_T, z_T, n_T)

def h_tep_matter_only(z, Om, eps_T, z_T, n_T, H0):
    """Pure TEP: matter-only H / gamma (NO dark energy)."""
    return H0 * np.sqrt(Om * (1+z)**3) / gamma_tep(z, eps_T, z_T, n_T)

def d_L(z, H_func, args):
    c = 299792.458
    integral, _ = quad(lambda zp: c / H_func(zp, *args), 0, z, limit=100)
    return (1 + z) * integral

def main():
    print("=" * 60)
    print("TEP-C0 Step 09: Pantheon+ Cosmological Analysis")
    print("REFINED gamma: 1 + eps_T * log(1+z) * f_T(z)")
    print("=" * 60)
    
    data_file = PROJECT_ROOT / "data" / "pantheon_plus" / "Pantheon+SH0ES.dat"
    if not data_file.exists():
        print("[ERROR] No Pantheon+ data found")
        return 1
    
    print(f"\nLoading {data_file.name}...")
    data = parse_pantheon_plus(data_file)
    print(f"  {len(data)} SNe parsed, z=[{min(d['z'] for d in data):.3f}, {max(d['z'] for d in data):.3f}]")
    
    z = np.array([d['z'] for d in data])
    mB = np.array([d['mB'] for d in data])
    mB_err = np.array([d['mB_err'] for d in data])
    c_arr = np.array([d['c'] for d in data])
    x1_arr = np.array([d['x1'] for d in data])
    
    # Standardization: mu = mB - M_B + alpha*x1 - beta*c
    # We fit M_B, alpha, beta as nuisance parameters
    
    def mu_from_dL(z_val, H_func, args):
        d = d_L(z_val, H_func, args)
        return 5 * np.log10(d) + 25
    
    def chi2_lcdm(p):
        Om, OL, M_B, alpha, beta = p
        if Om < 0.05 or Om > 1.0 or OL < 0 or OL > 1.5 or (Om+OL) > 1.5:
            return 1e10
        args = (Om, OL, 70.0)
        mu_model = np.array([mu_from_dL(zi, h_lcdm, args) for zi in z])
        mB_model = mu_model + M_B - alpha * x1_arr + beta * c_arr
        return np.sum(((mB - mB_model) / mB_err)**2)
    
    def chi2_tep_c0(p):
        Om, OL, eps_T, z_T, n_T, M_B, alpha, beta = p
        if Om < 0.05 or Om > 1.0 or OL < 0 or OL > 1.5 or (Om+OL) > 1.5:
            return 1e10
        if eps_T < 0 or eps_T > 10 or z_T < 0.01 or z_T > 10 or n_T < 0.1 or n_T > 10:
            return 1e10
        args = (Om, OL, eps_T, z_T, n_T, 70.0)
        mu_model = np.array([mu_from_dL(zi, h_tep_lcdm, args) for zi in z])
        mB_model = mu_model + M_B - alpha * x1_arr + beta * c_arr
        return np.sum(((mB - mB_model) / mB_err)**2)
    
    def chi2_tep_pure(p):
        Om, eps_T, z_T, n_T, M_B, alpha, beta = p
        if Om < 0.05 or Om > 1.0:
            return 1e10
        if eps_T < 0 or eps_T > 20 or z_T < 0.01 or z_T > 10 or n_T < 0.1 or n_T > 10:
            return 1e10
        args = (Om, eps_T, z_T, n_T, 70.0)
        mu_model = np.array([mu_from_dL(zi, h_tep_matter_only, args) for zi in z])
        mB_model = mu_model + M_B - alpha * x1_arr + beta * c_arr
        return np.sum(((mB - mB_model) / mB_err)**2)
    
    # Fit LCDM
    print("\n--- LCDM Baseline ---")
    p0_lcdm = [0.3, 0.7, -19.3, 0.14, 3.1]
    bounds_lcdm = [(0.05, 1.0), (0.0, 1.5), (-21.0, -17.0), (0.0, 0.5), (0.0, 5.0)]
    r_l = minimize(chi2_lcdm, p0_lcdm, method='L-BFGS-B', bounds=bounds_lcdm, options={'maxiter': 1000})
    Om_l, OL_l, MB_l, a_l, b_l = r_l.x
    print(f"  LCDM: Om={Om_l:.4f}, OL={OL_l:.4f}, M_B={MB_l:.4f}, alpha={a_l:.4f}, beta={b_l:.4f}")
    print(f"        chi2={r_l.fun:.2f}")
    
    # Fit TEP-C0 (LCDM + gamma)
    print("\n--- TEP-C0 (LCDM + gamma) ---")
    p0_tep = [0.3, 0.7, 0.5, 1.0, 1.0, -19.3, 0.14, 3.1]
    bounds_tep = [(0.05, 1.0), (0.0, 1.5), (0.0, 10.0), (0.01, 10.0), (0.1, 5.0),
                  (-21.0, -17.0), (0.0, 0.5), (0.0, 5.0)]
    r_t = minimize(chi2_tep_c0, p0_tep, method='L-BFGS-B', bounds=bounds_tep, options={'maxiter': 1000})
    Om_t, OL_t, eps_t, zt_t, nt_t, MB_t, a_t, b_t = r_t.x
    print(f"  TEP-C0: Om={Om_t:.4f}, OL={OL_t:.4f}, eps={eps_t:.4f}, z_T={zt_t:.4f}, n_T={nt_t:.4f}")
    print(f"          M_B={MB_t:.4f}, alpha={a_t:.4f}, beta={b_t:.4f}")
    print(f"          chi2={r_t.fun:.2f}")
    
    # Fit TEP pure (matter-only + gamma)
    print("\n--- TEP Pure (matter + gamma) ---")
    p0_pure = [0.3, 1.0, 0.5, 2.0, -19.3, 0.14, 3.1]
    bounds_pure = [(0.05, 1.0), (0.0, 20.0), (0.01, 10.0), (0.1, 10.0),
                   (-21.0, -17.0), (0.0, 0.5), (0.0, 5.0)]
    r_p = minimize(chi2_tep_pure, p0_pure, method='L-BFGS-B', bounds=bounds_pure, options={'maxiter': 1000})
    Om_p, eps_p, zt_p, nt_p, MB_p, a_p, b_p = r_p.x
    print(f"  TEP pure: Om={Om_p:.4f}, eps={eps_p:.4f}, z_T={zt_p:.4f}, n_T={nt_p:.4f}")
    print(f"            M_B={MB_p:.4f}, alpha={a_p:.4f}, beta={b_p:.4f}")
    print(f"            chi2={r_p.fun:.2f}")
    
    # Model comparison
    n = len(data)
    bic_l = r_l.fun + 5 * np.log(n)
    bic_t = r_t.fun + 8 * np.log(n)
    bic_p = r_p.fun + 7 * np.log(n)
    
    print(f"\n--- Model Comparison (BIC) ---")
    print(f"  LCDM:     chi2={r_l.fun:.2f}, BIC={bic_l:.2f} (5 params)")
    print(f"  TEP-C0:   chi2={r_t.fun:.2f}, BIC={bic_t:.2f} (8 params)")
    print(f"  TEP pure: chi2={r_p.fun:.2f}, BIC={bic_p:.2f} (7 params)")
    
    print(f"\n  TEP-C0 gamma(z): ", end="")
    for z_test in [0.01, 0.1, 0.5, 1.0, 2.0]:
        print(f"z={z_test:.2f}->{gamma_tep(z_test, eps_t, zt_t, nt_t):.4f} ", end="")
    print()
    
    print(f"  TEP pure gamma(z): ", end="")
    for z_test in [0.01, 0.1, 0.5, 1.0, 2.0]:
        print(f"z={z_test:.2f}->{gamma_tep(z_test, eps_p, zt_p, nt_p):.4f} ", end="")
    print()
    
    # Save results
    results_dir = PROJECT_ROOT / "results"
    results_dir.mkdir(exist_ok=True)
    results = {
        "step": "09_pantheon_analysis",
        "status": "COMPLETE",
        "gamma_formulation": "REFINED: gamma = 1 + eps_T * log(1+z) * f_T(z)",
        "data": f"Pantheon+SH0ES.dat ({n} SNe)",
        "models": {
            "lcdm": {"Om": Om_l, "OL": OL_l, "M_B": MB_l, "alpha": a_l, "beta": b_l,
                     "chi2": r_l.fun, "bic": bic_l, "converged": r_l.success},
            "tep_c0": {"Om": Om_t, "OL": OL_t, "eps": eps_t, "z_T": zt_t, "n_T": nt_t,
                       "M_B": MB_t, "alpha": a_t, "beta": b_t,
                       "chi2": r_t.fun, "bic": bic_t, "converged": r_t.success},
            "tep_pure": {"Om": Om_p, "eps": eps_p, "z_T": zt_p, "n_T": nt_p,
                         "M_B": MB_p, "alpha": a_p, "beta": b_p,
                         "chi2": r_p.fun, "bic": bic_p, "converged": r_p.success},
        }
    }
    with open(results_dir / "step09_pantheon_results.json", "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\nResults saved to {results_dir / 'step09_pantheon_results.json'}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
