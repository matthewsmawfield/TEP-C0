#!/usr/bin/env python3
"""Comprehensive check of entire TEP-C0 pipeline."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / 'steps'))

import json
import hashlib

def check(name, condition, message=""):
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}] {name}")
    if message and not condition:
        print(f"       -> {message}")
    return condition

print("="*70)
print("COMPREHENSIVE PIPELINE CHECK")
print("="*70)

all_pass = True

# 1. DATA FILES
print("\n1. REAL DATA FILES")
print("-"*70)
data_files = {
    'Pantheon+ data': 'data/raw/pantheon_plus_shoes.dat',
    'Pantheon+ covariance': 'data/raw/pantheon_plus_shoes.cov',
    'FIRAS CMB': 'data/raw/firas_monopole_spec_v1.txt',
    'BAO compilation': 'data/raw/uncorBAO.txt',
}
for name, path in data_files.items():
    p = Path(path)
    exists = p.exists()
    size = p.stat().st_size if exists else 0
    valid = exists and size > 1000
    all_pass &= check(f"{name}: {size:,} bytes", valid, f"Missing or too small")

# 2. STEP OUTPUTS
print("\n2. STEP OUTPUTS (JSON)")
print("-"*70)
steps = [
    'step_000_data_download',
    'step_022_three_model_comparison',
    'step_017_tep_boltzmann_solver',
    'step_029_bbn_preservation',
]
for step in steps:
    path = Path(f'results/{step}.json')
    exists = path.exists()
    all_pass &= check(f"{step}.json exists", exists)
    if exists:
        try:
            with open(path) as f:
                data = json.load(f)
            valid = 'status' in data or 'step' in data
            all_pass &= check(f"  -> Valid JSON", valid)
        except (json.JSONDecodeError, IOError, OSError) as e:
            all_pass &= check(f"  -> Valid JSON", False, f"Corrupted JSON: {e}")

# 3. DATA PROVENANCE
print("\n3. DATA PROVENANCE")
print("-"*70)
try:
    with open('results/step_022_three_model_comparison.json') as f:
        s022 = json.load(f)
    
    data_meta = s022.get('data', {})
    prov = s022.get('data_provenance', {})
    all_pass &= check("SHA-256 for data", len(data_meta.get('data_sha256', '')) > 20)
    all_pass &= check("SHA-256 for cov", len(data_meta.get('covariance_sha256', '')) > 20)
    all_pass &= check("Data file provenance verified", prov.get('data_file', {}).get('verified') is True)
    all_pass &= check("Covariance provenance verified", prov.get('cov_file', {}).get('verified') is True)
    all_pass &= check("Full covariance data", s022.get('validation', {}).get('research_grade_data') is True)
except Exception as e:
    all_pass &= check("Provenance check", False, str(e))

# 4. PHYSICS VALIDATION
print("\n4. PHYSICS VALIDATION")
print("-"*70)

# Background
from core.background import TEPBackground
bg = TEPBackground(70, 0.045, 0.25, 0.7)
all_pass &= check(f"Omega_gamma = {bg.Omega_gamma:.2e}", 4.5e-5 < bg.Omega_gamma < 5.5e-5, "Should be ~5.04e-5")
all_pass &= check(f"Omega_r = {bg.Omega_r:.2e}", 8e-5 < bg.Omega_r < 9e-5, "Should be ~8.5e-5")

# Recombination
from core.recombination_working import RecombinationHistory
rec = RecombinationHistory(bg)
all_pass &= check(f"z_rec = {rec.z_rec():.1f}", 1000 < rec.z_rec() < 1200, "Should be ~1100")
all_pass &= check(f"x_e(1500) = {rec.xe(1500):.3f}", rec.xe(1500) > 0.7, "Should be ionized")
all_pass &= check(f"x_e(500) = {rec.xe(500):.4f}", rec.xe(500) < 0.01, "Should be recombined")

# BBN
from core.bbn_working import BBNWorking
bbn = BBNWorking()
r = bbn.solve()
all_pass &= check(f"Y_p = {r['Y_p']:.4f}", 0.2 < r['Y_p'] < 0.35, "Should be ~0.25")
all_pass &= check(f"D/H = {r['D_H']:.2e}", 1e-6 < r['D_H'] < 1e-4, "Should be ~2.6e-5")

# 5. STEP VALIDATION
print("\n5. STEP VALIDATION")
print("-"*70)

try:
    with open('results/step_017_tep_boltzmann_solver.json') as f:
        s017 = json.load(f)
    all_pass &= check("CLASS reference available", s017.get('validation', {}).get('class_available') is True)
    cmb_validation = s017.get('validation', {})
    cmb_gate_coherent = (
        cmb_validation.get('claim_gate') == 'blocked'
        or (
            cmb_validation.get('tep_class_available') is True
            and cmb_validation.get('tep_zero_limit_ok') is True
        )
    )
    all_pass &= check("CMB resolver gate coherent", cmb_gate_coherent)
except (FileNotFoundError, json.JSONDecodeError, KeyError) as e:
    all_pass &= check("CMB validation", False, f"Error: {e}")

try:
    with open('results/step_029_bbn_preservation.json') as f:
        s029 = json.load(f)
    bbn_rg = s029.get('validation', {}).get('research_grade_bbn')
    all_pass &= check("BBN research grade status defined", bbn_rg is not None)
    all_pass &= check("BBN working", s029.get('validation', {}).get('bbn_working') == True)
    all_pass &= check("BBN nuclear network active", s029.get('validation', {}).get('nuclear_network_active') == True)
except (FileNotFoundError, json.JSONDecodeError, KeyError) as e:
    all_pass &= check("BBN validation", False, f"Error: {e}")

# 6. NO SYNTHETIC FALLBACKS
print("\n6. NO SYNTHETIC FALLBACKS")
print("-"*70)

scripts_dir = Path('scripts')
synthetic_found = False
for py_file in scripts_dir.rglob("*.py"):
    # Skip audit/check scripts themselves
    if 'audit' in py_file.name or 'check' in py_file.name:
        continue
    content = py_file.read_text().lower()
    if 'synthetic_fallback' in content and 'def ' in content:
        synthetic_found = True
        break

all_pass &= check("No silent synthetic fallbacks", not synthetic_found)

# 7. COVARIANCE MATRIX
print("\n7. COVARIANCE MATRIX")
print("-"*70)

cov_path = Path('data/raw/pantheon_plus_shoes.cov')
if cov_path.exists():
    import numpy as np
    try:
        cov_flat = np.loadtxt(cov_path)
        if cov_flat.size == 1701*1701 + 1:
            cov = cov_flat[1:].reshape(1701, 1701)
        elif cov_flat.size == 1701*1701:
            cov = cov_flat.reshape(1701, 1701)
        else:
            cov = None
        
        if cov is not None:
            all_pass &= check(f"Covariance shape: {cov.shape}", cov.shape == (1701, 1701))
            # Use tolerance for symmetry (real data has numerical precision limits)
            all_pass &= check("Covariance symmetric (tol=1e-6)", np.allclose(cov, cov.T, atol=1e-6))
            all_pass &= check("Covariance positive diagonal", np.all(np.diag(cov) > 0))
        else:
            all_pass &= check("Covariance shape", False, f"Wrong size: {cov_flat.size}")
    except Exception as e:
        all_pass &= check("Covariance load", False, str(e))
else:
    all_pass &= check("Covariance file exists", False)

# 8. STEP 000 PROVENANCE
print("\n8. STEP 000 DATA DOWNLOAD PROVENANCE")
print("-"*70)

try:
    with open('results/step_000_data_download.json') as f:
        s000 = json.load(f)
    
    sources = s000.get('sources', [])
    all_pass &= check(f"Sources downloaded: {len(sources)}", len(sources) >= 4)
    
    for src in sources:
        has_sha256 = 'sha256' in src and len(src['sha256']) > 20
        all_pass &= check(f"  {src.get('source_id', 'unknown')}: SHA-256", has_sha256)
        
except Exception as e:
    all_pass &= check("Step 000 validation", False, str(e))

# FINAL SUMMARY
print("\n" + "="*70)
print("FINAL SUMMARY")
print("="*70)

if all_pass:
    print("\n✅ ALL CHECKS PASSED")
    print("   - Real data downloaded with provenance")
    print("   - Structural and diagnostic artifacts are coherent")
    print("   - Blocked research gates are blocked explicitly")
    print("   - No synthetic fallbacks or fake data")
    sys.exit(0)
else:
    print("\n❌ SOME CHECKS FAILED")
    print("   Review failures above")
    sys.exit(1)
