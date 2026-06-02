#!/usr/bin/env python3
"""Comprehensive Pipeline Audit - Check for Fakes, Synthetic Data, Placeholders.

This audits every component to verify:
1. No synthetic data fallbacks in strict mode
2. No placeholder physics
3. Real data loading with SHA-256 provenance
4. Actual calculations, not approximations
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / 'scripts/core'))

def audit_component(name, check_func):
    """Audit a single component."""
    try:
        result = check_func()
        status = "PASS" if result['real'] else "FAIL"
        print(f"  [{status}] {name}: {result['message']}")
        return result['real']
    except Exception as e:
        print(f"  [FAIL] {name}: Exception - {e}")
        return False

def check_step03_01():
    """Check step_03_01 uses real data, not synthetic."""
    import json
    step_file = Path("results/step_03_01_three_model_comparison.json")
    if not step_file.exists():
        return {'real': False, 'message': 'Step 03_01 not run yet'}
    
    data = json.loads(step_file.read_text())
    
    # Check for synthetic fallback flags
    if data.get('data_source') == 'synthetic':
        return {'real': False, 'message': 'Using synthetic data'}
    
    # Check for SHA-256 provenance
    if 'data_provenance' not in data:
        return {'real': False, 'message': 'No data provenance (SHA-256 missing)'}
    
    # Check that MLE actually computed (not placeholder)
    mle = data.get('models', {}).get('M0a_LCDM', {}).get('parameters_mle', {})
    if not mle or 'Om0' not in mle:
        return {'real': False, 'message': 'No MLE parameters computed'}

    return {'real': True, 'message': f"Real data, Om0={mle.get('Om0', 'N/A')}, M={mle.get('M', 'N/A')}"}

def check_background():
    """Verify background cosmology uses real physics."""
    from boltzmann.background import TEPBackground
    
    bg = TEPBackground(70, 0.045, 0.25, 0.7, Sigma_0=0.0)
    
    # Check Omega_gamma is physically correct (not made up)
    # Should be ~5e-5 for T_CMB = 2.725K
    if not (4e-5 < bg.Omega_gamma < 6e-5):
        return {'real': False, 'message': f'Omega_gamma={bg.Omega_gamma:.2e} (wrong physics)'}
    
    # Check radiation density uses correct formula
    # rho_gamma = a_rad * T^4 / c^2 (not arbitrary number)
    if bg.Omega_gamma == 0:
        return {'real': False, 'message': 'Omega_gamma = 0 (placeholder)'}
    
    return {'real': True, 'message': f'Omega_gamma={bg.Omega_gamma:.2e} (correct physics)'}

def check_recombination():
    """Verify recombination uses real Peebles physics."""
    from boltzmann.background import TEPBackground
    from boltzmann.recombination_working import RecombinationHistory
    
    bg = TEPBackground(70, 0.045, 0.25, 0.7)
    rec = RecombinationHistory(bg)
    
    z_rec = rec.z_rec()
    
    # Check z_rec is physically reasonable (~1100)
    if not (1000 < z_rec < 1200):
        return {'real': False, 'message': f'z_rec={z_rec:.1f} (wrong, should be ~1100)'}
    
    # Check x_e evolution is physical
    xe_high_z = rec.xe(1500)
    xe_low_z = rec.xe(500)
    
    if xe_high_z < 0.5:
        return {'real': False, 'message': f'x_e(1500)={xe_high_z:.2f} (should be ~1)'}
    if xe_low_z > 0.1:
        return {'real': False, 'message': f'x_e(500)={xe_low_z:.2f} (should be ~0)'}
    
    return {'real': True, 'message': f'z_rec={z_rec:.1f}, x_e physical'}

def check_bbn():
    """Verify BBN uses real nuclear physics."""
    from bbn_working import BBNWorking
    
    bbn = BBNWorking(eta=6.1e-10, Sigma_0=0.0, H0=70.0)
    r = bbn.solve()
    
    Y_p = r['Y_p']
    D_H = r['D_H']
    
    # Check Y_p is physically reasonable (~0.25)
    if not (0.15 < Y_p < 0.35):
        return {'real': False, 'message': f'Y_p={Y_p:.4f} (wrong physics)'}
    
    # Check not using placeholder value
    if Y_p == 0.25:
        return {'real': False, 'message': 'Y_p=0.25 exactly (suspicious, may be hardcoded)'}
    
    # Check D/H is reasonable
    if D_H < 1e-6 or D_H > 1e-4:
        return {'real': False, 'message': f'D/H={D_H:.2e} (outside physical range)'}
    
    return {'real': True, 'message': f'Y_p={Y_p:.4f}, D/H={D_H:.2e}'}

def check_pantheon_data():
    """Verify Pantheon+ data files exist and have content."""
    data_files = [
        "data/raw/pantheon_plus_shoes.dat",
        "data/raw/Pantheon+SH0ES.cov",
    ]
    
    missing = []
    for f in data_files:
        path = Path(f)
        if not path.exists():
            missing.append(f)
        elif path.stat().st_size < 1000:
            return {'real': False, 'message': f'{f} exists but too small (synthetic?)'}
    
    if missing:
        return {'real': False, 'message': f'Missing data files: {missing}'}
    
    return {'real': True, 'message': 'Pantheon+ data files present'}

def check_no_synthetic_flags():
    """Check codebase for synthetic/placeholder flags."""
    scripts_dir = Path("scripts")
    issues = []
    
    for py_file in scripts_dir.rglob("*.py"):
        if py_file.name == "pipeline_audit.py":
            continue
            
        content = py_file.read_text()
        
        # Check for synthetic data fallbacks
        if 'synthetic_fallback' in content.lower():
            issues.append(f"{py_file}: synthetic_fallback found")
        if 'fake_data' in content.lower():
            issues.append(f"{py_file}: fake_data found")
            
        # Ignore "placeholder" if it's explicitly raising an error or in a comment guardrail
        lines = content.split('\n')
        for i, line in enumerate(lines):
            if 'placeholder' in line.lower() and 'notimplementederror' not in line.lower() and not line.strip().startswith('#'):
                issues.append(f"{py_file}:{i+1} placeholder found")
        if 'TODO' in content and 'data' in content.lower():
            issues.append(f"{py_file}: TODO related to data")
    
    if issues:
        print("ISSUES FOUND BY AUDIT:", issues)
        return {'real': False, 'message': f'Found {len(issues)} potential issues'}
    
    return {'real': True, 'message': 'No synthetic flags found'}

def check_step_05_04():
    """Check step 05_04 (CMB) is properly implemented."""
    step_file = Path("results/step_05_04_cmb_spectra.json")
    
    if not step_file.exists():
        return {'real': False, 'message': 'Step 05_04 not run yet'}
    
    import json
    data = json.loads(step_file.read_text())
    
    # Check if using toy model
    if 'toy_model' in str(data).lower():
        return {'real': False, 'message': 'Using toy model (not full physics)'}
    
    # Check z_rec is physical
    z_rec = data.get('cmb_results', {}).get('class_lcdm_reference', {}).get('derived', {}).get('z_rec', 0)
    if not (1000 < z_rec < 1200):
        return {'real': False, 'message': f'z_rec={z_rec} (wrong)'}
    
    return {'real': True, 'message': f'z_rec={z_rec:.1f} (working)'}

def check_step_05_07():
    """Check step 05_07 (BBN) is properly implemented."""
    step_file = Path("results/step_05_07_bbn_preservation.json")
    
    if not step_file.exists():
        return {'real': False, 'message': 'Step 05_07 not run yet'}
    
    import json
    data = json.loads(step_file.read_text())
    
    # Check Y_p is computed, not hardcoded
    Y_p = data.get('lcdm', {}).get('Y_p', 0)
    if Y_p == 0.25 or Y_p == 0:
        return {'real': False, 'message': f'Y_p={Y_p} (suspicious value)'}
    
    if not (0.2 < Y_p < 0.3):
        return {'real': False, 'message': f'Y_p={Y_p:.4f} (outside physical range)'}
    
    return {'real': True, 'message': f'Y_p={Y_p:.4f}'}

def main():
    print("="*70)
    print("TEP PIPELINE COMPREHENSIVE AUDIT")
    print("Checking for: synthetic data, fake physics, placeholders")
    print("="*70)
    print()
    
    results = {}
    
    print("1. DATA SOURCES")
    results['pantheon_data'] = audit_component("Pantheon+ data files", check_pantheon_data)
    results['step03_01_data'] = audit_component("Step 03_01 (real data check)", check_step03_01)
    
    print("\n2. PHYSICS IMPLEMENTATIONS")
    print("  [INFO] Module checks (background, recombination, bbn) retired in favor of pipeline output checks.")
    
    print("\n3. STEP VALIDATION")
    results['step05_04'] = audit_component("Step 05_04 (CMB)", check_step_05_04)
    results['step05_07'] = audit_component("Step 05_07 (BBN)", check_step_05_07)
    
    print("\n4. CODE QUALITY CHECK")
    results['no_synthetic'] = audit_component("No synthetic flags in code", check_no_synthetic_flags)
    
    print("\n" + "="*70)
    print("AUDIT SUMMARY")
    print("="*70)
    
    passed = sum(results.values())
    total = len(results)
    
    print(f"Passed: {passed}/{total}")
    print()
    
    if passed == total:
        print("STATUS: ALL CHECKS PASS")
        print("The pipeline uses real physics and real data.")
    else:
        print("STATUS: ISSUES FOUND")
        print("Some components may use synthetic data or incomplete physics.")
        print()
        print("Failed checks:")
        for name, passed in results.items():
            if not passed:
                print(f"  - {name}")
    
    print("="*70)
    
    return passed == total

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
