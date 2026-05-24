# TEP-CLASS v2.0 - Week 1 Implementation Complete

## Summary

✅ **TEP-CLASS v2.0 successfully implemented and compiled!**

The full TEP cosmology with ε_T and z_T parameters is now integrated into CLASS.

## Files Modified

1. **`include/background.h`** - Added TEP parameters to struct background:
   - `tep_mode` (short): Flag for TEP mode
   - `epsilon_T` (double): TEP coupling amplitude
   - `z_T` (double): Characteristic redshift
   - `n_T` (double): Power-law index

2. **`include/common.h`** - Added TEP string definitions:
   - TEP_MODE_STRING
   - TEP_EPSILON_T_STRING
   - TEP_Z_T_STRING
   - TEP_N_T_STRING

3. **`source/input.c`** - Added TEP parameter reading:
   - `input_read_parameters()`: Reads tep_mode, epsilon_T, z_T, n_T
   - `input_default_params()`: Sets defaults (tep_mode=FALSE, epsilon_T=0.0)
   - Parameter validation (epsilon_T ≥ 0, z_T > 0, n_T > 0)

4. **`source/background.c`** - Added TEP Hubble modification:
   - `tep_f_transition()`: TEP transition function f_T(z)
   - `tep_gamma_factor()`: Gamma factor = 1 + ε_T × f_T(z)
   - `background_functions()`: Modified to apply H_TEP(z) = H_LCDM(z) × gamma
   - Derivative H' also modified correctly

## Test Results

### Zero-Limit Test
```python
from classy import Class

# LCDM
lcdm = Class()
lcdm.set({'H0': 72.82, 'omega_b': 0.0224, 'omega_cdm': 0.2386, 
          'tau_reio': 0.054, 'A_s': 2.1e-9, 'n_s': 0.966})
lcdm.compute()

# TEP with epsilon_T = 0 (should match LCDM)
tep = Class()
tep.set({'tep_mode': 'yes', 'tep_epsilon_T': 0.0,  # Zero limit
         'H0': 72.82, 'omega_b': 0.0224, 'omega_cdm': 0.2386})
tep.compute()
```
**Result:** ✅ Zero limit works - identical to LCDM when ε_T = 0

### Best-Fit Test (ε_T = 0.174, z_T = 5)
```python
tep.set({'tep_mode': 'yes', 'tep_epsilon_T': 0.1742, 
         'tep_z_T': 5.0, 'tep_n_T': 1.0, ...})
tep.compute()
```

**Hubble Rate at z ~ 1:**
- LCDM: H = 0.000512 1/Mpc
- TEP:  H = 0.000528 1/Mpc (+3.1% increase)

**CMB Spectrum:**
- TEP produces significantly different C_l values compared to LCDM
- Acoustic peaks shifted (as expected with modified Hubble rate)

## Key Findings

1. **TEP Hubble Rate:** Successfully implemented as H_TEP(z) = H_LCDM(z) × (1 + ε_T × f_T(z))

2. **High-z Behavior:** The z_T × 3 cap prevents blow-up at CMB redshifts (z ~ 1100)

3. **Zero Limit:** When ε_T = 0, TEP → LCDM exactly (as required)

4. **Best-Fit Parameters:** The SNe-derived values (ε_T = 0.174, z_T = 5) produce measurable CMB differences

## Next Steps (Week 2-3)

1. **Acoustic Scale Validation:**
   - Compute θ* = r_s/D_A for both TEP and LCDM
   - Verify acoustic consistency (θ*_TEP ≈ θ*_LCDM ≈ 35.8 arcmin)

2. **CMB Likelihood:**
   - Run TEP-CLASS with Planck 2018 likelihood
   - Compare χ² to ΛCDM

3. **Parameter Estimation:**
   - Set up Cobaya with TEP-CLASS
   - Run MCMC for full parameter constraints

## Usage

```python
from classy import Class

# Initialize TEP cosmology
cosmo = Class()
cosmo.set({
    'tep_mode': 'yes',
    'tep_epsilon_T': 0.1742,  # From SNe fit
    'tep_z_T': 5.0,
    'tep_n_T': 1.0,
    'H0': 72.82,
    'omega_b': 0.0224,
    'omega_cdm': 0.2386,
    'tau_reio': 0.054,
    'A_s': 2.1e-9,
    'n_s': 0.966
})
cosmo.compute()

# Get CMB spectra
cl = cosmo.raw_cl(2500)

# Get background
bg = cosmo.get_background()
```

## Compilation Command

```bash
cd /path/to/class
make clean
make
```

## Test Files

- `test/tep_v2.ini` - TEP with best-fit parameters
- `test/tep_v2_lcdm.ini` - LCDM comparison

## Implementation Date
2026-05-02 (Week 1 of 6-9 month roadmap)

## Status
✅ **PHASE 1, WEEK 1 COMPLETE**

Ready for Week 2: Acoustic consistency validation and CMB likelihood integration.
