# TEP-CLASS v2.0 Patch Instructions

## Overview
This patch adds TEP support to CLASS using ε_T and z_T parameters from SNe fit.

## Files to Modify

### 1. include/background.h

Add to `struct background` (after line ~300):

```c
  /* TEP parameters */
  short tep_mode;           /**< Flag for TEP mode (_TRUE_ or _FALSE_) */
  double epsilon_T;       /**< TEP coupling amplitude */
  double z_T;              /**< Characteristic redshift for TEP transition */
  double n_T;              /**< Power-law index for TEP transition */
```

### 2. include/common.h

Add TEP strings for input parsing (after other #define strings):

```c
/* TEP parameters */
#define TEP_MODE_STRING "tep_mode"
#define TEP_EPSILON_T_STRING "tep_epsilon_T"
#define TEP_Z_T_STRING "tep_z_T"
#define TEP_N_T_STRING "tep_n_T"
```

### 3. source/input.c

#### 3a. Add to `input_read_parameters` function

After reading other background parameters (around line ~500):

```c
  /* TEP parameters */
  class_read_string(TEP_MODE_STRING, string1);
  if (strcmp(string1, "yes") == 0) {
    pba->tep_mode = _TRUE_;
  } else {
    pba->tep_mode = _FALSE_;
  }
  
  if (pba->tep_mode == _TRUE_) {
    class_read_double(TEP_EPSILON_T_STRING, pba->epsilon_T);
    class_read_double(TEP_Z_T_STRING, pba->z_T);
    class_read_double(TEP_N_T_STRING, pba->n_T);
    
    /* Validate TEP parameters */
    class_test(pba->epsilon_T < 0, 
               errmsg, "epsilon_T must be non-negative");
    class_test(pba->z_T <= 0, 
               errmsg, "z_T must be positive");
    class_test(pba->n_T <= 0, 
               errmsg, "n_T must be positive");
  } else {
    /* Set defaults when TEP is disabled */
    pba->epsilon_T = 0.0;
    pba->z_T = 1.0;
    pba->n_T = 1.0;
  }
```

#### 3b. Add to `input_default_params` function

```c
  pba->tep_mode = _FALSE_;
  pba->epsilon_T = 0.0;
  pba->z_T = 1.0;
  pba->n_T = 1.0;
```

### 4. source/background.c

#### 4a. Add TEP function declarations at top of file

```c
/* TEP helper functions */
double tep_f_transition(struct background *pba, double z) {
    if (pba->tep_mode == _FALSE_ || pba->epsilon_T == 0.0) {
        return 0.0;
    }
    
    /* Cap z contribution to prevent blow-up at very high z */
    double z_effective = z;
    if (pba->z_T > 0 && z > pba->z_T * 3.0) {
        z_effective = pba->z_T * 3.0;
    }
    
    double ratio = z_effective / pba->z_T;
    double exponent = pow(ratio, pba->n_T);
    
    return 1.0 - exp(-exponent);
}

double tep_gamma_factor(struct background *pba, double z) {
    if (pba->tep_mode == _FALSE_ || pba->epsilon_T == 0.0) {
        return 1.0;
    }
    double f_T = tep_f_transition(pba, z);
    return 1.0 + pba->epsilon_T * f_T;
}
```

#### 4b. Modify `background_functions` 

Find where H is computed (around line ~2000), add TEP modification:

```c
  /* Original Hubble rate computation */
  pvecback[pba->index_bg_H] = ... /* existing code */
  
  /* TEP modification if enabled */
  if (pba->tep_mode == _TRUE_ && pba->epsilon_T != 0.0) {
    double z = 1.0/a - 1.0;
    double H_original = pvecback[pba->index_bg_H];
    double gamma = tep_gamma_factor(pba, z);
    pvecback[pba->index_bg_H] = H_original * gamma;
    
    /* Recompute conformal Hubble rate */
    pvecback[pba->index_bg_H * pba->index_bg_a] = a * pvecback[pba->index_bg_H];
    
    /* Recompute any other derived quantities */
    /* ... */
  }
```

#### 4c. Modify `background_initial_conditions` or add TEP init

In `background_init`, after reading parameters:

```c
  /* TEP initialization */
  if (pba->tep_mode == _TRUE_) {
    fprintf(stdout, "TEP mode enabled:\n");
    fprintf(stdout, "  epsilon_T = %.6f\n", pba->epsilon_T);
    fprintf(stdout, "  z_T = %.6f\n", pba->z_T);
    fprintf(stdout, "  n_T = %.6f\n", pba->n_T);
    fprintf(stdout, "  f_T(z=1) = %.6f\n", tep_f_transition(pba, 1.0));
    fprintf(stdout, "  f_T(z=1100) = %.6f\n", tep_f_transition(pba, 1100.0));
  }
```

### 5. test/tep.ini (Create new test file)

```ini
# TEP test parameter file
# Uses best-fit parameters from Pantheon+ SNe

root = output/tep_test_

# TEP parameters
tep_mode = yes
tep_epsilon_T = 0.1742
tep_z_T = 5.0
tep_n_T = 1.0

# Background parameters
H0 = 72.82
omega_b = 0.0224
omega_cdm = 0.2386

# ... rest of standard CLASS parameters
```

### 6. python/classy.pyx (if using Python wrapper)

Add TEP parameters to the classy wrapper:

```python
# In the set() method, add TEP parameters
if 'tep_mode' in k:
    self.ba.tep_mode = 1 if v == 'yes' else 0
if 'tep_epsilon_T' in k:
    self.ba.epsilon_T = float(v)
if 'tep_z_T' in k:
    self.ba.z_T = float(v)
if 'tep_n_T' in k:
    self.ba.n_T = float(v)
```

## Compilation Instructions

```bash
cd /path/to/class_public

# Apply the patch (or make edits manually)
# Then compile:
make clean
make

# Test with TEP:
./class test/tep.ini

# Compare with LCDM:
./class test/lcdm.ini
```

## Validation Tests

### Test 1: Zero Limit
```ini
tep_mode = yes
tep_epsilon_T = 0.0
```
Should give identical results to LCDM.

### Test 2: SNe Best-Fit
```ini
tep_mode = yes
tep_epsilon_T = 0.1742
tep_z_T = 5.0
```
Should run successfully and produce modified CMB spectrum.

### Test 3: Acoustic Consistency
Check that θ* is preserved despite 31% sound horizon shift.

## Expected Timeline
- Week 1: Files 1-3 (header and input modifications)
- Week 2: Files 4-6 (background.c and testing)
- Week 3: Validation and debugging
- Week 4: CMB spectrum production

## Notes
- The z_T * 3 cap is a temporary fix for high-z behavior
- Full TEP theory might require different high-z behavior
- This parametrization is optimized for z < 2 (SNe range)
- CMB acoustic peaks may need additional tuning
