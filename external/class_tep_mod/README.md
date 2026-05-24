# CLASS TEP Modification

Research-grade TEP implementation in CLASS (Cosmic Linear Anisotropy Solving System).

## Modification Strategy

The TEP modifies the Hubble rate as:
```
H_TEP(z) = H_LCDM(z) × Gamma_TEP(z)
```

where:
```
Gamma_TEP(z) = exp(Sigma_0 × c/H0 × ln(1+z))
```

This modification must be applied in `background.c` where the Hubble rate is computed.

## Files Modified

- `source/background.c`: Add TEP Hubble modification
- `include/common.h`: Add TEP parameters to background structure
- `source/input.c`: Parse TEP parameters from input file

## Implementation Status

- [x] Design specification
- [x] `background.c` background-expansion hook in `external/class`
- [x] `background.h` TEP parameter storage in `external/class`
- [x] `input.c` parsing for `tep_mode`, `tep_Sigma_0`, and `tep_epsilon_T`
- [x] Python resolver detects whether the active `classy` build exposes TEP parameters
- [x] Rebuilt local `classy` runtime is available at `external/class_tep_mod/python_runtime`
- [x] Validated the `tep_Sigma_0=0` limit against vanilla CLASS spectra in Step 017/028
- [ ] Add perturbation-sector derivation beyond background transport if C0 requires non-background TEP variables
- [ ] Run Planck 2018 likelihood or documented compressed CMB covariance

## Usage

```bash
cd external/class_tep_mod
./apply_tep_patch.sh /path/to/class_public
make clean && make
```

Then use the modified classy:
```python
from classy import Class
cosmo = Class()
cosmo.set({'tep_mode': 'yes', 'tep_Sigma_0': 0.001, 'tep_epsilon_T': 1.0})
```

For the C0 resolver, run with:

```bash
TEP_CLASS_PYTHONPATH=external/class_tep_mod/python_runtime python3 scripts/steps/step_028_cmb_full_spectra.py
```
