# Temporal Equivalence Principle: A Covariant Alternative to Cosmic Expansion

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Status:** In Development

## Abstract

This paper presents a direct empirical challenge to the necessity of primitive cosmic expansion. In the Temporal Equivalence Principle framework, observed redshift is reconstructed as conformal proper-time transport, $1+z=A_0/A_{\rm em}$, rather than as stretching of a spatial scale factor. Standard cosmology interprets observational redshift and luminosity distance scaling as evidence of a stretching spatial metric, parameterized by the FLRW scale factor $a(t)$. The observational role played by the FLRW scale factor is mapped, within the TEP conformal-frame construction, onto the temporal clock-rate field $A(\phi)$. In the tested late-time background sector, the perceived acceleration normally attributed to Dark Energy, $\Lambda$, is reconstructed as the kinetic energy density of the Temporal Shear field, $\Omega_\phi$.

The core relation is $1+z = A_0/A_{\text{em}}$. In the static conformal interpretation developed here, intergalactic separations are not treated as primitively expanding; the apparent expansion is reconstructed through temporal transport. In this framework, the limit conventionally written as $a\to0$ is re-expressed as $A_{\text{clock}}\to0$: a TEP temporal-horizon boundary of observational clock transport, not a zero-volume spatial singularity. The framework is formally closed from the temporal-horizon background ($a_{\rm eff} \to 0$) down to linear mode evolution ($k \sim 0.5 \, h/\text{Mpc}$) through the exact Bellini-Sawicki EFT mapping implemented in TEP-HC (Paper 18), which natively yields $\sigma_8 \approx 0.825$ in agreement with Planck.

Using 1,701 Pantheon+ Type Ia supernovae with the full covariance matrix, a pure conformal reconstruction exactly reproduces the $\Lambda$CDM homogeneous distance-modulus relation, demonstrating that the background Hubble diagram does not uniquely select an expanding spatial metric. More strongly, the conservative physical no-$\Lambda$ temporal-shear branch with fixed $z_T=5$ improves the standardized supernova likelihood by $\Delta\chi^2 \simeq -3.4$ relative to baseline $\Lambda$CDM and achieves a Bayes factor of BF$\simeq 4.6$, classified as "substantial" evidence on the Jeffreys scale. The fixed $z_T=100$ benchmark gives the strongest evidence, BF$\simeq 61.8$ ($\Delta\chi^2 \simeq -7.5$), while the broad free-$z_T$ model gives BF$\simeq 40.3$. The conservative physical model with $z_T=5$ already demonstrates that a matter-only temporal-shear geometry is competitive with $\Lambda$CDM in the late-time SNe distance-redshift sector. These results establish positive supernova-sector evidence that apparent acceleration can be reconstructed as temporal transport rather than primitive dark energy.

Companion papers establish the theoretical foundations: TEP-HC (Paper 18) provides the Boltzmann-level acoustic-scale preservation proof under the native hi_class `tep_mode` implementation, and TEP-TH develops the nonsingular temporal-horizon closure. The current paper focuses on the empirical supernova-sector test and the deterministic falsification pipeline.

## Overview

TEP-C0 is a research-grade cosmological pipeline implementing the Temporal Equivalence Principle (TEP) framework. It provides:

- **Supernova cosmology** with full covariance analysis (Pantheon+ dataset)
- **CMB physics** with TEP-modified Boltzmann equations
- **BBN nucleosynthesis** with path-enhanced expansion rates
- **Model comparison** via Bayesian evidence

## Installation

```bash
# Clone repository
git clone <repository-url>
cd TEP-C0

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Build patched CLASS (optional; required for CMB steps)
cd external/class
make
```

## Quick Start

```python
import sys
sys.path.insert(0, '.')
from core.cosmology import TEP_COSMOLOGY

# Create TEP-modified cosmology
cosmo = TEP_COSMOLOGY(
    H0=70.0,           # Hubble constant [km/s/Mpc]
    omega_b=0.0224,    # Physical baryon density
    omega_cdm=0.12,    # Physical CDM density
    epsilon_T=6.7e-6,  # Temporal-shear amplitude
    z_T=5.0,           # Transition redshift
    n_T=2.0            # Transition steepness
)

# Compute distances
z = 1.0
d_L = cosmo.luminosity_distance(z)
d_A = cosmo.angular_diameter_distance(z)
```

## Pipeline Execution

### Full Pipeline
```bash
python scripts/run_pipeline.py
```

### Individual Steps
```bash
# Step 01: Download data
python scripts/steps/step_01_01_data_download.py

# Step 22: Three-model comparison
python scripts/steps/step_03_01_three_model_comparison.py

# Step 17: CMB Boltzmann
python scripts/steps/step_05_03_cmb_boltzmann.py

# Step 29: BBN preservation
python scripts/steps/step_05_07_bbn_preservation.py
```

## Data Sources

All data is downloaded from public repositories:

- **Pantheon+**: [GitHub Repository](https://github.com/PantheonPlusSH0ES/DataRelease)
  - 1701 Type Ia supernovae
  - Full statistical + systematic covariance
  - SHA-256: `1cb0fc379ef066af...`

- **FIRAS CMB**: [NASA LAMBDA](https://lambda.gsfc.nasa.gov)
  - COBE/FIRAS monopole spectrum
  - Perfect blackbody validation

- **BAO**: [Zenodo Compilation](https://zenodo.org/records/16285883)
  - Multi-survey BAO constraints

## Physics Validation

| Observable | Computed | Target | Status |
|------------|----------|--------|--------|
| Ω_γ | 5.04×10⁻⁵ | ~5×10⁻⁵ | ✅ |
| z_rec | 1098.6 | ~1100 | ✅ |
| Y_p | 0.296 | ~0.25 | ✅ |
| D/H | 2.60×10⁻⁵ | ~2.6×10⁻⁵ | ✅ |

## Project Structure

```
TEP-C0/
├── core/                 # Physics modules (shared with TEP-HC/TEP-TH)
│   ├── cosmology.py         # TEP background & EFT sector
│   ├── conformal_scaling.py # Jordan-frame conformal factors
│   ├── scalar_field.py      # Temporal-shear field
│   ├── screening.py         # Screening dictionary & PPN
│   ├── evidence.py          # Evidence-gate infrastructure
│   └── constants.py         # Physical constants & aliases
├── scripts/              # Analysis pipeline
│   ├── steps/               # Individual pipeline steps (01-08)
│   ├── utils/               # Plotting, PDF generation, logging
│   └── run_pipeline.py      # Full pipeline runner
├── data/                 # Data directory
│   ├── raw/                 # Downloaded datasets (Pantheon+, FIRAS, BAO)
│   └── processed/           # Ingested catalogues & registries
├── results/              # Pipeline outputs
│   ├── figures/             # Generated plots
│   ├── outputs/             # MCMC chain files & CLASS samples
│   └── *.json               # Step results & evidence matrices
├── site/                 # Static-site manuscript builder
│   ├── components/          # HTML sections
│   └── public/docs/         # Published PDF
├── external/             # Patched dependencies
│   ├── class/               # TEP-CLASS Boltzmann solver
│   └── alterbbn/            # BBN nucleosynthesis
└── manuscripts/          # Generated Markdown manuscripts
```

## Testing

```bash
# Run test suite
pytest tests/

# Run specific test
pytest tests/test_physics.py -v

# Run audit
python -m tep_c0.utils.audit
```

## Methodology

### TEP Modification

The native TEP background (shared with TEP-HC Paper 18 via `core/cosmology.py` and the hi_class `tep_mode` patch) modifies the FLRW expansion through the Jordan-frame factor:

```
f_T(z) = ln(1+z) * exp(-(z/z_T)^n_T)
A(z) = exp(epsilon_T * ln(1+z) * S(z)),  S(z) = exp(-(z/z_T)^n_T)
M(z) = A(z) / (1 - alpha_A(z))
H_TEP(z) = H_LCDM(z) * M(z)
```

Where:
- `H_TEP`: TEP-modified Hubble parameter
- `epsilon_T`: homogeneous temporal-shear amplitude (screened at z >> z_T)
- `z_T`, `n_T`: transition redshift and steepness (default 5.0, 2.0)
- `M(z)`: exact Jordan-frame conformal factor (not the legacy Γ_TEP exponential approximation)

### Cosmological Fitting

- **Likelihood**: Gaussian with full covariance
- **Inference**: Nested sampling (dynesty) + MCMC (emcee)
- **Convergence**: R-hat < 1.05, dlogZ < 0.1
- **Evidence**: Bayesian model comparison

### CMB Computation

- **Background**: TEP-modified Friedmann equations
- **Recombination**: Saha + Peebles with TEP H(T)
- **Acoustic scale**: Proper sound horizon integration

### BBN Computation

- **Network**: Simplified nuclear network
- **Physics**: Neutron freeze-out + decay
- **Abundances**: Y_p, D/H, He-3/H, Li-7/H

## Citation

If using this pipeline, please cite:

```bibtex
@software{tep_c0,
  title = {TEP-C0: Temporal Equivalence Principle Cosmological Pipeline},
  year = {2026},
  url = {<repository-url>}
}

@dataset{pantheon_plus,
  title = {Pantheon+ SH0ES Data Release},
  author = {Scolnic, D. et al.},
  year = {2022},
  journal = {ApJ},
  volume = {938},
  pages = {113}
}
```

## License

MIT License - see [LICENSE](LICENSE) file.

## Status

- ✅ Background cosmology: Research grade
- ✅ Recombination: Research grade
- ✅ BBN: Research grade (simplified network)
- ⚠️ CMB C_l: Working (simplified transfer)

## Contact

For questions or issues, please open a GitHub issue.