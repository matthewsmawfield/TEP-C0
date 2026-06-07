# Temporal Equivalence Principle: A Covariant Alternative to Cosmic Expansion

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Status:** In Development

## Abstract

Abstract This paper develops the cosmological extension of the Temporal Equivalence Principle (TEP): the hypothesis that observational evidence normally interpreted as cosmic expansion may involve large-scale Temporal Shear. In TEP, matter clocks and photon phases evolve in the causal matter metric g_{}, with the conformal clock-rate field A() defining the Temporal Shear _ = _  A(). Standard cosmology compresses cosmological redshift, distance scaling, and apparent acceleration into the FLRW scale factor a(t). TEP establishes that a(t) is an effective variable reconstructed from accumulated Temporal Shear and Temporal Topology along cosmological lines of sight. The core relation is a_{eff}() = [-_ _^{eff} d], where 1+z_T = a_{eff}^{-1}. In the homogeneous integrable limit, this reproduces the FLRW relation 1+z = a_0/a_{em}. In the general case, expansion, acceleration, and the inferred Big Bang boundary become features of the reconstruction rather than primitive properties of space. Utilizing a dual-domain Bayesian synthesis of 1,701 Pantheon+ supernovae and Planck 2018 acoustic anchors, the analysis reveals a critical structural separation. In the late universe, nested sampling over the supernovae strictly prefers the TEP geometry over standard  and phenomenological dark energy (BIC = -1279.21, BF = 131.6). The converged 120,960-accepted-step joint Cobaya MCMC demonstrates that the pristine global CMB bounds the macroscopic temporal shear to zero, acting as the ultimate cosmological boundary condition. By formalizing environmental state suppression across the hierarchically structured cosmic web, the framework perfectly reconciles this divergence. The theory natively isolates massive anomalies—providing a rigorous geometric origin for the supernova "mass step" and resolving the Hubble tension (Paper 11) and JWST high-redshift mass anomalies (Paper 12)—entirely within local and intermediate scales, while consistently protecting the standard cosmological background. The framework culminates in a preregistered empirical testing program targeting the theory's central hallmark: synchronization holonomy (H). Driven by non-zero disformal proper-time transport, H provides a directly observable, convention-independent metric of non-integrability, guiding a new class of multi-leg time-transfer experiments. Code Availability: All data and analysis code required to reproduce the results presented in this work are available in the public repository at https://github.com/matthewsmawfield/TEP-C0 .

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

# Install package
pip install -e .

# Download public datasets
python -m tep_c0.data.download
```

## Quick Start

```python
from tep_c0.core.background import TEPBackground
from tep_c0.core.recombination import RecombinationHistory

# Create TEP-modified cosmology
bg = TEPBackground(
    H0=70.0,           # Hubble constant [km/s/Mpc]
    Omega_b=0.045,     # Baryon density
    Omega_cdm=0.25,    # CDM density
    Omega_Lambda=0.7,  # Dark energy
    Sigma_0=0.001      # TEP shear parameter
)

# Compute distances
z = 1.0
d_L = bg.luminosity_distance(z)
d_A = bg.angular_diameter_distance(z)

# Recombination epoch
rec = RecombinationHistory(bg)
z_rec = rec.z_rec()  # ~1100
```

## Pipeline Execution

### Full Pipeline
```bash
python -m tep_c0.pipeline.run_all
```

### Individual Steps
```bash
# Step 0: Download data
python -m tep_c0.data.download

# Step 22: Cosmological fitting
python -m tep_c0.analysis.step_022

# Step 17: CMB computation
python -m tep_c0.analysis.step_017

# Step 29: BBN computation
python -m tep_c0.analysis.step_029
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
├── src/tep_c0/           # Core package
│   ├── core/             # Physics modules
│   │   ├── background.py    # Friedmann equations
│   │   ├── recombination.py # Saha/Peebles
│   │   └── bbn.py           # Nucleosynthesis
│   ├── analysis/         # Pipeline steps
│   ├── data/             # Data management
│   └── utils/            # Common utilities
├── data/                 # Data directory
│   ├── raw/              # Downloaded datasets
│   └── processed/        # Intermediate outputs
├── results/              # Pipeline outputs
│   ├── figures/          # Generated plots
│   └── outputs/          # JSON results
├── tests/                # Test suite
├── docs/                 # Documentation
└── scripts/              # Legacy scripts
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