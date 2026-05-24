# TEP-C0: Temporal Equivalence Principle Cosmological Pipeline

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

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

The TEP framework modifies the standard FLRW expansion:

```
H_TEP(z) = H_LCDM(z) × Γ_TEP(z)
Γ_TEP(z) = exp(Sigma_0 × c/H_0 × ln(1+z))
```

Where:
- `H_TEP`: TEP-modified Hubble parameter
- `Sigma_0`: TEP shear amplitude
- `Γ_TEP`: Path enhancement factor

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
