# Temporal Equivalence Principle: A Covariant Alternative to Cosmic Expansion

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20370143.svg)](https://doi.org/10.5281/zenodo.20370143)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

**Author:** Matthew Lukin Smawfield  
**Version:** v0.5 (Athens)  
**First published:** 5 July 2026 · **Last updated:** 24 October 2026  
**Status:** Preprint  
**DOI:** [10.5281/zenodo.20370143](https://doi.org/10.5281/zenodo.20370143)  
**Website:** [https://mlsmawfield.com/tep/c0/](https://mlsmawfield.com/tep/c0/)  
**Paper Series:** TEP Series · Group F (VI. Cosmological Transport, Horizon Boundary & Standard Sirens) · Paper 26 (`F26`)

## Abstract

<!-- Generated from site/components/1_abstract.html; update the manuscript source, then run tools/sync_current_corpus_readmes.py. -->

This paper presents a direct empirical challenge to the necessity of primitive cosmic expansion. In the Temporal Equivalence Principle framework, observed redshift is reconstructed as conformal proper-time transport, $1+z=A_0/A_{\text{em}}$ rather than as stretching of a spatial scale factor. Standard cosmology interprets observational redshift and luminosity distance scaling as evidence of a stretching spatial metric, parameterized by the Friedmann-Lemaître-Robertson-Walker (FLRW) scale factor $a(t)$. The observational role played by the FLRW scale factor is mapped, within the TEP conformal-frame construction, onto the temporal clock-rate field $A(\phi)$. In the tested late-time background sector, the perceived acceleration normally attributed to dark energy, $\Lambda$, is reconstructed as the kinematic clock-drift rate contribution of the Temporal Shear field, $\Omega_\phi$, in the static-manifold distance-redshift Jacobian (rather than a $w = +1$ Friedmann fluid driving metric expansion $a(t)$).

The core relation is $1+z = A_0/A_{\text{em}}$. In the static conformal interpretation developed here, intergalactic separations are not treated as primitively expanding; the apparent expansion is reconstructed through temporal transport. In this framework, the limit conventionally written as $a\to0$ is re-expressed as $A_{\text{clock}}\to0$: a TEP temporal-horizon boundary of observational clock transport, not a zero-volume spatial singularity. The temporal-horizon background and linear-mode closure are supplied by companion papers TEP-TH and TEP-HC; C0 imports and cross-checks the TEP-HC linear-growth output while focusing on the empirical supernova-sector test.

Using 1,701 Pantheon+ Type Ia supernova measurements (1,543 distinct events) with the full covariance matrix and the official peculiar-velocity-corrected Hubble-diagram redshifts, a pure conformal reconstruction exactly reproduces the $\Lambda$CDM homogeneous distance-modulus relation, demonstrating that the background Hubble diagram does not uniquely select an expanding spatial metric. More strongly, the conservative physical no-$\Lambda$ temporal-shear branch with fixed line-of-sight turnover $z_{\rm los}=5$ (using the acoustic-sector scale as a conservative reference) improves the standardized supernova likelihood by $\Delta\chi² \simeq -3.81$ relative to baseline $\Lambda$CDM and achieves a Bayes factor of approximately 5.8, nominally in the "substantial" band of the Jeffreys scale — reported alongside its prior-sensitivity grid, and its reduction to BF $=0.52$–$0.96$ after low-redshift removal ($z>0.05$–$0.01$): the improvement localizes to the 77 SH0ES calibrators and the $z>1$ tail, so the SNe-only preference is evidence for a transport-shaped deviation whose physical attribution is tested by the joint and environmental channels. Under the canonical suppression exponent $n_T=2$ the free-$z_{\rm los}$ branch marginalizes the turnover scale, with its likelihood maximum at the admissible turnover $z_{\rm los}\simeq5.3$ and its posterior volume extending over the nearly flat unsuppressed tail, and returns BF $\approx$ 5.4, while the fixed $z_{\rm los}=100$ benchmark gives a comparable BF $\approx$ 5.4 ($\Delta\chi² \simeq -3.69$); the admissible branches thus saturate at a moderate, internally consistent level rather than at extreme amplification, and both standard dark-energy extensions (wCDM, CPL) return weaker evidence at matched or greater freedom. The conservative physical model with $z_{\rm los}=5$—which turns over before the acoustic epoch and is therefore consistent with the joint CMB bound on the homogeneous shear component—demonstrates that a matter-only temporal-shear geometry is competitive with $\Lambda$CDM in the late-time SNe distance-redshift sector. A single-history joint fit to Pantheon+, DESI DR2 BAO and the Planck acoustic angle (Step 03-14) selects the native Jordan-frame TEP background of Paper 18, which is preferred over $\Lambda$CDM by $\Delta\chi²=-8.1$ ($\Delta{\rm AIC}=-4.1$), with the preference growing from DESI DR1 ($-6.3$) to DR2; the M1 kernel serves as the supernova-sector distance law. The same framework supplies a galaxy-mass-locked host-environment diagnostic from the nested clock-rate structure of the scalar sector: the transferred Cepheid response predicts the established massive-host-brighter orientation at the literature-step order, while overshooting the residual in the already host-corrected Pantheon+ product; freeing the response amplitude on the continuous per-supernova depth coordinate returns a $2.9\sigma$ detection in the predicted direction (deeper wells brighter) at roughly a tenth of the naive transfer normalization. The environmental prediction is tested directly: separating the corpus's two nested channels — ambient host/density depth, which enhances the transport amplitude, and within-host radial depth, which screens it — yields a two-channel modulation at permutation significance $p=0.036$ (1000 realisations) that survives redshift-projection, anchor-removal, band-pass, survey-offset, covariance-perturbation, leave-one-survey-out, and out-of-sample transfer controls, and re-derives with the correct signs inside each resolved host-association scheme independently — including positive transfer to held-out SDSS and DES; PS1MD transfer is negative; after jointly modeling the light-curve-quality amplitude structure that competes with it, the environmental attribution reaches $\Delta\chi²\simeq-9.13$ ($p=0.015$ against 1000 whole-event permutations within survey-composition and redshift strata), concentrated in the transport-sensitive bulk. The predicted depth-response coefficient is measured directly on the per-supernova depth coordinate: $\hat\kappa_{\rm SN}=+3.2\times10^{4}\pm1.1\times10^{4}$ mag on Pantheon+ ($2.9\sigma$, predicted sign, calibrator-block localized, permutation $p<0.005$), replicated on the independent DES-Dovekie catalog at $+8.9\times10^{4}\pm2.3\times10^{4}$ ($3.9\sigma$, concentrated at $4.8\sigma$ in the high-mass half) — the joint direction is positive at $>4\sigma$ across two compilations with independent standardization, at a normalization consistent with channel-specific response projection.

Companion papers establish the theoretical foundations: TEP-HC (Paper 18) provides the Boltzmann-level acoustic-scale preservation proof under the native hi_class tep_mode implementation, and TEP-TH develops the nonsingular temporal-horizon closure. The current paper focuses on the empirical supernova-sector test and the deterministic falsification pipeline.

Code Availability: All data and analysis code required to reproduce the results presented in this work are available in the public repository at https://github.com/matthewsmawfield/TEP-C0.

Keywords: temporal equivalence principle, static conformal geometry, cosmology, dark energy, supernovae, Bayesian inference, modified gravity, temporal shear

Current manuscript source: [1_abstract.html](site/components/1_abstract.html).

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
    ω_b=0.0224,    # Physical baryon density
    ω_cdm=0.12,    # Physical CDM density
    ε_T=6.7e-6,  # Temporal-shear amplitude
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
| Ω_γ | 5.04 × 10⁻⁵ | ~5 × 10⁻⁵ | ✅ |
| z_rec | 1098.6 | ~1100 | ✅ |
| Y_p | 0.296 | ~0.25 | ✅ |
| D/H | 2.60 × 10⁻⁵ | ~2.6 × 10⁻⁵ | ✅ |

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
A(z) = exp(ε_T * ln(1+z) * S(z)),  S(z) = exp(-(z/z_T)^n_T)
M(z) = A(z) / (1 - α_A(z))
H_TEP(z) = H_LCDM(z) * M(z)
```

Where:
- `H_TEP`: TEP-modified Hubble parameter
- `ε_T`: homogeneous temporal-shear amplitude (screened at z >> z_T)
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

CC-BY-4.0 - see [LICENSE](LICENSE) file.

## Status

- ✅ Background cosmology: Research grade
- ✅ Recombination: Research grade
- ✅ BBN: Research grade (simplified network)
- ⚠️ CMB C_l: Working (simplified transfer)

## Contact

For questions or issues, please open a GitHub issue.
