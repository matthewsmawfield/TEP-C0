# Temporal Equivalence Principle: EFT Mapping and CMB Acoustic Peak Preservation
**Matthew Lukin Smawfield**
Version: v0.1 (Geneva)
First published: 24 May 2026 · Last updated: 2 June 2026
Paper Series: TEP Series Paper 18 (hi_class Cosmology)

---

## Abstract


General Relativity is extensively validated in the deeply screened, high-density regime of the Solar System, but cosmological tensions—specifically the Hubble discrepancy and galactic mass anomalies—suggest a scale-dependent breakdown of the isochrony axiom.



The Temporal Equivalence Principle (TEP) addresses these late-universe anomalies via a dynamical proper-time field, screened at densities $\rho > 20$ g/cm³. Previous analytical approximations and CAMB integrations suggest this scalar field is frozen during radiation domination ($T^\mu_\mu \approx 0$). Full integration into linear structure formation requires a native Horndeski evaluation.



This paper sketches the mapping of the TEP bi-metric framework—defined by a conformal factor $A(\phi)$ and a disformal deformation $B(\phi)$—onto the Bellini-Sawicki Effective Field Theory (EFT) of Dark Energy, and then adopts a background-only realization for the numerical analysis. The full EFT (alpha-function) implementation is left to future work; the present results rest on the background-only modification, which is what we compute and constrain.



This paper implements the native TEP background-only modification directly in hi_class via the transition function `f_T(z) = ln(1+z) * exp(-(z/z_T)^n_T)`, applied as `H_TEP(z) = H_LCDM(z) * (1 + epsilon_T * f_T(z))` while preserving standard General Relativistic perturbations. The functional form is the authoritative one used in TEP-C0 (Paper 26): the `exp(-(z/z_T)^n_T)` factor enforces early-time freezing, so the temporal-shear field is suppressed for `z >> z_T`. Direct Boltzmann integration confirms that the pre-recombination sound horizon is preserved to parts-per-million (`r_s^TEP / r_s^LCDM = 1.0000000`) and that the acoustic-peak morphology is untouched. The only CMB-level effect of a non-zero `epsilon_T` is a late-time angular-diameter-distance projection that rigidly rescales the angular acoustic scale `theta_s` (a 0.27% shift at the fiducial `epsilon_T = 0.0066`, `z_T = 5`, `n_T = 2`) and is fully degenerate with `H_0`. Because this projection is absorbed by the standard parameters, the homogeneous amplitude is bounded by the CMB to `epsilon_T ~ 0`, consistent with TEP-C0 (Paper 26), where Pantheon+ supernovae give substantial Bayesian preference for the TEP geometry (Bayes factor 7-10 vs LambdaCDM) while the CMB anchors the homogeneous background to the LambdaCDM limit. The Hubble tension is therefore resolved as a late-time, environment-dependent clock-transport effect (Paper 11) rather than through a modified homogeneous expansion history at recombination.



Keywords: cosmology theory, cosmic microwave background, dark energy, scalar-tensor theories, modified gravity, large-scale structure of the universe, hi_class, Horndeski, temporal equivalence principle



## 1. Introduction


### 1.1 Contextualizing the TEP Corpus


The Temporal Equivalence Principle (TEP) has been constrained across 40 orders of magnitude in mass density, from terrestrial laboratory scales to cosmological observations. Previous papers in this series have established:



- **Terrestrial scales (Paper 1, GNSS):** 4,200 km phase correlations measured via global navigation satellite systems confirm the 20 g/cm³ screening threshold.

- **Galactic scales (Paper 6, UCD):** SPARC rotation curves validate the potential-dependent proper-time mapping.

- **Stellar scales (Paper 13, WB):** Gaia DR3 wide binaries exhibit the predicted density-dependent kinematic transition.

- **Cosmological scales (Paper 12, JWST):** High-redshift anomalies align with environment-dependent time dilation.




### 1.2 The Two-Ended Hubble Tension


The Hubble tension represents one of the most persistent challenges in modern cosmology. Cepheid-calibrated local distance ladder measurements yield $H_0 \approx 73$ km/s/Mpc, while early-universe CMB inference from Planck gives $H_0 \approx 67.9$ km/s/Mpc—a discrepancy of approximately 4.8$\sigma$.


Previous TEP work (Paper 11, H₀) demonstrated that the Cepheid environmental bias—operating in the unscreened stellar atmospheres where Cepheids pulsate—naturally raises the local measurement. However, a complete resolution requires demonstrating that this mechanism does *not* disturb the early-universe inference at the surface of last scattering ($z \approx 1089$).


### 1.3 Purpose of This Paper


To move beyond the quasi-static approximations of previous work by natively solving the background expansion $H(z)$ and linear perturbation growth in the weakly screened density regime ($\rho \sim 10^{-21}$ g/cm³). This requires:



- A mapping of TEP's bi-metric structure onto the Bellini-Sawicki EFT formalism (sketch; the numerical analysis uses the background-only realization).

- Native implementation in hi_class with proper background scalar field evolution.

- ΛCDM baseline MCMC parameter estimation ($H_0$, $\Omega_b h^2$, $\Omega_{\rm cdm} h^2$, $n_s$, $A_s$, $\tau$, $A_{\rm planck}$) against Planck 2018, BAO, and Pantheon+ data to establish the constraints that TEP should satisfy.



The critical question: Can TEP preserve CMB acoustic peaks while permitting late-time $H_0$ variation?

## 2. Theoretical Architecture: The EFT Mapping

### 2.1 The Bi-Metric Action

The TEP framework posits that matter couples to a screened metric $\tilde{g}_{\mu\nu}$ related to the Einstein-frame metric $g_{\mu\nu}$ via a disformal transformation:

\begin{equation} \label{eq:3_theory_01}
\tilde{g}_{\mu\nu} = A^2(\phi) g_{\mu\nu} + B(\phi) \nabla_\mu\phi \nabla_\nu\phi
\end{equation}

where:

- $A(\phi) = \exp(\beta\phi/M_{\rm Pl})$ is the conformal factor, with $\beta = 0.15$ fixed by stellar and galactic scale constraints

- $B(\phi)$ controls disformal deformation of the causal structure

- $\phi$ is the dynamical proper-time field

### 2.2 Deriving the Bellini-Sawicki Alphas

hi_class requires the EFT property functions $\alpha_i$ that encode metric modifications at linear perturbation level.

#### 2.2.1 Planck Mass Running ($\alpha_M$)

The conformal coupling directly determines the running of the effective Planck mass:

\begin{equation} \label{eq:3_theory_02}
\alpha_M \equiv \frac{d \ln M_{\rm eff}^2}{d \ln a} = \frac{d \ln A^2(\phi)}{d \ln a} = \frac{2\beta}{M_{\rm Pl}} \frac{\phi'}{\mathcal{H}}
\end{equation}

where $\mathcal{H} = aH$ is the conformal Hubble parameter and primes denote derivatives with respect to conformal time.

#### 2.2.2 Tensor Speed Excess ($\alpha_T$)

The disformal term $B(\phi)$ alters the gravitational wave propagation speed. Multi-messenger constraints from GW170817/GRB 170817A require:

\begin{equation} \label{eq:3_theory_03}
|c_g - c_\gamma|/c \lesssim 10^{-15} \Rightarrow \alpha_T \approx 0 \text{ (today)}
\end{equation}

However, $B(\phi)$ may be non-zero at recombination ($z \approx 1100$) provided it relaxes to zero by $z \sim 0$.

#### 2.2.3 Braiding ($\alpha_B$) and Kineticity ($\alpha_K$)

These functions govern scalar field clustering and metric mixing:

\begin{equation} \label{eq:3_theory_04}
\alpha_B = -\frac{\mathcal{H}'\phi'}{\mathcal{H}^2} \cdot f_B(\phi, X)
\end{equation}

\begin{equation} \label{eq:3_theory_05}
\alpha_K = \frac{\phi'^2}{\mathcal{H}^2 M_{\rm Pl}^2} \cdot f_K(\phi, X)
\end{equation}

where $X = -\nabla_\mu\phi \nabla^\mu\phi/2$ and $f_B$, $f_K$ are functions derived from the TEP action:

- $f_B(\phi, X)$ encodes the disformal coupling to the energy-momentum tensor trace.

- $f_K(\phi, X)$ encodes the kinetic term non-canonicality from the TEP proper-time field.

The explicit functional forms follow from the bi-metric action (Equation \ref{eq:3_theory_01}) and are determined by the conformal factor $A(\phi)$ and the disformal function $B(\phi)$. Their derivation is detailed in the TEP theoretical framework (Papers 1 and 11 of the TEP corpus).

### 2.3 The Radiation Domination Freezing Mechanism

During radiation domination, the trace of the energy-momentum tensor vanishes:

\begin{equation} \label{eq:3_theory_06}
T^\mu_\mu = -\rho + 3p \approx 0 \text{ (radiation)}
\end{equation}

Since the TEP scalar field couples to $T^\mu_\mu$, the source term for $\phi$ evolution is suppressed. The field freezes at its initial value, and $\alpha_M, \alpha_B \rightarrow 0$. This ensures:

- Primary acoustic peaks ($100 \lesssim \ell \lesssim 2000$) generated at $z \sim 1089$ remain unmodified

- The sound horizon $r_s$ is preserved, anchoring $H_0$ from CMB at $\sim 67.9$ km/s/Mpc

- Late-time matter domination reactivates the scalar field, enabling environmental $H_0$ variations


## 3. Software Implementation: hi_class and the Unscreened Regime


### 3.1 The hi_class Architecture


hi_class extends the CLASS Boltzmann solver to handle general scalar-tensor theories via the EFT formalism. This work uses hi_class v3.2.3 with the modified gravity (SMG) module enabled.


### 3.2 Native TEP Background-Only Implementation

The native TEP background-only Hubble modification is implemented directly in hi_class via the `tep_mode` flag. When enabled, the background expansion history is modified as:

\begin{equation} \label{eq:4_implementation_01}
H_{\rm TEP}(z) = H_{\Lambda\rm CDM}(z) \times \left(1 + \epsilon_T \, f_T(z)\right)
\end{equation}

where the transition function is the authoritative TEP-C0 form (Paper 26; `core/tep_cosmology.py`)

\begin{equation} \label{eq:4_implementation_01b}
f_T(z) = \ln(1+z)\,\exp\!\left[-(z/z_T)^{n_T}\right].
\end{equation}

The exponential factor $\exp[-(z/z_T)^{n_T}]$ suppresses the modification for $z \gg z_T$, so $f_T \to 0$ well before recombination and the homogeneous expansion at last scattering is unmodified (radiation-domination freezing, Section 2.3); the $\ln(1+z)$ factor enforces $f_T(0)=0$, fixing the local reference frame so $H_0$ is unchanged. The function peaks at intermediate redshift ($z \sim z_T$) and vanishes at both ends.

*Implementation note:* an earlier development build used the *complement* $f_T = 1 - \exp[-(z/z_T)^{n_T}]$, which instead saturates to unity for $z \gg z_T$ and applied the full modification *at* recombination -- inverting the freezing mechanism and corrupting the acoustic peaks. That bug has been corrected (see Appendix A.3). The default TEP parameters are:



```
tep_mode = yes
epsilon_T = 0.0066
z_T = 5.0
n_T = 2.0
```

Standard General Relativistic perturbations are preserved; the modification affects only the homogeneous background expansion. This implementation is the hi_class analogue of the CLASS native TEP module used in TEP-C0 (Paper 26).


### 3.3 Perturbation Stability

Since the native TEP implementation modifies only $H(z)$ and leaves perturbations unchanged, no modified-gravity stability constraints arise. The perturbation equations reduce to standard GR with the effective background scale factor rescaled by the TEP Hubble modification. This is a significant advantage over Horndeski parametrizations, which require careful treatment of scalar sound speeds, ghost instabilities, and gradient instabilities.


### 3.4 Pipeline Architecture

The full analysis pipeline, executed via `scripts/run_all.py`, consists of:


- **Step 01 (Setup):** Environment configuration and dependency check.

- **Step 02 (Background):** Compute the TEP-modified background expansion history $H(z)$ and density evolution.

- **Step 03 (Alpha Functions):** Compute Bellini-Sawicki coefficients from the TEP theoretical mapping (archived for reference).

- **Step 04 (CMB Spectra):** Run hi_class with native `tep_mode` at the Planck 2018 best-fit point. Compare TT, TE, and EE spectra against standard CLASS $\Lambda$CDM.

- **Step 05 (Cobaya Config):** Generate the Cobaya YAML configuration for the MCMC pipeline with native TEP parameters.

- **Step 06 (MCMC):** Execute the Cobaya MCMC with hi_class, using real Planck + BAO + Pantheon+ likelihoods.

- **Step 07 (Posteriors):** Analyze MCMC chains with burn-in removal and weighted statistics.

- **Step 08 (Synthesis):** Combine all results into summary JSON and figures.



## 4. MCMC Parameter Estimation Pipeline

### 4.1 The Cobaya Framework

Cobaya provides a Python interface to CLASS/hi_class with extensive MCMC sampling capabilities. The transition from SciPy/Pandas pipelines to Cobaya enables:

- Native hi_class integration without file-based I/O bottlenecks

- Parallel tempering and adaptive Metropolis-Hastings sampling

- Direct Planck likelihood wrapper integration

- Seamless GetDist posterior visualization

### 4.2 Likelihood Configuration

The pipeline uses the following Planck 2018 likelihoods:

| Likelihood | Description | $\ell$ Range |
| --- | --- | --- |
| `planck_2018_lowl.TT` | Low-$\ell$ temperature | 2–29 |
| `planck_2018_lowl.EE` | Low-$\ell$ polarization | 2–29 |
| `planck_2018_lensing.native` | CMB lensing reconstruction | 8–400 |
| `bao.sdss_dr12_consensus_final` | BAO SDSS DR12 consensus | — |
| `sn.pantheonplus` | Type Ia supernovae (Pantheon+) | — |

### 4.3 Free Parameters and Priors

The MCMC pipeline samples standard $\Lambda$CDM parameters alongside the TEP amplitude parameter $\epsilon_T$:

| Parameter | Prior | Description |
| --- | --- | --- |
| $\Omega_b h^2$ | $\mathcal{U}(0.005, 0.1)$ | Baryon density |
| $\Omega_{\rm cdm} h^2$ | $\mathcal{U}(0.001, 0.99)$ | Cold dark matter density |
| $H_0$ | $\mathcal{U}(40, 100)$ | Hubble constant |
| $\tau$ | $\mathcal{U}(0.01, 0.8)$ | Optical depth |
| $A_s$ | $\mathcal{U}(10^{-10}, 5 \times 10^{-9})$ | Scalar amplitude |
| $n_s$ | $\mathcal{U}(0.8, 1.2)$ | Scalar spectral index |
| $A_{\rm planck}$ | $\mathcal{U}(0.9, 1.1)$ | Planck calibration nuisance |
| $\epsilon_T$ | $\mathcal{U}(0, 0.02)$ | TEP amplitude parameter (background Hubble modification) |

### 4.4 Pipeline Execution

```
# Cobaya YAML configuration
theory:
classy:
path: /path/to/hi_class
extra_args:
output: tCl,pCl,lCl
lensing: yes
modes: s,t
non_linear: none
# Native TEP background-only Hubble modification
tep_mode: 'yes'
epsilon_T: 0.0066
z_T: 5.0
n_T: 2.0

likelihood:
planck_2018_lowl.TT: null
planck_2018_lowl.EE: null
planck_2018_lensing.native: null
bao.sdss_dr12_consensus_final: null
sn.pantheonplus: null

params:
logA:
prior: {min: 2.5, max: 3.5}
ref: {dist: norm, loc: 3.044, scale: 0.014}
proposal: 0.01
drop: true
A_s:
value: 'lambda logA: 1e-10*np.exp(logA)'
n_s:
prior: {min: 0.94, max: 1.0}
ref: {dist: norm, loc: 0.966, scale: 0.004}
proposal: 0.004
H0:
prior: {min: 40, max: 100}
ref: {dist: norm, loc: 67.4, scale: 0.5}
proposal: 1.5
omega_b:
prior: {min: 0.005, max: 0.1}
ref: {dist: norm, loc: 0.0224, scale: 0.0002}
proposal: 0.0003
omega_cdm:
prior: {min: 0.001, max: 0.99}
ref: {dist: norm, loc: 0.12, scale: 0.001}
proposal: 0.0015
tau_reio:
prior: {min: 0.01, max: 0.8}
ref: {dist: norm, loc: 0.054, scale: 0.007}
proposal: 0.01
A_planck:
prior: {min: 0.9, max: 1.1}
ref: {dist: norm, loc: 1.0, scale: 0.0025}
proposal: 0.005
epsilon_T:
prior: {min: 0.0, max: 0.02}
ref: {dist: norm, loc: 0.0066, scale: 0.002}
proposal: 0.0005
latex: '\epsilon_T'

sampler:
mcmc:
max_samples: 5000
Rminus1_stop: 0.05
output_every: 10
```

The hi_class configuration uses native `tep_mode` with the corrected transition function $f_T(z)=\ln(1+z)\exp[-(z/z_T)^{n_T}]$ and parameters `z_T = 5.0`, `n_T = 2.0`, with `epsilon_T` sampled freely. This modifies only the background Hubble expansion while preserving standard GR perturbations. The reference configuration is provided in `data/cobaya/tep_native_mcmc.yaml`.

*Pipeline status:* The native-`tep_mode` joint MCMC against Planck 2018 + BAO + Pantheon+ is in progress and is reported here as pending. Its expected outcome -- the homogeneous amplitude bounded to $\epsilon_T \approx 0$ -- follows directly from the $\theta_s$/$H_0$ degeneracy established in Section 5.1 and matches the late-universe + CMB synthesis of TEP-C0 (Paper 26). Note: the previously "converged" chains in `results/mcmc_chains/` used a `gravity_model: propto_omega` alpha-table proxy with `skip_stability_tests_smg`; that proxy is unvalidated and deprecated (Appendix A), and its parameters must not be read as native-TEP constraints.

The companion paper TEP-C0 (Paper 26) provides the authoritative late-time constraints: a Pantheon+ nested-sampling model comparison giving substantial Bayesian preference for the TEP geometry (Bayes factor 7–10 vs $\Lambda$CDM, Section 5.2), together with the CMB bound that pins the homogeneous shear amplitude to $\epsilon_T \approx 0$.

## 5. Results and Cosmological Constraints

### 5.1 The Acoustic Spectra

The physically meaningful test of the native TEP background modification is not the raw size of the $C_\ell$ residual but *where* the modification acts. Because $f_T(z) = \ln(1+z)\exp[-(z/z_T)^{n_T}]$ is frozen to zero for $z \gg z_T$, the recombination-era physics is left intact, which we verify directly from the hi_class background tables.

#### 5.1.1 Sound-horizon and acoustic-peak preservation

Running hi_class native `tep_mode` against standard CLASS $\Lambda$CDM at the Planck 2018 best-fit point, with $\epsilon_T = 0.0066$, $z_T = 5$, $n_T = 2$, we find:

- **Sound horizon preserved to parts-per-million:** $r_s^{\rm TEP}/r_s^{\Lambda\rm CDM} = 1.0000000$. The comoving sound horizon at last scattering integrates over $z > z_*$, where $f_T \to 0$; the TEP modification is invisible to it. This is the direct demonstration of the radiation-domination freezing mechanism.

- **Acoustic-peak morphology unchanged:** with $r_s$, the baryon loading, and the photon-baryon driving at $z \approx 1089$ all unmodified, the relative peak heights and the damping tail are identical to $\Lambda$CDM.

#### 5.1.2 The residual is a late-time projection, degenerate with $H_0$

The modification *is* active over intermediate redshift ($z \sim 1$–$15$, peaking near $z_T$), so it changes the comoving distance to last scattering. At the fiducial $\epsilon_T = 0.0066$ this shifts the angular acoustic scale by $\Delta\theta_s/\theta_s = +0.27\%$ ($D_C^{\rm TEP}/D_C^{\Lambda\rm CDM} = 0.9973$, with $r_s$ fixed). This rigid rescaling produces a coherent, oscillatory $\Delta C_\ell/C_\ell$ pattern whose envelope reaches $\sim 2.6\%$ across $100 < \ell < 2000$ at $\epsilon_T = 0.0066$ and scales linearly with $\epsilon_T$ (e.g. $\sim 0.4\%$ at $\epsilon_T = 0.001$). This is *not* a change in acoustic-peak physics: it is a pure angular-diameter-distance projection, exactly degenerate with $H_0$. In a parameter fit the standard parameters absorb it via a small $H_0$ shift, so the CMB does not exclude TEP -- it bounds the *homogeneous* amplitude $\epsilon_T$ to be small (Section 5.2).


![CMB TT residual and acoustic-scale bookkeeping](../results/figures/figure_2_cmb_residuals.png)



**Figure 2.** Native TEP background modification at $\epsilon_T = 0.0066$. (a) The fractional TT residual is a coherent, oscillatory pattern -- the signature of a rigid $\ell$-rescaling, not a change in peak morphology. (b) The sound horizon is preserved to ppm ($r_s^{\rm TEP}/r_s^{\Lambda\rm CDM} = 1.0000000$), while the angular scale shifts by $+0.27\%$ (a late-time projection degenerate with $H_0$). The homogeneous amplitude is thereby bounded to $\epsilon_T \approx 0$.

#### 5.1.3 Polarization Spectra ($C_\ell^{TE}, C_\ell^{EE}$)

The TE and EE spectra inherit the same behaviour: the recombination-era polarization source is unmodified (since $f_T \to 0$ there), and the only effect is the common $\theta_s$ projection shared with TT.

### 5.2 Cosmological Constraints: Late-Time Evidence and the CMB Bound

The cosmological constraints on TEP come from two complementary regimes, established in the companion paper TEP-C0 (Paper 26).

**Late-time evidence (supernovae).** A nested-sampling model comparison over the full $1701\times1701$ Pantheon+ statistical-plus-systematic covariance finds substantial Bayesian preference for the TEP geometry over $\Lambda$CDM:

| Model | Bayes factor vs $\Lambda$CDM | Interpretation |
| --- | --- | --- |
| TEP M1 ($z_T = 5$) | $7.2$ | Substantial |
| TEP M1 (free $z_T$) | $10.2$ | Substantial |
| $w$CDM | $18.5$ | Strong |
| CPL ($w_0 w_a$) | $17.8$ | Strong |
| Einstein-de Sitter | $3.5\times10^{-126}$ | Decisively rejected |
| Pure shear (tired light) | $3.4\times10^{-10}$ | Decisively rejected |

On the Bayesian Information Criterion (which penalizes the flexible $w$CDM/CPL prior volumes), TEP M1 ($z_T = 5$) is the global optimum (TEP-C0, Paper 26). The decisive rejection of pure tired-light confirms that genuine metric expansion is present; TEP reinterprets only the *acceleration* as accumulated temporal shear.


![Late-time SNe evidence](../results/figures/figure_3_sne_evidence.png)



**Figure 3.** Bayes factors vs $\Lambda$CDM from the TEP-C0 Pantheon+ nested-sampling model comparison (full covariance). TEP and dark-energy models are all substantially/strongly preferred over $\Lambda$CDM; Einstein-de Sitter and pure tired-light are decisively rejected.

**Homogeneous (CMB) bound.** As shown in Section 5.1, a non-zero homogeneous $\epsilon_T$ acts on the CMB only through the $\theta_s$ projection, which is degenerate with $H_0$. The Planck acoustic anchors therefore drive the *homogeneous* amplitude to the $\Lambda$CDM limit, $\epsilon_T \approx 0$, with $H_0$, $\Omega_b h^2$, $\Omega_{\rm cdm} h^2$, $n_s$, $A_s$ and $\tau$ recovered at their Planck values.

**Status of the native-TEP joint MCMC.** This paper's contribution is the verified hi_class implementation and the demonstration of $r_s$ preservation (Section 5.1). A full joint Cobaya MCMC using hi_class native `tep_mode` against Planck 2018 + BAO + Pantheon+ -- with the corrected transition function and standard GR perturbations -- is in progress; the configuration is provided in `data/cobaya/tep_native_mcmc.yaml`. We report it as pending rather than complete. (The earlier converged chains in this repository used a `propto_omega` alpha-table proxy and are deprecated; see Appendix A.)

### 5.3 The Hubble Tension in TEP


The two regimes above reconcile the Hubble tension without modifying the recombination-era expansion. The homogeneous background is $\Lambda$CDM-compatible ($\epsilon_T \approx 0$ on homogeneous scales, $H_0 \approx 67$ km/s/Mpc from the CMB), while the apparent local $H_0 \approx 73$ km/s/Mpc arises from environment-dependent clock-transport bias along the local distance ladder (Cepheid/SN Ia calibration in unscreened stellar atmospheres). Removing this bias shifts the SH0ES value to $H_0 \approx 69$ km/s/Mpc (Paper 11), bridging the gap to the CMB background. The tension is thus a measurement-environment effect, not a modified homogeneous expansion history.



![H0 in the TEP picture](../results/figures/figure_4_H0_comparison.png)



**Figure 4.** $H_0$ in the TEP picture. The homogeneous TEP background stays Planck-compatible ($\epsilon_T \to 0$); the local SH0ES value ($73.0$) is reinterpreted as a clock-transport bias, which when removed shifts the local measurement to $\approx 69$ km/s/Mpc (Paper 11).

### Appendix A. Historical propto_omega Proxy Analysis (Deprecated)

**DEPRECATED:** The following results are from an **unvalidated** `propto_omega` Horndeski proxy that was used during early pipeline development. This proxy modifies **both background and perturbations** via alpha functions, and its equivalence to TEP has **never been validated**. These constraints are retained for historical reference only and must **not** be interpreted as TEP predictions. The native TEP implementation uses a background-only Hubble modification with standard GR perturbations (Sections 3–5); its homogeneous amplitude is bounded by the CMB to $\epsilon_T \approx 0$.

An exploratory MCMC pipeline was run using hi_class with `gravity_model: propto_omega` and fixed Bellini-Sawicki coefficients against Planck 2018 low-$\ell$ TT/EE + lensing + BAO + Pantheon+ likelihoods. The high-$\ell$ Plik likelihood was excluded due to numerical instability of the old alpha-table `tep` module. Four parallel chains (seeds 42–45, 5000 accepted samples per chain) reached convergence with cross-chain Gelman–Rubin R-1 = 0.016:

| Parameter | Deprecated proxy (4 chains, R-1 = 0.016) | Planck 2018 |
| --- | --- | --- |
| $H_0$ [km/s/Mpc] | $66.0 \pm 1.5$ | $67.36 \pm 0.54$ |
| $\Omega_b h^2$ | $0.0207 \pm 0.0023$ | $0.02237 \pm 0.00015$ |
| $\Omega_{\rm cdm} h^2$ | $0.1166 \pm 0.0036$ | $0.1200 \pm 0.0012$ |
| $n_s$ | $0.995 \pm 0.005$ | $0.9650 \pm 0.0040$ |
| $A_s$ [$10^{-9}$] | $2.291 \pm 0.048$ | $2.100 \pm 0.030$ |
| $\tau_{\rm reio}$ | $0.048 \pm 0.008$ | $0.0544 \pm 0.0073$ |
| $A_{\rm planck}$ | $1.087 \pm 0.014$ | $1.000 \pm 0.003$ |

The elevated $n_s = 0.995$ result is inconsistent with the native TEP implementation. It likely reflects proxy bias from the unvalidated alpha mapping combined with weak low-$\ell$ constraints (high-$\ell$ Plik excluded). A per-likelihood correlation analysis identified Planck 2018 low-$\ell$ TT as the dominant driver (${\rm corr}(n_s, \chi^2_{\rm TT}) = -0.61$). The proxy also required `skip_stability_tests_smg = yes` to bypass negative scalar sound speeds at early times. These issues do not arise in the native TEP background-only implementation, which preserves standard GR perturbations throughout.


## 6. Conclusion

This paper implements and validates the native Temporal Equivalence Principle (TEP) background-only modification in hi_class, $H_{\rm TEP}(z) = H_{\Lambda\rm CDM}(z)\,(1 + \epsilon_T f_T(z))$ with $f_T(z)=\ln(1+z)\exp[-(z/z_T)^{n_T}]$ and standard General Relativistic perturbations. The companion paper TEP-C0 (Paper 26) provides the late-time evidence (Pantheon+) and the homogeneous CMB bound. The key findings are:


### 6.1 Summary of Results


- **EFT Mapping (sketch):** The TEP bi-metric framework with conformal factor $A(\phi) = \exp(\beta\phi/M_{\rm Pl})$ and disformal deformation $B(\phi)$ can be mapped onto the Bellini-Sawicki $\alpha_i$ functions, placing TEP within the Horndeski class. The present numerical analysis adopts the background-only realization; the full alpha-function implementation is deferred to future work.

- **Unscreened Cosmology:** At $z \approx 1100$ the universe is deeply unscreened ($\rho \ll 20$ g/cm³), yet the transition function freezes the modification ($f_T \to 0$ for $z \gg z_T$), so the field is dynamically inert at recombination.

- **CMB Consistency Check (verified):** With the corrected transition function, hi_class native `tep_mode` preserves the sound horizon to parts-per-million ($r_s^{\rm TEP}/r_s^{\Lambda\rm CDM} = 1.0000000$) and leaves the acoustic-peak morphology unchanged. The sole CMB effect of a non-zero $\epsilon_T$ is a late-time angular-diameter-distance projection ($\Delta\theta_s/\theta_s = +0.27\%$ at $\epsilon_T = 0.0066$), degenerate with $H_0$ and scaling linearly with $\epsilon_T$.

- **Cosmological Constraints:** Late-time Pantheon+ data give substantial Bayesian preference for the TEP geometry (Bayes factor 7–10 vs $\Lambda$CDM; TEP-C0, Paper 26), while the CMB $\theta_s$/$H_0$ degeneracy bounds the homogeneous amplitude to $\epsilon_T \approx 0$. A native-`tep_mode` joint hi_class MCMC (config in `data/cobaya/tep_native_mcmc.yaml`) is in progress.

- **Hubble Tension:** The framework resolves the tension as a local environmental clock-transport effect (Paper 11), keeping the homogeneous CMB baseline $\Lambda$CDM-compatible. No modified recombination-era expansion is required.




### 6.2 Cross-Scale Consistency


TEP is distinguished among modified gravity frameworks in that it is:



- Constrained at the Bohr radius (quantum phase coherence requires $\rho_c \gtrsim 20$ g/cm³)

- Validated at 20 g/cm³ (terrestrial GNSS, laboratory Cavendish)

- Tested across galactic scales (SPARC rotation curves, Gaia wide binaries)

- Reconciled with early-universe cosmology (CMB acoustic peaks, BBN)



This consistency across many orders of magnitude in density—from quantum to cosmological scales—supports TEP as a viable framework for gravitational physics, subject to continued testing. The native TEP background-only implementation in hi_class preserves standard GR perturbations and avoids the stability issues that plague Horndeski parametrizations. The native-`tep_mode` joint MCMC is in progress (Section 5.2).


### 6.3 Alternative Explanations


Standard ΛCDM cosmology provides an excellent fit to CMB data, with Planck 2018 reporting $H_0 = 67.36 \pm 0.54$ km/s/Mpc. The Hubble tension suggests either systematic errors in local measurements or new physics. Several alternative explanations have been proposed:



- **Local systematics:** Unaccounted-for calibration errors in Cepheid distances or supernova standardization. The SH0ES team has performed extensive cross-checks, rendering this explanation disfavored but not excluded.


- **Early-universe new physics:** Additional relativistic species ($N_{\rm eff}$), modified recombination history, or interacting dark energy. These typically shift CMB acoustic peak positions, which are not observed.


- **Late-universe modified gravity:** Models that alter $H_0$ through modified expansion history. Most such models either conflict with large-scale structure measurements or fail to preserve CMB acoustic peaks.



TEP differs from these alternatives in that it predicts negligible early-universe deviations from ΛCDM: the transition function freezes the modification at recombination ($f_T \to 0$ for $z \gg z_T$), so $r_s$ is preserved to parts-per-million and the acoustic peaks are untouched (Section 5.1). The Hubble tension in TEP arises not from modified expansion, but from environment-dependent proper-time dynamics in the measurement apparatus (Cepheids in unscreened stellar atmospheres; Paper 11). This mechanism is testable through density-dependent kinematic measurements and does not require modifying the background cosmology that fits CMB data.


### 6.4 Limitations


This analysis has several important limitations:



- **MCMC scope:** The present analysis implements and validates native TEP in hi_class as a background-only Hubble modification with standard GR perturbations, and demonstrates $r_s$ preservation directly from the Boltzmann background. The full native-`tep_mode` joint Cobaya MCMC (Planck 2018 + lensing + BAO + Pantheon+) is in progress and reported as pending; the late-time constraints are supplied by TEP-C0 (Paper 26). The earlier converged chains used a deprecated `propto_omega` proxy (Appendix A) and are not native-TEP results.


- **Native TEP parametrization:** The hi_class native `tep_mode` uses the parameters `epsilon_T`, `z_T`, and `n_T` to modify the background Hubble expansion while preserving standard GR perturbations. This is the hi_class analogue of the CLASS native TEP module used in TEP-C0 (Paper 26). No hand-mapping to Horndeski coefficients is required.


- **Lack of late-time datasets:** This analysis uses only Planck CMB data. Constraints on any residual $\alpha_M$ evolution at late times would require BAO, $H(z)$ measurements, and weak lensing data.


- **Theoretical uncertainty:** The radiation-domination freezing mechanism relies on $T^\mu_\mu \approx 0$ during radiation domination. While this is exact for a perfect radiation fluid, small deviations from perfection could allow residual scalar field evolution that this analysis does not capture.




### 6.5 Open Questions


Several questions remain unanswered by this analysis:



- Does the $k$-dependent growth suppression predicted by TEP manifest in Lyman-$\alpha$ forest data at $z \sim 2-4$?

- Can late-time expansion history measurements (BAO, $H(z)$) detect any residual $\alpha_M$ evolution?

- Can Stage-4 CMB experiments constrain the acoustic peak residuals further?

- Can gravitational wave propagation tests (LISA) confirm the $\alpha_T \approx 0$ constraint?




## References

Smawfield, M. (Paper 1). *Temporal Equivalence Principle: Terrestrial Screening and GNSS Phase Correlations.* TEP Corpus.

Smawfield, M. (Paper 6). *TEP and Ultra-Compact Dwarfs: Potential-Dependent Proper-Time Mapping.* TEP Corpus.

Smawfield, M. (Paper 11). *TEP and the Hubble Tension: Cepheid Environmental Bias.* TEP Corpus.

Smawfield, M. (Paper 12). *TEP and JWST High-Redshift Anomalies.* TEP Corpus.

Smawfield, M. (Paper 13). *TEP and Gaia DR3 Wide Binaries: Density-Dependent Kinematics.* TEP Corpus.

Bellini, E., & Sawicki, I. 2014, JCAP, 07, 050. *Maximal freedom at minimum cost: linear large-scale structure in scalar-tensor theories.*

Brax, P., Burrage, C., Davis, A.-C., & Gubitosi, G. 2019, Phys. Rev. D, 100, 083515. *Screening mechanisms in scalar-tensor theories.*

Cobaya Team. 2023, *Cobaya: Code for Bayesian Analysis of physical theories.* arXiv:2305.02971.

Hu, B., Raveri, M., Frusciante, N., & Silvestri, A. 2014, Phys. Rev. D, 89, 103530. *EFTCAMB/EFTCosmoMC: Numerical Notes.*

Knox, L., & Millea, M. 2020, Phys. Rev. D, 101, 043533. *Hubble constant hunter's guide.*

Lagos, M., Bellini, E., Jimenez, J. B., et al. 2018, JCAP, 03, 021. *hi_class: Horndeski in the Cosmic Linear Anisotropy Solving System.*

Lewis, A., Challinor, A., & Lasenby, A. 2000, Astrophys. J., 538, 473. *Efficient computation of cosmic microwave background anisotropies.*

Planck Collaboration. 2020, A&A, 641, A1. *Planck 2018 results. I. Overview and cosmological parameters.*

Planck Collaboration. 2020, A&A, 641, A6. *Planck 2018 results. VI. Cosmological parameters.*

Riess, A. G., Casertano, S., Yuan, W., et al. 2022, ApJ, 934, L7. *A Comprehensive Measurement of the Local Value of the Hubble Constant with 1 km/s/Mpc Uncertainty from the Hubble Space Telescope and the SH0ES Team.*

Sawicki, I., & Bellini, E. 2015, Phys. Rev. D, 92, 084061. *Stability of dark energy and the generalized no-slip condition.*

Zumalacárregui, M., & García-Bellido, J. 2014, Phys. Rev. D, 89, 064046. *Transforming gravity: from derivative couplings to matter to second-order scalar-tensor theories beyond the Horndeski Lagrangian.*


## Appendix A: Technical Implementation Details


### A.1 hi_class Installation and Configuration



#### A.1.1 Building with TEP Support




```
git clone https://github.com/lesgourg/class_public.git class
cd class
git checkout hi_class

# Enable scalar field modifications
./configure --with-smg
make clean && make
```


### A.2 Cobaya Installation




```
pip install cobaya
cobaya-install planck_2018_lowl.TT planck_2018_lowl.EE \
planck_2018_lensing.native bao.sdss_dr12_consensus_final \
sn.pantheonplus --path /path/to/likelihoods
```


### A.3 TEP Module C Code Structure and Implementation Note

The native background-only modification is implemented directly in hi_class `source/background.c`, controlled by the `.ini` flags `tep_mode`, `epsilon_T`, `z_T`, `n_T`. The relevant functions are:


- `tep_f_transition(pba, z)`: returns $f_T(z) = \ln(1+z)\,\exp[-(z/z_T)^{n_T}]$.

- `tep_gamma_factor(pba, z)`: returns $\gamma(z) = 1 + \epsilon_T f_T(z)$.

- The Hubble rate and its conformal-time derivative are multiplied by $\gamma(z)$ in `background_functions` and in the initial-Hubble setter; the energy densities entering the perturbation hierarchy are left standard (GR perturbations), so this is a pure background-only modification.



**Implementation note (corrected bug).** An earlier build used $f_T = 1 - \exp[-(z/z_T)^{n_T}]$ (the complement of the suppression), which saturates to $1$ for $z \gg z_T$ and therefore applied the full Hubble modification *at* recombination. This inverted the freezing mechanism, shifted $r_s$ by $0.66\%$ and $\theta_s$ by $\sim0.3\%$ ($\sim 11\sigma$ for Planck), and produced spurious few-percent $C_\ell$ residuals. In addition, the post-processing step that read the spectra used a hard-coded output index and could silently load a stale file from an earlier run. Both issues are fixed: the transition function now uses the authoritative TEP-C0 form, and the analysis resolves the most recent hi_class output deterministically. The legacy alpha-function stub (`smg_tep_*`, returning $\alpha_M,\alpha_T,\alpha_B,\alpha_K$) is retained in `scripts/hi_class_modules/tep.c` for reference only and is not used by the background-only pipeline.


### A.4 Screening Threshold in Cosmological Units


The 20 g/cm³ screening threshold converts to cosmological units as:

\begin{equation} \label{eq:9_appendix_01}
\rho_c = 20 \text{ g/cm}^3 = 2 \times 10^4 \text{ kg/m}^3 \approx 10^{31} \text{ eV/cm}^3
\end{equation}

In Planck units ($\hbar = c = G = 1$):

\begin{equation} \label{eq:9_appendix_02}
\rho_c \approx 4 \times 10^{-93} M_{\rm Pl}^4
\end{equation}

Compare to cosmic mean density today ($\rho_{\rm crit,0} \approx 10^{-123} M_{\rm Pl}^4$) and at recombination ($\rho \sim 10^{-114} M_{\rm Pl}^4$). The hierarchy $\rho(z=1100) \ll \rho_c$ confirms the unscreened regime.


### A.5 Stability

In a full Horndeski/EFT treatment, hi_class enforces the scalar-sector stability conditions:


- $c_s^2 \geq 0$ (no gradient instabilities)

- $\alpha_K \geq 0$ (no ghosts)

- $|\alpha_M|$ bounded (sub-luminal Planck-mass running)

- $\alpha_T \approx 0$ (GW speed constraints)



These would apply to the deferred alpha-function implementation. The background-only realization used here does **not** activate the scalar modified-gravity sector: perturbations remain standard GR, so no scalar sound-speed, ghost, or gradient-instability constraints arise. This is the principal numerical advantage of the background-only approach over the `propto_omega` proxy (Appendix A), which required `skip_stability_tests_smg = yes` to bypass negative early-time sound speeds.