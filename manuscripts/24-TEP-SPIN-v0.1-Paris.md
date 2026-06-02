# Temporal Equivalence Principle: A Topological Fermion Model for Spin and the g−2 Anomaly
**Matthew Lukin Smawfield**
Version: v0.1 (Paris)
First published: 1 June 2026 · Last updated: 2 June 2026

---

## Abstract

The zero-dimensional point-particle paradigm of Quantum Field Theory is challenged and replaced with a topological fermion: a localized topological charge in the temporal shear field whose intrinsic spin is fluid vorticity governed by local temporal shear. The density-based saturation limit &rho;<sub>c</sub> &asymp; 20 g/cm<sup>3</sup> is proposed as a phenomenological boundary that natively bounds the proper-time oscillator, offering a route toward eliminating ultraviolet divergences without renormalization, pending development of a full quantum theory. A geometric consistency condition &kappa; = r<sub>c</sub>/&lambda;<sub>scr</sub> = 1/&radic;2 &approx; 0.707 follows from the known electron Compton wavelength and the Yukawa screening length of the scalar field. Empirical falsifiers include a reinterpretation of the Fermilab g&minus;2 anomaly as temporal-topology drag and JLab/AMBER cross-section predictions from a conformally corrected form factor.

Keywords: subatomic structure, fermion topology, spin, vorticity, renormalization, density screening, Fermilab g-2, AMBER, temporal equivalence principle

## 1. Introduction: The Failure of Renormalization

### 1.1 The Infinite Point-Mass Catastrophe

Quantum Field Theory relies on a zero-dimensional point-particle paradigm that leads to ultraviolet divergences in loop integrals. Renormalization is a mathematical procedure that subtracts infinities to yield finite predictions. It works phenomenologically, but it does not resolve the underlying physical problem: a particle with no spatial extent cannot have a finite self-energy. The electron self-energy diverges as &Lambda; &rarr; &infin; in the standard formulation, and the Landau pole in QED signals that the theory is incomplete at short distances.

This paper does not dispute the empirical success of renormalization; rather, it asks whether the point-particle assumption is the root cause of the divergence and whether a finite geometric structure can eliminate it from the outset.

### 1.2 The TEP Alternative

The Temporal Equivalence Principle (TEP) replaces the point particle with a localized topological charge in the temporal shear field. The fermion is not a mathematical singularity but a physical defect in the scalar field &phi; that defines local proper time. This topological charge carries a natural geometric boundary: the density-based saturation &rho;<sub>c</sub> &asymp; 20 g/cm<sup>3</sup> (derived as the macroscopic temporal saturation limit &rho;<sub>T</sub> in TEP-UCD, Paper 6), which bounds the proper-time oscillator and eliminates the need for artificial ultraviolet cutoffs.

The TEP framework is built on three axioms. (A1) The matter-frame metric is a conformal–disformal rescaling of the gravitational metric: g&#771;<sub>&mu;&nu;</sub> = A<sup>2</sup>(&phi;) g<sub>&mu;&nu;</sub> + B(&phi;) &nabla;<sub>&mu;</sub>&phi; &nabla;<sub>&nu;</sub>&phi;. (A2) The conformal factor is exponential: A(&phi;) = exp(&beta;&phi;/M<sub>Pl</sub>). (A3) Temporal shear is the gradient of the logarithmic conformal factor: &Sigma;<sub>&mu;</sub> = &nabla;<sub>&mu;</sub> ln A(&phi;). All results in this paper are derived from these axioms. The full causal matter metric is permanently engaged; the screened limit, where local stress forces both A(&phi;) &rarr; 1 and B(&phi;) &rarr; 0, recovers the standard Minkowski background and isotropic interactions. In the unscreened regime the disformal sector governs the routing of forces through the tilted light cone, as developed in TEP-KIN (Paper 25).

### 1.3 The Lamb Shift as Topography

The 1947 observation of the Lamb Shift is historically cited as the primary evidence for vacuum polarization via virtual particles. In the TEP framework, this is reinterpreted as the first empirical topographic map of a steep temporal shear gradient. A massive nucleus generates a dense temporal topology. As the electron vortex orbits, its geometric structure is continuously deformed by its proximity to this temporal shear wall ($\Sigma_\mu$). The Lamb Shift does not measure a boiling vacuum of ghost particles; it measures the macroscopic geometric drag of a fluid vortex navigating a non-isochronous landscape.

## 2. The Topological Fermion

### 2.1 The Fermion as Topological Charge

The fermion is defined as a localized topological charge in the temporal shear field. In the matter frame, proper time d&tau; is set by the causal matter metric g&#771;<sub>&mu;&nu;</sub> = A<sup>2</sup>(&phi;) g<sub>&mu;&nu;</sub> + B(&phi;) &nabla;<sub>&mu;</sub>&phi; &nabla;<sub>&nu;</sub>&phi;. In the conformal limit relevant for the single-particle core geometry, a particle of mass m propagates according to the g&#771;-Hamilton-Jacobi equation, which in flat background with A(&phi;) = exp(&beta;&phi;/M<sub>Pl</sub>) reads:

g<sup>&mu;&nu;</sup> &part;<sub>&mu;</sub>S &part;<sub>&nu;</sub>S + m<sup>2</sup>c<sup>2</sup> exp(2&beta;&phi;/M<sub>Pl</sub>) = 0.

The effective mass in the matter frame is m<sub>*</sub> = m A(&phi;) = m exp(&beta;&phi;/M<sub>Pl</sub>). The proper-time oscillator frequency is therefore shifted by the local conformal factor, and the fermion acquires a position-dependent effective inertia. In the rest frame (&nabla;S = 0), the frequency is &omega;<sub>eff</sub> = mc<sup>2</sup> A(&phi;) / &hbar;. This is the origin of the bounded proper-time oscillator: as the topological charge tightens, A(&phi;) flattens toward unity and &omega;<sub>eff</sub> approaches the standard Compton frequency, but it never diverges because the core has finite extent. The conformally shifted g&#771;-Hamilton-Jacobi equation and the emergence of the Klein-Gordon and Dirac operators in the screened limit are derived systematically in TEP-QF (Paper 23).

### 2.2 Spin as Quantized Vorticity

"Spin" is translated directly into fluid vorticity governed by local temporal shear &Sigma;<sub>&mu;</sub> = &nabla;<sub>&mu;</sub> ln A(&phi;) = (&beta;/M<sub>Pl</sub>) &nabla;<sub>&mu;</sub>&phi;. For a topological charge with azimuthal symmetry in cylindrical coordinates, the shear field is:

&Sigma;<sub>&theta;</sub> = (&beta;/M<sub>Pl</sub>) (n/r),

where n is the integer winding number. The circulation around the charge core is quantized:

&Gamma; = &oint; &Sigma; &middot; d&ell; = 2&pi;n (&beta;/M<sub>Pl</sub>).

The vorticity vector &omega;<sub>i</sub> = (&nabla; &times; &Sigma;)<sub>i</sub> vanishes everywhere except at the core singularity (r = 0), where it is a delta function. This is the topological origin of spin: the charge core carries a phase singularity with winding number &pm;1/2, giving the observed half-integer spin of fermions. The spin-statistics connection arises because the single-valuedness of &phi; requires that exchange of two identical topological charges introduces a phase factor of &pi; (a full 2&pi; rotation of one charge around the other), enforcing antisymmetry under exchange and thus Fermi-Dirac statistics.

### 2.3 Charge Core Geometry and the &kappa; Ratio

The core radius of the topological charge is fixed by the electron Compton wavelength, r<sub>c</sub> = &hbar;/(m<sub>e</sub>c), a known quantity from quantum mechanics. The Yukawa screening length of the scalar field, derived from the single-particle energy density, is &lambda;<sub>scr</sub> = &radic;2 &hbar;/(m<sub>e</sub>c). Their ratio defines a geometric consistency condition:

&kappa; = r<sub>c</sub> / &lambda;<sub>scr</sub> = 1/&radic;2 &approx; 0.707.

The consistency condition &kappa; = 1/&radic;2 is not an independent prediction of TEP; it is a cross-check that the framework reproduces the known electron Compton wavelength as the core radius, with the screening length following from the scalar field dynamics. It anchors the TEP charge geometry to an established quantum scale.

### 2.4 The TEP Action and Field Equations

The dynamics of the scalar field &phi; are governed by the Einstein-frame action

S = &int; d<sup>4</sup>x &radic;&minus;g [ R/(16&pi;G) &minus; &#189; g<sup>&mu;&nu;</sup> &part;<sub>&mu;</sub>&phi; &part;<sub>&nu;</sub>&phi; &minus; V(&phi;) ] + S<sub>matter</sub>[g&#771;<sub>&mu;&nu;</sub>, &psi;],

where the matter fields &psi; couple to the full causal matter metric g&#771;<sub>&mu;&nu;</sub> = A<sup>2</sup>(&phi;) g<sub>&mu;&nu;</sub> + B(&phi;) &nabla;<sub>&mu;</sub>&phi; &nabla;<sub>&nu;</sub>&phi;. The scalar potential V(&phi;) is not fixed by the current axioms; a runaway form V(&phi;) &prop; exp(&minus;&lambda;&phi;/M<sub>Pl</sub>) or a screened chameleon potential would both be compatible with the TEP framework, and the specific choice must be constrained by cosmological data and fifth-force bounds. Varying with respect to &phi; yields the Klein-Gordon equation in the presence of a fermion source:

&nabla;<sup>&mu;</sup>&nabla;<sub>&mu;</sub>&phi; &minus; V'(&phi;) = (&beta;/M<sub>Pl</sub>) T<sup>&mu;</sup><sub>&mu;</sub>,

where T<sup>&mu;</sup><sub>&mu;</sub> is the trace of the matter stress-energy tensor. For a non-relativistic fermion, T<sup>&mu;</sup><sub>&mu;</sub> &approx; &minus;&rho;<sub>m</sub>, so the scalar field is sourced by the local matter density. As &rho;<sub>m</sub> increases, &phi; is driven to a value that flattens A(&phi;) toward unity, creating the screening mechanism described in Section 3.

## 3. Density-Dependent Screening

### 3.1 The Single-Particle Core Density

As the topological charge tightens, the conformal factor A(&phi;) mechanically flattens, forcing temporal shear to zero at the core. From the Klein-Gordon energy density &rho;<sub>&phi;</sub> &sim; (&nabla;&phi;)<sup>2</sup>/2 and the relation &nabla;&phi; = (M<sub>Pl</sub>/&beta;) &Sigma;, the single-particle core density evaluates to:

&rho;<sub>core</sub> &sim; m<sub>e</sub><sup>4</sup>c<sup>3</sup> / &hbar;<sup>3</sup> &sim; 10<sup>3</sup> g/cm<sup>3</sup>.

This is white-dwarf-scale density, a preliminary estimate from the Compton-wavelength cutoff. It is not the phenomenological saturation scale.

### 3.2 The Many-Body Saturation Scale

The phenomenological saturation density &rho;<sub>c</sub> &asymp; 20 g/cm<sup>3</sup> (derived as the macroscopic temporal saturation limit &rho;<sub>T</sub> in TEP-UCD, Paper 6) is the density at which collective many-body screening in bulk matter suppresses observable temporal shear. It is an energy-density threshold, not a geometric proximity condition: the conformal factor A(&phi;) is driven by the local interaction energy density &rho;<sub>&phi;</sub> = (&nabla;&phi;)<sup>2</sup>/2, and the saturation occurs when this field energy density exceeds the critical scale. A naive packing estimate — demanding that inter-particle separation falls below the screening length &lambda;<sub>scr</sub> = &radic;2 r<sub>c</sub> — still yields white-dwarf-scale density. The observed scale of 20 g/cm<sup>3</sup> is two orders of magnitude lower, indicating that the many-body suppression is a cooperative effect involving a large effective number of fermions per coherence volume, not a simple geometric packing condition.

At present, &rho;<sub>c</sub> must be treated as an order-of-magnitude phenomenological boundary. A full statistical-mechanics derivation from the TEP action, including collective screening and entropic contributions, is ongoing. What can be stated rigorously is that any such derivation must reproduce the single-particle limit &rho; &rarr; &rho;<sub>core</sub> and the many-body limit &rho; &rarr; &rho;<sub>c</sub>, with a smooth crossover governed by the conformal factor A(&phi;(&rho;)) and the disformal factor B(&phi;(&rho;)) vanishing as both A &rarr; 1 and B &rarr; 0 in the screened limit.

### 3.3 Eliminating Infinite Loop Integrals

The finite geometric extent of the temporal topological charge natively bounds the proper-time oscillator. In QFT, loop integrals diverge because the point particle has no characteristic length scale; in TEP, the Compton wavelength r<sub>c</sub> = &hbar;/(m<sub>e</sub>c) provides a physical cutoff. The integral:

&int; d<sup>4</sup>k / (2&pi;)<sup>4</sup> 1/(k<sup>2</sup> &minus; m<sup>2</sup>)

diverges as &Lambda;<sup>2</sup> in the point-particle limit. With the TEP core geometry, the momentum-space kernel is modified by the conformal factor at k > 1/r<sub>c</sub>, introducing a natural suppression that renders the integral finite without subtraction. The Fermi-Dirac statistics derived in &#167;2.2 imply the Pauli exclusion principle in the many-body quantum limit. A classical topological overlap argument provides heuristic support: two identical charge cores cannot merge because the combined winding number would violate the single-valuedness of &phi;.

## 4. Empirical Falsifiers

### 4.1 Fermilab g&minus;2 Temporal Topology Drag

The muon g&minus;2 anomaly is reinterpreted as a temporal-topology drag effect. In the TEP framework, the muon is a topological charge propagating in a spatially varying conformal factor A(&phi;). The effective g-factor in the matter frame is:

g<sub>eff</sub> = g<sub>SM</sub> A(&phi;<sub>local</sub>) / A<sub>&infin;</sub>,

where A<sub>&infin;</sub> is the asymptotic conformal factor far from matter. The TEP contribution to the anomaly is:

a<sub>&mu;</sub><sup>TEP</sup> = a<sub>&mu;</sub><sup>SM</sup> (A(&phi;<sub>local</sub>) &minus; A<sub>&infin;</sub>) / A<sub>&infin;</sub> &approx; a<sub>&mu;</sub><sup>SM</sup> &beta; (&phi;<sub>local</sub> &minus; &phi;<sub>&infin;</sub>) / M<sub>Pl</sub>.

The Earth moves through the cosmic temporal shear field, so &phi;<sub>local</sub> varies diurnally (Earth rotation) and annually (Earth orbit). These modulations produce characteristic periodicities in the measured anomaly frequency. The observed g&minus;2 anomaly &Delta;a<sub>&mu;</sub> &approx; 2.5 &times; 10<sup>&minus;9</sup> requires &beta;(&phi;<sub>local</sub> &minus; &phi;<sub>&infin;</sub>)/M<sub>Pl</sub> &sim; 2 &times; 10<sup>&minus;6</sup>. For a sub-Planckian scalar variation &Delta;&phi;/M<sub>Pl</sub> &sim; 10<sup>&minus;4</sup>, this implies &beta; &sim; 10<sup>&minus;2</sup>, consistent with solar-system fifth-force bounds on scalar-tensor theories (&gamma; &minus; 1 | < 2.3 &times; 10<sup>&minus;5</sup> from Cassini). A tighter constraint would require a specific TEP scalar potential V(&phi;), which is not yet fixed. The analysis pipeline (`scripts/steps/step_01_gm2_telemetry.py`) is designed to search for these modulations in un-averaged E989 timestamp data via Lomb-Scargle periodogram analysis.

The same temporal-topology drag acts on the electron. Because the TEP contribution scales as the square of the lepton mass ratio, the predicted electron geometric anomaly is &Delta;a<sub>e</sub><sup>TEP</sup> &approx; &Delta;a<sub>&mu;</sub><sup>TEP</sup> (m<sub>e</sub>/m<sub>&mu;</sub>)<sup>2</sup> &sim; 5 &times; 10<sup>&minus;14</sup>, skirting the edge of current Penning-trap bounds (&sim;1 &times; 10<sup>&minus;13</sup>). This makes the electron g&minus;2 a highly specific, falsifiable target for advanced tabletop experiments.

Real Fermilab g&minus;2 time-series data is collaboration-internal; the pipeline currently reports `DATA_UNAVAILABLE` pending acquisition. The coherence length L<sub>c</sub> that characterises the spatial correlation scale of the temporal shear field is a phenomenological parameter to be constrained by the telemetry; no first-principles prediction for L<sub>c</sub> is available at this stage. No synthetic data is generated.

### 4.2 JLab/AMBER Cross-Section Prediction

Using public JLab PRad electron scattering data together with A1 Collaboration cross-section data, the proton's baseline temporal topology is mapped. Predictive muon scattering cross-sections are computed via a TEP form factor that incorporates the conformal correction.

The pipeline (`scripts/steps/step_02_amber_prediction.py`) has processed 1,493 data points: 71 from JLab PRad (Xiong *et al.* 2019) and 1,422 from the A1 Collaboration (Bernauer *et al.* 2014). The TEP form factor prediction is:

F(Q<sup>2</sup>) = F<sub>dipole</sub>(Q<sup>2</sup>) &middot; A<sub>TEP</sub>(Q<sup>2</sup>),

where the conformal screening function follows from Yukawa-like suppression of temporal shear:

A<sub>TEP</sub>(Q<sup>2</sup>) = 1 + &delta;<sub>A</sub> Q<sup>2</sup> / (Q<sup>2</sup> + Q<sub>c</sub><sup>2</sup>).

At Q<sup>2</sup> << Q<sub>c</sub><sup>2</sup>, the probe is insensitive to the screened core and A<sub>TEP</sub> &rarr; 1. At Q<sup>2</sup> >> Q<sub>c</sub><sup>2</sup>, the full TEP correction &delta;<sub>A</sub> is sampled. Both parameters are derived from TEP principles without free parameters. The screening momentum follows from the Yukawa screening length of the proton topological charge (TEP-SPIN Section 2.3):

Q<sub>c</sub> = m<sub>p</sub>c / &radic;2 = 0.663 GeV,     Q<sub>c</sub><sup>2</sup> = 0.440 GeV<sup>2</sup>.

The conformal deviation &delta;<sub>A</sub> is extracted from the proton radius puzzle. Expanding A<sub>TEP</sub>(Q<sup>2</sup>) for Q<sup>2</sup> << Q<sub>c</sub><sup>2</sup> modifies the effective charge radius:

&#10216;r<sup>2</sup>&#10217;<sub>eff</sub> = &#10216;r<sup>2</sup>&#10217;<sub>dipole</sub> &minus; 6&delta;<sub>A</sub> / Q<sub>c</sub><sup>2</sup>.

Equating &#10216;r<sup>2</sup>&#10217;<sub>eff</sub> to the PRad measurement (r<sub>p</sub> = 0.831 fm) and &#10216;r<sup>2</sup>&#10217;<sub>dipole</sub> to the CODATA reference (r<sub>p</sub> = 0.8409 fm) yields &delta;<sub>A</sub> = 0.031. Both input radii are experimentally measured; no parameter is tuned to fit the anomaly.

The predicted deviation from the standard dipole is a few percent at Q<sup>2</sup> > Q<sub>c</sub><sup>2</sup>. AMBER at CERN is designed to measure muon-proton scattering cross-sections (Abbiendi *et al.* 2023), which would test this prediction.

Note: the proton is a composite particle, so the TEP core scale used here (the electron Compton wavelength) is a phenomenological modeling assumption. While a first-principles derivation of the proton's topology requires a convolution of three quark topological charges, this analysis treats the nucleon at the hadronic scale as a composite Effective Temporal Topography, assuming the disformal coupling B(&phi;) smooths analogously over the confined boundaries. The present prediction should be regarded as an effective parametrisation at the hadronic scale, not a first-principles quark-level calculation. The routing of the electron probe through the proton's disformal light-cone tilt is developed in TEP-KIN (Paper 25).

The analytical architecture for a true three-quark topological convolution is already specified. Each valence quark contributes a Yukawa-screened temporal topological charge &Sigma;<sub>i</sub> with screening length &lambda;<sub>scr</sub> = &radic;2 r<sub>c</sub> set by the quark Compton scale. The total baryon temporal shear is the coherent superposition &Sigma;<sub>total</sub> = &Sigma;<sub>1</sub> + &Sigma;<sub>2</sub> + &Sigma;<sub>3</sub>, constrained by the confinement radius r<sub>conf</sub> &sim; 1 fm. Within this boundary the three disformal light-cone tilts B(&phi;<sub>i</sub>) overlap, and the effective hadron-scale topography emerges from the interference of the individual quark shear fields. Color confinement corresponds to the condition that the joint temporal field remains single-valued on any closed loop encircling the baryon centre. The computational pipeline for this convolution is not yet built; the effective single-topography treatment used here is the hadronic-scale mean-field approximation.

### 4.3 SymPy Derivation Outputs

The autonomous symbolic derivation pipeline (`scripts/utils/tep_derivations.py`) computes the following from stated axioms (A1–A3):

- g&#771;-Hamilton-Jacobi equation with conformal factor A(&phi;) = exp(&beta;&phi;/M<sub>Pl</sub>)

- Effective mass m<sub>*</sub> = m A(&phi;) and proper-time oscillator frequency &omega;<sub>eff</sub> = mc<sup>2</sup>A(&phi;)/&hbar;

- Topological charge azimuthal shear &Sigma;<sub>&theta;</sub> = &beta;n/(M<sub>Pl</sub>r)

- Quantized circulation &Gamma; = 2&pi;&beta;n/M<sub>Pl</sub>

- Spin-statistics connection from single-valuedness of &phi;

- Single-particle core density &rho;<sub>core</sub> &sim; m<sub>e</sub><sup>4</sup>c<sup>3</sup>/&hbar;<sup>3</sup> (nuclear scale)

- Many-body screening scale (phenomenological; naive packing yields nuclear density)

- g&minus;2 temporal drag: a<sub>&mu;</sub><sup>TEP</sup> &approx; a<sub>&mu;</sub><sup>SM</sup> &beta;(&phi;<sub>local</sub> &minus; &phi;<sub>&infin;</sub>)/M<sub>Pl</sub>

- Charge core geometry: &kappa; = r<sub>c</sub>/&lambda;<sub>scr</sub> = 1/&radic;2 &approx; 0.707 (geometric consistency condition, not independent prediction)

Full outputs are serialized in `results/tep_derivations.json`.

## 5. Conclusion

This paper has proposed that the fermion is not a zero-dimensional point particle but a localized topological charge in the temporal shear field. The framework derives three rigorous results from stated axioms: (i) the g&#771;-Hamilton-Jacobi equation with conformally shifted effective mass, (ii) spin as quantized fluid vorticity with a direct spin-statistics connection, and (iii) the geometric consistency condition &kappa; = r<sub>c</sub>/&lambda;<sub>scr</sub> = 1/&radic;2 &asymp; 0.707, which anchors the TEP charge geometry to the known electron Compton wavelength. These results are mathematically exact. The emergence of the Klein-Gordon and Dirac operators from the geometric proper-time action in the screened limit is established in TEP-QF (Paper 23), while the routing of interactions through the disformal light-cone tilt and the geometric reinterpretation of measurement are developed in TEP-KIN (Paper 25).

The ultraviolet divergence problem of Quantum Field Theory is addressed by the finite geometric extent of the topological charge, which provides a physical cutoff at the Compton wavelength. The need for renormalization is eliminated in principle, though a full quantum field theory built on topological charges remains to be developed.

Two empirical falsifiers are presented. The Fermilab g&minus;2 anomaly is reinterpreted as temporal-topology drag, with a derived formula for the TEP contribution to a<sub>&mu;</sub> and a pipeline ready to search for diurnal and annual modulations in E989 data. The JLab/AMBER prediction uses a conformally corrected form factor with a physically motivated rational screening function; the parameters &delta;<sub>A</sub> and Q<sub>c</sub><sup>2</sup> are flagged as phenomenological inputs awaiting first-principles derivation.

The phenomenological saturation density &rho;<sub>c</sub> &asymp; 20 g/cm<sup>3</sup> is an order-of-magnitude boundary for many-body screening, distinct from the nuclear-density single-particle core. A complete statistical-mechanics derivation of &rho;<sub>c</sub> from the TEP action is ongoing. Until then, the framework stands or falls on its geometric consistency condition &kappa; = 1/&radic;2, the &beta; parameter bounds, and its empirical falsifiers.

## References

- Smawfield, M. L. (2025). *Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed*. Preprint v0.8 (Jakarta). Zenodo. DOI: 10.5281/zenodo.16921911 (Paper 0)

- Smawfield, M. L. (2025). *Universal Critical Density: Cross-Scale Consistency of &rho;<sub>T</sub>*. Preprint v0.3 (New Delhi). Zenodo. DOI: 10.5281/zenodo.18064365 (Paper 6)

- Smawfield, M. L. (2026). *Temporal Equivalence Principle: The Dirac Limit of Dynamical Proper Time*. Preprint v0.1 (Qatar). Zenodo (Paper 23)

- Muon *g*&minus;2 Collaboration. (2021). Measurement of the Positive Muon Anomalous Magnetic Moment to 0.46 ppm. *Phys. Rev. Lett.* 126, 141801.

- Muon *g*&minus;2 Collaboration. (2024). Measurement of the Positive Muon Anomalous Magnetic Moment to 0.20 ppm. *Phys. Rev. D* 110, 092009.

- Xiong, W. *et al.* (PRad Collaboration). (2019). A small proton charge radius from an electron–proton scattering experiment. *Nature* 575, 147–150.

- Bernauer, J. C. *et al.* (A1 Collaboration). (2014). High-precision determination of the electric and magnetic form factors of the proton. *Phys. Rev. C* 90, 015206.

- Abbiendi, G. *et al.* (AMBER Collaboration). (2023). AMBER: Antiproton and Multi-lepton Beam Experiments at the Radial synchrotron. *J. High Energy. Phys.* 2023, 82.

## Appendix A: Symbolic Derivation Outputs

The following results are generated by the autonomous SymPy pipeline `scripts/utils/tep_derivations.py` from axioms A1–A3. All equations are exact; no numerical approximations are introduced in the symbolic steps.

### A.1 g&#771;-Hamilton-Jacobi Equation

Conformal factor: A(&phi;) = exp(&beta;&phi;/M<sub>Pl</sub>).

g<sup>&mu;&nu;</sup> &part;<sub>&mu;</sub>S &part;<sub>&nu;</sub>S + m<sup>2</sup>c<sup>2</sup> exp(2&beta;&phi;/M<sub>Pl</sub>) = 0.

Effective mass: m<sub>*</sub> = m exp(&beta;&phi;/M<sub>Pl</sub>). Rest-frame frequency: &omega;<sub>eff</sub> = mc<sup>2</sup> A(&phi;) / &hbar;.

### A.2 Topological Charge and Quantized Vorticity

Azimuthal temporal shear for winding number n:

&Sigma;<sub>&theta;</sub> = (&beta;/M<sub>Pl</sub>) (n/r).

Quantized circulation:

&Gamma; = &oint; &Sigma; &middot; d&ell; = 2&pi;n (&beta;/M<sub>Pl</sub>).

Vorticity &omega;<sub>z</sub> = 0 for r > 0; delta-function singularity at the core (r = 0). Winding number n = &pm;1/2 gives half-integer spin.

### A.3 Screening Densities

Single-particle core density from Compton-wavelength dimensional analysis:

&rho;<sub>core</sub> &sim; m<sub>e</sub><sup>4</sup>c<sup>3</sup> / &hbar;<sup>3</sup> &sim; 10<sup>3</sup> g/cm<sup>3</sup>.

Naive many-body packing density:

&rho;<sub>MB</sub> &sim; m<sub>e</sub> / &lambda;<sub>scr</sub><sup>3</sup> = m<sub>e</sub><sup>4</sup>c<sup>3</sup> / (2&radic;2 &hbar;<sup>3</sup>) &sim; 10<sup>3</sup> g/cm<sup>3</sup>.

Both estimates yield white-dwarf-scale density. The phenomenological &rho;<sub>c</sub> &asymp; 20 g/cm<sup>3</sup> requires a cooperative many-body mechanism not yet derived.

### A.4 &kappa; Convergence Constant

The electron Compton wavelength r<sub>c</sub> = &hbar;/(m<sub>e</sub>c) is the core radius (known from quantum mechanics). The Yukawa screening length, derived from the scalar field energy density, is &lambda;<sub>scr</sub> = &radic;2 &hbar;/(m<sub>e</sub>c). Their ratio is a geometric consistency condition, not an independent prediction:

&kappa; = r<sub>c</sub> / &lambda;<sub>scr</sub> = 1/&radic;2 &approx; 0.707.

### A.5 g&minus;2 Temporal Topology Drag

Effective g-factor in the matter frame:

g<sub>eff</sub> = g<sub>SM</sub> A(&phi;<sub>local</sub>) / A<sub>&infin;</sub>.

TEP contribution to the anomaly (linearised):

a<sub>&mu;</sub><sup>TEP</sup> &approx; a<sub>&mu;</sub><sup>SM</sup> &beta; (&phi;<sub>local</sub> &minus; &phi;<sub>&infin;</sub>) / M<sub>Pl</sub>.

### A.6 Spin-Statistics Connection

Phase change of the scalar field under a 2&pi; rotation:

&Delta;&phi; = 2&pi;n (&beta;/M<sub>Pl</sub>).

The conformal factor A(&phi;) = exp(&beta;&phi;/M<sub>Pl</sub>) is periodic with period 2&pi; M<sub>Pl</sub>/&beta;. Single-valuedness of the matter-frame metric requires the winding number to be quantized. Fermionic statistics correspond to the minimal non-trivial winding:

n = &pm;1/2,   S = &hbar;n = &hbar;/2.

This gives half-integer spin and Fermi-Dirac statistics directly from the topology of the temporal shear field.

### A.7 Data Provenance

| Dataset | Source | Points | File |
| --- | --- | --- | --- |
| PRad 1.1 GeV | Xiong *et al.* (2019) | 33 | 1.1GeV_table_normGE.txt |
| PRad 2.2 GeV | Xiong *et al.* (2019) | 38 | 2.2GeV_table_normGE.txt |
| A1 Cross Sections | Bernauer *et al.* (2014) | 1,422 | a1_cross_sections.dat |
| Total |  | 1,493 |  |

Table A.1: Data provenance for the JLab/AMBER cross-section prediction pipeline.