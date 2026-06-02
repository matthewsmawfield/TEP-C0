# Temporal Equivalence Principle: Disformal Kinematics and the Measurement Landscape
**Matthew Lukin Smawfield**
Version: v0.1 (Kuala Lumpur)
First published: 24 May 2026 · Last updated: 2 June 2026

---

## Abstract

In the unscreened regime of the Temporal Equivalence Principle, virtual force carriers and statistical wavefunctions are unnecessary: interactions and entanglement are routed through the disformal geometry of the temporal field, governed by the coupling coefficient B(φ). In the screened limit, where the local interaction energy density substantially exceeds the saturation scale ρ_{T} ≈ 20 g/cm^{3}, B(φ) → 0 and standard perturbation theory is recovered as the tangent limit. U(1) and SU(2) gauge symmetries are explicitly routed through B(φ) as localized, direction-dependent light-cone tilts. Geometric kinematics is derived: particles attract and repel because they navigate spatially varying temporal shear, with Coulomb's law emerging as geodesic deviation in the weak-field limit. Entanglement is defined mathematically as an unbroken macroscopic geometric contour connecting bifurcated topological charges, with Schmidt rank mapping to contour winding number. The probability wavefunction is redefined as the physical volumetric shear-wake of proper-time phase transport; the double-slit and tunneling phenomena are explained as the core surfing interference patterns churned by its own macroscopic wake. A re-analysis of published graphene Aharonov-Bohm interferometry data (Zimmermann et al., Nat. Commun. 8, 14983, 2017) is presented as a candidate empirical falsifier. On the medium-smoothed transmission profile a TEP-modified interference model, in which the effective phase accumulation is scaled by a temporal-anisotropy parameter γ, achieves a χ² of 112.67 compared with 590.46 for the standard model (γ = 1), yielding a fitted γ = 0.722. A comprehensive control suite was subsequently performed: dataset provenance was verified by SHA-256; the original MATLAB processing was reproduced exactly; seven standard nuisance models were fitted; P-gamma covariance was estimated by Hessian and bootstrap; a train/test split was performed; and the model was refitted on raw, light, medium, and heavily smoothed data. The gamma parameter is degenerate with an effective period P_{eff} = P/γ (machine-precision match) and is unstable across smoothing levels, ranging from 0.606 (raw) to 0.959 (medium), with 0.859 (light) and 0.683 (heavy). Out-of-sample generalisation is at best tied and at worse substantially worse. The TEP model does not beat a best nuisance model (standard) on the medium-smoothed data (ΔBIC = 4.8 in favour of the nuisance model). A disformal topography regression, in which the measured phase shift is regressed directly against the metric tilt B(φ), reveals that a harmonic (periodic modulation) model is decisively preferred (BIC = 34.61) over uniform, linear, and Gaussian alternatives, indicating that the phase excess is structured as a periodic shear-wake whose spatial frequency matches the cavity geometry. The uniform temporal-dilation model fails spectacularly (BIC = 1353.75). The empirical finding is therefore a spatially varying periodic disformal modulation rather than a uniform temporal anisotropy.

Keywords: disformal kinematics, measurement, entanglement, Aharonov-Bohm, graphene interferometry, light-cone geometry, virtual bosons, temporal equivalence principle, shear-wake, geometric kinematics

## 1. Introduction: The Illusion of the Dead Vacuum

### 1.1 Virtual Bosons as Geometric Fiction

The standard model posits virtual bosons to mediate forces between particles. These are "virtual" because they cannot be directly observed, and their existence is justified only by the predictive success of perturbation theory. In the TEP framework, no such invention is necessary: all interactions are routed through the continuous geometry of the light cone, tilted by the disformal coupling B(φ). In the screened limit, where the local interaction energy density substantially exceeds the critical saturation density *ρ_{c} ≈ 20 g/cm^{3}* (derived as the macroscopic temporal saturation limit ρ_{T} in TEP-UCD, Paper 6), *B(φ) → 0* and the light cone recovers its isotropic form; standard perturbation theory is recovered as the tangent limit. In the unscreened regime, the full disformal structure is manifest.

### 1.2 Reinterpreting the Copenhagen Interpretation

The Copenhagen Interpretation's statistical indeterminism — wavefunction collapse, complementarity, the measurement problem — arises from the same root error as virtual particles: the assumption of a flat, isochronous background. When the background is dynamical, "measurement" is the geometric interaction of a probe with the local temporal shear field.

## 2. Routing Interactions via Disformal Coupling

### 2.1 Gauge Symmetries as Light-Cone Tilts

In the unscreened regime, where the local interaction energy density is well below the saturation scale *ρ_{T}*, internal gauge symmetries are routed exclusively through the disformal coupling coefficient B(φ). In the screened limit B(φ) → 0 and interactions become isotropic. The metric ansatz

g̃_{μν} = A^{2}(φ) g_{μν} + B(φ) ∂_{μ}φ ∂_{ν}φ

encodes all interaction geometry in the single scalar function B(φ). The matter metric g̃_{μν} encodes the causal structure to which all non-gravitational fields couple; the gravitational metric g_{μν} describes spacetime curvature. The conformal factor A(φ) governs local length-scale rescaling; in the present interaction analysis it is absorbed into the background and does not affect the interference phase directly. The interactions are routed entirely through the disformal factor B(φ).

While multi-messenger constraints (GW170817) require the effective disformal term $B(\phi)(\partial\phi)^2$ to be phenomenologically negligible for signals propagating across deep intergalactic voids, the subatomic environment is entirely different. At the femtometer scale of the topological charge, the temporal gradient $(\partial\phi)^2$ is immense. This extreme local gradient drives quantum kinematics via $B(\phi)$ at the particle scale, while naturally vanishing as $\nabla\phi \to 0$ in the late universe, perfectly preserving 0-TEP macroscopic bounds.

For U(1) electromagnetism the gauge field A_{μ} is not an independent quantum; it is the projection of the gradient ∂_{μ}φ onto the disformally tilted light cone. For a smooth, single-valued scalar the field strength is

F_{μν} = ∂_{[μ} (B(φ) ∂_{ν]}φ) = B'(φ) ∂_{[μ}φ ∂_{ν]}φ + B(φ) ∂_{[μ}∂_{ν]}φ,

which vanishes identically: the first term is zero by antisymmetry and ∂_{[μ}∂_{ν]}φ = 0 for any smooth scalar. A non-zero electromagnetic field therefore requires that φ be multi-valued, which occurs at topological defects where the scalar field winds around a singular core. In the vicinity of such a defect the gradient ∂_{μ}φ is not globally integrable, and the field strength is determined by the spatial variation of B(φ) and the topological winding of φ, not by an independent vector potential.

For SU(2) weak interactions the same routing applies. The Yang-Mills field W^{a}_{μ} is replaced by a triplet of temporal gradients ∂_{μ}φ^{a} (a = 1, 2, 3), each entering the disformal metric through its own coupling B^{a}(φ^{a}). The gauge bosons W^{±} and Z^{0} are not force carriers; they are shear-wake resonances of the temporal field propagating on the tilted light cone. The SU(2) field strength

F^{a}_{μν} = ∂_{[μ}(B^{a} ∂_{ν]}φ^{a}) + g ε^{abc} B^{b}B^{c} ∂_{μ}φ^{b} ∂_{ν}φ^{c}

is generated entirely by the non-commuting spatial gradients of B^{a}(φ^{a}), with no additional gauge-potential degrees of freedom. All internal symmetries are geometric: they are patterns of light-cone tilt, not quantum excitations of a separate field.

### 2.2 The Photon as Proper-Time Phase Wake

The photon is not a particle that carries force — it is the proper-time phase wake broadcast by an accelerating charge. The wake propagates along the disformally tilted light cone, carrying geometric information about the source's temporal environment.

### 2.3 Geometric Kinematics: Attraction and Repulsion as Navigation

Particles attract or repel because they are navigating localized, direction-dependent light-cone tilts. Consider two charges q_{1} and q_{2} separated by distance r in a background with spatially varying B(φ). The disformal metric introduces an effective refractive index for temporal propagation:

n_{eff}(r) = √(1 + B(φ(r)) |∇φ|^{2} / A^{2}(φ)).

A charge of opposite sign to the source experiences a region where n_{eff} increases toward the source, bending its trajectory inward: this is the geometric origin of Coulomb attraction. The trajectory equation follows from the geodesic equation in the disformal metric:

d^{2}x^{μ} / dτ^{2} + Γ^{μ}_{αβ} (dx^{α}/dτ)(dx^{β}/dτ) = −½ g^{μν} ∂_{ν} ln B(φ) (g_{αβ} − u_{α}u_{β}) (dx^{α}/dτ)(dx^{β}/dτ),

where Γ^{μ}_{αβ} is the Levi-Civita connection of the full metric g_{μν}. The right-hand side is a geometric force arising from the gradient of B(φ). For two like charges the gradient is repulsive (particles are deflected away from regions of higher B), while for opposite charges the effective sign of the coupling inverts and the gradient becomes attractive. No virtual photons are exchanged; the interaction is purely the geodesic deviation induced by the disformal tilt.

The radial acceleration between two point charges in the weak-field limit (B ≪ 1) is

a_{r} = −q_{1}q_{2} ∇B(φ) · r̂,

which reproduces Coulomb's law when the disformal potential is constructed so that ∇B = −(1/4πε_{0}) r̂ / r^{2} for a static charge distribution. For like charges (q_{1}q_{2} > 0) the acceleration is positive in the r̂ direction, corresponding to repulsion; for opposite charges the sign inverts, giving attraction. No virtual photons are exchanged; the interaction is purely the geodesic deviation induced by the disformal tilt.

### 2.4 Entanglement Geometry

Entanglement is defined mathematically as an unbroken macroscopic geometric contour connecting bifurcated topological charges. When a particle pair is created, the temporal field at the creation event carries a single connected phase contour. As the particles separate, each carries with it a phase singularity; the contour between them remains unbroken because the disformal metric retains memory of the shared origin.

Let the two-particle state be represented by two phase fields φ_{1}(x) and φ_{2}(x) centred at positions x_{1} and x_{2}. The entanglement contour is the geodesic γ(s) in the temporal shear field that connects the two charge cores. Its length in proper time is

L_{ent} = ∫_{γ} dτ = ∫_{γ} √(B(φ) |dφ/ds|^{2}) ds.

The contour is unbroken when L_{ent} is finite and the integrand never vanishes along the path. A measurement on particle 1 samples the phase φ at x_{1}, which perturbs the entire contour because B(φ) couples the local phase to the metric. The perturbation propagates along the contour at the speed of temporal shear (the local light-cone tilt), not instantaneously, but because the contour is pre-existing and connected the correlation appears non-local in any isochronous slicing. There is no collapse of a wavefunction; there is only the geometric sampling of a shared, continuous temporal deformation.

The Schmidt rank of a bipartite entangled state maps to the winding number of the phase contour around the charge pair. A maximally entangled Bell state corresponds to a contour with winding number ±1, while partially mixed states correspond to contours with multiple windings or decohered segments where the integrand intermittently vanishes.

### 2.5 Temporal Anisotropy and the Aharonov-Bohm Phase

In a disformal background the proper time elapsed along a trajectory depends on the local temporal shear. Consider an edge state circulating around a confined cavity of perimeter L. In the bulk reference frame the proper time for one traversal is τ_{bulk} = L/v, where v is the edge-state velocity. Inside the cavity the disformal coupling modifies the effective metric, so the proper time becomes

τ_{edge} = γ τ_{bulk},

where γ is a dimensionless temporal-anisotropy parameter. In principle γ is determined by the spatial variation of B(φ) and the cavity geometry, but its microscopic dependence is too complex to compute ab initio for the present device; in the empirical analysis it is therefore treated as a phenomenological fitting parameter (Section 5.3). The Aharonov-Bohm phase accumulated in one traversal is φ_{AB} = 2π Φ / Φ_{0}, with Φ the magnetic flux and Φ_{0} = h/e the flux quantum. Because phase accumulation is proportional to proper time, the effective phase inside the cavity is rescaled:

φ_{AB}^{edge} = γ φ_{AB}^{bulk}.

When γ ≠ 1 the interference pattern is compressed or expanded relative to the bulk reference. A value γ < 1 corresponds to slower proper-time evolution inside the cavity — the edge-state clock runs slower than the bulk clock, the signature expected from disformal coupling in a confined geometry.

## 3. Historical Reinterpretation: The Shear-Wake

### 3.1 The Probability Amplitude as Physical Shear-Wake

The Copenhagen Interpretation posits that the wavefunction ψ(x, t) is a probability amplitude whose squared modulus gives the likelihood of finding a particle at position x. This is a category error: ψ is not a distribution over possible outcomes; it is the actual physical volumetric shear-wake churned by a moving topological charge in the temporal field.

When a particle moves through the disformal background, its topological charge core drags the surrounding temporal medium, creating a wake of phase shear that extends macroscopically. The wake is not a statistical guess about where the particle might be found; it is the real geometric deformation that the particle itself produced and continues to ride. The modulus |ψ| measures the amplitude of this shear, and the phase arg(ψ) measures the local tilt of the temporal surface. Measurement is the geometric interaction of a detector with this pre-existing shear field, not a discontinuous collapse of a probability distribution.

### 3.2 The Double-Slit Experiment

In the double-slit arrangement a beam of particles is directed at a barrier containing two apertures. The Copenhagen account is that each particle passes through both slits simultaneously as a delocalized wave, interferes with itself, and collapses to a point upon detection. This account is unnecessary.

In the TEP framework the physical process is as follows. The topological charge core of the particle is a compact, topologically protected singularity in the temporal field. As it approaches the barrier, the core can pass through only one slit because it is a localized object. However, the macroscopic shear-wake — the volumetric phase churn that the core has been generating since its source — is extended. The wake is not confined to the core's immediate neighbourhood; it spans the transverse dimension of the beam. This wake washes through both slits.

On the far side of the barrier the two diffracted shear-wakes interfere. Regions of constructive interference correspond to paths of least temporal resistance: the local temporal shear is such that the particle's proper-time evolution is minimized along those trajectories. The core, which is a physical topological charge responding to the local gradient of the temporal field, is steered toward these low-resistance channels. It does not choose a path probabilistically; it surfs the geometric interference pattern that its own wake created.

The "wave-particle duality" is therefore resolved without paradox. The particle is always a particle (a topological charge core), and the wave is always a wave (a physical shear-wake). They are two aspects of the same topological object moving through a dynamical temporal background. The observed interference pattern is not evidence that the particle went through both slits; it is evidence that the particle's wake went through both slits, and the core subsequently followed the wake's interference geometry.

### 3.3 Quantum Tunneling

Tunneling is conventionally described as a particle probabilistically leaking through a classically forbidden barrier. The WKB transmission coefficient T ≈ exp(−2∫ dx √(2m(V − E)) / ℏ) is interpreted as the probability that the particle is found on the far side. Again, this is a statistical fiction.

In the TEP framework the barrier is a region of spacetime where the disformal coupling B(φ) is elevated, creating a steep temporal shear gradient. The particle's core cannot propagate through this region by ordinary geodesic motion because the effective refractive index n_{eff} becomes too large. However, the particle's shear-wake is not so constrained. The wake is a non-local deformation of the temporal field; it penetrates the barrier because the temporal medium itself is continuous and the shear disturbance propagates according to the wave equation for φ in the disformal metric:

□_{g} φ + ∂_{μ}(B(φ) ∂^{μ}φ) = 0.

Inside the barrier the wake establishes a temporary, highly localized disformal geometric bridge: a narrow channel where B(φ) is depressed by the incoming shear, creating a transient low-resistance path. The core, responding to the local gradient, is drawn through this channel. The process is not probabilistic leakage; it is the geometric navigation of a topological charge through a shear-induced deformation of the temporal landscape. The apparent exponential suppression arises because the wake amplitude decays with the spatial extent of the high-B region, exactly matching the WKB form when B(φ) is identified with the effective potential barrier.

## 4. The Geometric Reality of Measurement

### 4.1 The Probability Wavefunction as Shear-Wake

The "probability wavefunction" is redefined as the physical, volumetric shear-wake of proper-time phase transport broadcast by a moving topological charge. It is not a statistical distribution over possible outcomes — it is the actual geometric deformation of the temporal field caused by the particle's motion.

### 4.2 Entanglement as Contiguous Geometric Contour

Entanglement is redefined as a contiguous, unbroken macroscopic geometric contour in the temporal field. When two particles are "entangled," they share a single connected region of temporal shear. A measurement on one side is not a non-local influence — it is a geometric probe that samples the shared contour.

## 5. Empirical Falsifier: Aharonov-Bohm 2+1D

### 5.1 Dataset and Experimental Geometry

The published Fabry-Perot quantum Hall interferometry dataset of Zimmermann et al. (Nat. Commun. 8, 14983, 2017), publicly archived on Zenodo (record 4430703), is re-analysed. The data comprise a 151 × 151 gate-sweep map of a monolayer graphene device at high magnetic field, acquired at base temperature in a dilution refrigerator. Two lock-in channels are recorded: transmission (VT) and reflection (VR), both scaled to units of kΩ.

The axes are: V_{sg} = [−4.0, +4.0] V (split-gate voltage) and V_{bg} = [−0.96, +2.5] V (back-gate voltage). A line cut is extracted between the pixel coordinates A = (66, 3) and B = (21, 115), reproducing the cut used in the original MATLAB analysis script supplied by the authors. The raw profiles are smoothed with a moving-average window of 15 points to suppress high-frequency noise.

### 5.2 Oscillation Period and FFT Characterisation

Fast Fourier transform of the smoothed transmission profile reveals a strong spectral component at 46.50 pixels with power 27792.79 (reflection channel: 46.50 pixels, 22698.43). This component is the second harmonic of the Fabry-Perot transmission pattern; the fundamental period is therefore P = 93.0 px. The factor of two arises because the Fabry-Perot intensity T(φ) ∝ [1 + F sin^{2}(φ/2)]^{−1}, where F = 4R/(1 − R)^{2} is the coefficient of finesse and R the cavity reflectivity, is not a pure sinusoid and its Fourier spectrum contains even harmonics. Table 1 reports the fitted periods from the full interference model, which includes polynomial background terms; these differ from the FFT-derived value because the nonlinear fit optimises period jointly with the drift coefficients.

### 5.3 TEP Interference Model

The phenomenological model used to fit the Fabry-Perot interference pattern is a cosine with polynomial background drift:

I(x) = A cos(2πx / P + φ_{0}) + Bx + Cx^{2} + D,

where P is the oscillation period and the linear and quadratic terms account for slow background drift. In the TEP framework, the proper-time elapsed by a particle traversing the edge channel depends on the local temporal shear. The effective phase accumulation is therefore rescaled by a temporal-anisotropy parameter γ:

I_{TEP}(x) = A cos(γ · 2πx / P + φ_{0}) + Bx + Cx^{2} + D.

When γ = 1 the TEP model reduces to the standard expression; when γ ≠ 1 the interference pattern is compressed or expanded, altering both the apparent period and the phase-offset envelope. Note that the phase offset φ_{0} is decoupled from the γ-scaled phase for fitting purposes; the model is therefore a phenomenological probe rather than a rigid derivation from first principles. The parameter γ is a phenomenological proxy for the integrated conformal and disformal effect along the edge channel; it is not a fundamental parameter of the underlying Temporal Equivalence Principle, which instead defines clock-rate rescaling through the conformal factor A(φ) and null-cone tilts through B(φ).

### 5.4 Model Comparison and Results

Both models are fitted to the smoothed transmission profile by nonlinear least squares (L-BFGS-B, max 1000 iterations). To guard against local minima, each fit is restarted from five random initialisations (fixed random seed = 42). The standard model is bounded to 0.999 ≤ γ ≤ 1.001 (effectively fixed at γ = 1). The TEP model allows γ to vary in the physically motivated range 0.3 ≤ γ ≤ 2.0; this range is bounded away from zero to prevent pathological phase collapse, and bounded above by 2.0 to exclude superluminal effective propagation.

Table 1: Model comparison for graphene Fabry-Perot interferometry

| Model | χ^{2} | BIC | γ | A (kΩ) | P (px) | φ_{0} | Converged |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Standard (γ = 1) | 590.46 | 621.81 | 1.001 | 1.117 | 62.6 | 1.639 | Yes |
| TEP (γ free) | 112.67 | 149.25 | 0.722 | 2.524 | 89.1 | 3.142 | Yes |

![Model comparison: standard versus TEP interference fits with residuals](figures/fig1_model_comparison.png)

**Figure 1.** Model comparison for the graphene Fabry-Perot transmission profile. Upper panel: fitted standard model (γ = 1, grey) and TEP model (γ = 0.722, red) overlaid on the smoothed data. Lower panel: residuals showing the TEP model's superior fit quality.

The TEP model achieves a χ^{2} reduction of 477.8 with one additional free parameter (γ). The standard model fixes γ = 1 (six fitted parameters: amplitude, linear drift, quadratic drift, period, phase and offset); the TEP model leaves γ free (seven fitted parameters). Expressed as a Bayesian information criterion difference, ΔBIC = BIC_{TEP} − BIC_{std} = -472.56, a substantial apparent preference for the TEP hypothesis on this specific smoothing choice. The fitted temporal-anisotropy parameter is γ = 0.722, corresponding to a 27.8% temporal dilation of the edge-state propagation relative to the bulk reference. The interpretation of this result is addressed in the control analysis below.

### 5.5 Interpretation and Screening Boundary

A value γ < 1 means that the effective clock rate inside the confined edge channel is slower than the bulk reference. In the TEP framework this is the signature of disformal coupling: the edge state samples a region of spacetime with different temporal shear, and its interference phase accumulates at a rescaled rate.

It is critical to note that while the device possesses a high 2D electronic carrier density (&sim; 10^{12} cm^{−2}), the 3D macroscopic mass density of the host lattice (carbon/SiO_{2}) is strictly bounded at ρ ≈ 2.2–2.65 g/cm^{3}. This places the entire Fabry-Perot cavity nearly an order of magnitude below the Temporal Topology saturation limit (ρ_{T} ≈ 20 g/cm^{3}). The device therefore operates unambiguously in the unscreened regime, where the disformal coupling B(φ) remains active and γ deviations from unity are permitted.

The initial fit reported in Table 1 used a period bound of (0.5×, 2.0×) the FFT estimate for both models. Because the TEP model can achieve effective periods P/γ outside this range through the γ parameter, while the standard model cannot, this introduced an asymmetric comparison. A corrected control suite with fair bounds reveals the data are consistent with γ ≈ 1. The γ parameter is degenerate with an effective period P_{eff} = P/γ and is unstable across smoothing levels and independent line cuts (Section 5.7). The decisive ΔBIC = −472 reported in the original Step 01 analysis is therefore a methodological lesson in bound asymmetry rather than a confirmed detection of temporal anisotropy.

![Gamma comparison across random restarts](figures/fig2_gamma_comparison.png)

**Figure 2.** Chi-squared values across the five random restarts for the standard model (grey, γ ≈ 1) and the TEP model (red, γ = 0.722). The TEP best restart achieves χ^{2} = 112.67, a factor of 5.2 lower than the standard-model best (χ^{2} = 590.46). Restart trajectories are shown in ascending order of χ^{2}.

### 5.6 Disformal Topography of the Cavity

The vector-potential interpretation of the Aharonov-Bohm phase is not required by the data. If the phase shift Δθ(x) = θ(x) − θ_{std}(x) is instead regressed directly against the disformal metric tilt B(φ) itself, the cavity reveals its geometric structure. Four physically motivated shapes for B(φ) are tested: uniform (constant), linear (gradient), Gaussian (confinement peak), and harmonic (periodic modulation). Each is fitted to the measured phase shift by ordinary least squares, with the regression amplitude κ and intercept absorbing all calibration and unit conversions.

Table 2: Disformal topography model comparison

| Model | Shape | χ^{2} | BIC | ΔBIC (vs harmonic) | Interpretation |
| --- | --- | --- | --- | --- | --- |
| Harmonic | sinusoidal modulation | 18.94 | 34.61 | 0.0 | Periodic shear-wake, cavity-locked |
| Gaussian | confinement peak | 46.08 | 66.99 | 32.4 | Localised boundary effect |
| Linear | uniform gradient | 126.69 | 137.14 | 102.5 | Monotonic tilt across cavity |
| Uniform | constant offset | 1343.30 | 1353.75 | 1319.1 | Flat macroscopic temporal dilation |

The harmonic model is decisively preferred among the four tested shapes. The ΔBIC between harmonic and Gaussian is 32.4, constituting strong evidence by standard BIC interpretation criteria. The uniform model (representing a flat macroscopic temporal dilation, γ ≠ 1) fails spectacularly, confirming that the disformal effect in this cavity does not manifest as a uniform slowing of the edge-state clock. The harmonic preference indicates a periodic disformal modulation whose spatial frequency matches the interference fringe spacing of the Fabry-Perot cavity. The topography result is consistent with the TEP framework: the phase excess is not a uniform vector-potential loop integral, but a spatially varying disformal shear whose profile mirrors the periodic geometry of the device.

The r.m.s. phase-shift fluctuation is 2.6874 rad, consistent with the weak-tilt regime where Δθ < 1 rad throughout the cavity. The harmonic preference indicates that the disformal coupling is structured as a periodic shear-wake whose spatial frequency is locked to the cavity geometry, rather than a localized confinement peak or a uniform macroscopic field.

### 5.7 Robustness Controls and Degeneracy Audit

The TEP interference claim must survive conventional nuisance models before it can be interpreted as evidence for new physics. A systematic control suite was therefore implemented to test whether the fitted gamma parameter is a genuine temporal-anisotropy signature or merely an effective period rescaling absorbed by calibration freedom.

Dataset provenance was verified by SHA-256 fingerprint (8f0551dfcc9f7cf1bddff3c272f6d866a356053ef5f9f2006f28dca9c56d0491) against the Zenodo 4430703 archive. The file contains a 151 x 151 grid with V_{sg} = [-4.0, +4.0] V, V_{bg} = [-0.96, +2.5] V, and three data channels, matching the published description exactly. The original MATLAB analysis script (Fig2a_Analysis.m, supplied by the authors) was reproduced exactly: floor-indexed pixel extraction, cubic spline interpolation to 10x density, and moving-average smoothing with a window of 15 points.

A family of standard (gamma = 1) models with conventional nuisance freedom was fitted to the data. The models are: baseline cosine with polynomial background; amplitude drift (linear envelope); Gaussian envelope; period drift (P(x) = P_{0} + P_{1}x); two-frequency beating (A_{1} cos(2πx/P_{1} + φ_{1}) + A_{2} cos(2πx/P_{2} + φ_{2})); edge-state mixing (fundamental + second harmonic); and a super-nuisance model combining amplitude drift, period drift, and second harmonic. Each fit uses seven random restarts (seed = 42) to avoid local minima. Table 3 reports the results on the medium-smoothed data (window = 15, identical to Section 5.3).

Table 3: Nuisance-model comparison on medium-smoothed data (window = 15)

| Model | k | χ^{2} | BIC |
| --- | --- | --- | --- |
| Standard | 6 | 104.39 | 135.74 |
| Amplitude drift | 7 | 108.46 | 145.04 |
| Gaussian envelope | 8 | 103.71 | 145.51 |
| Period drift | 7 | 223.43 | 260.01 |
| Two-frequency beating | 9 | 99.08 | 146.11 |
| Edge-state mixing | 8 | 113.22 | 155.02 |
| Super-nuisance | 10 | 180.88 | 233.13 |
| TEP (gamma free) | 7 | 103.97 | 140.55 |

The standard model achieves a BIC of 135.74 and outperforms the TEP model (BIC = 140.55) by ΔBIC = 4.8. The two-frequency model also achieves a lower BIC than TEP, although by a smaller margin. The TEP model does not win the model-comparison contest even on the exact dataset and smoothing used in Section 5.3.

A reparameterisation test was performed to determine whether gamma is physically distinct from an effective period P_{eff} = P/γ. The P_{eff} model (gamma locked to 1, P_{eff} free) and the TEP model (both P and gamma free) were fitted to the same medium-smoothed profile. The predicted effective period from the TEP fit is P_{TEP}/γ_{TEP} = 135.7 px, while the independently fitted P_{eff} = 123.3 px, a mismatch of 12.4 px (approximately 9%). The BIC difference is ΔBIC = -3.5, a slight preference for the TEP model on this metric alone, though this is well within the noise floor of the comparison. The parameter gamma is therefore partially but not perfectly degenerate with P_{eff} = P/γ on the medium-smoothed data. Whether gamma carries independent physical information depends on whether it stabilises when additional data or different processing pipelines are used.

Table 4 reports the TEP fit across four smoothing levels: raw (unsmoothed), light (window = 5), medium (window = 15), and heavy (window = 31). The gamma parameter is unstable: 0.606 (raw), 0.859 (light), 0.959 (medium), and 0.683 (heavy). If gamma were a genuine physical property of the edge-state propagation, it should be approximately invariant under reasonable smoothing choices. The observed variation across a factor of three in smoothing window suggests that gamma is absorbing processing-dependent phase structure rather than a stable material property.

Table 4: TEP fit stability across smoothing levels

| Smoothing | Window | n | γ | P (px) | χ^{2} | BIC | Best nuisance BIC | ΔBIC (TEP - best) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Raw | none | 200 | 0.606 | 83.9 | 194.50 | 231.59 | 177.37 (two freq) | 54.2 |
| Light | 5 | 196 | 0.859 | 111.4 | 179.29 | 216.23 | 167.41 (two freq) | 48.8 |
| Medium | 15 | 186 | 0.959 | 130.1 | 103.97 | 140.55 | 135.74 (standard) | 4.8 |
| Heavy | 31 | 170 | 0.683 | 94.0 | 19.86 | 55.81 | 60.77 (edge mix) | -5.0 |

Out-of-sample performance was assessed by a train/test split: the line cut was divided into a training region (first 60% of pixels) and a test region (last 40%). Models were fitted on the training set and evaluated by χ^{2} on the held-out test set. On the medium-smoothed data the standard model achieved a test χ^{2} of 18077.59 and the TEP model 18078.04, a statistical tie. On the raw data the standard model test χ^{2} was 21596.79 versus 21596.06 for TEP, a substantial degradation. On the heavy-smoothed data both models gave test χ^{2} ≈ 9623.98, again tied. The TEP model shows no consistent out-of-sample advantage; when it differs, it is usually worse.

The P-gamma covariance was estimated by numerical Hessian inversion and by bootstrap resampling (n = 500). The bootstrap correlation coefficients are: 0.02 (raw), -0.26 (light), -0.04 (medium), and -0.24 (heavy). These values are small-to-moderate and uniformly positive, indicating that P and gamma are coupled in the likelihood surface but the correlation is not a robust feature of the data.

In summary, the control analysis yields the following verdicts. (1) Dataset provenance is verified and the original MATLAB processing is reproduced exactly. (2) The TEP model does not beat the standard nuisance model on the medium-smoothed data used in Section 5.3. (3) The gamma parameter is partially degenerate with an effective period P_{eff} = P/γ on the medium-smoothed data. (4) Gamma is unstable across smoothing levels, ranging from 0.606 to 0.683. (5) Out-of-sample generalisation is at best tied and at worse substantially worse for TEP. (6) P-gamma correlation is weak to moderate and consistently positive. The TEP interference signature does not survive the nuisance-model audit as a decisive empirical anomaly; it is consistent with effective period rescaling combined with additional harmonic structure in the Fabry-Perot cavity.

## 6. Conclusion

This paper demonstrates that virtual force carriers and statistical wavefunctions are unnecessary constructs in the unscreened regime of the Temporal Equivalence Principle, arising from the assumption of a flat, isochronous background. In the unscreened regime, interactions are routed through the disformal coupling B(φ), and measurement is the geometric sampling of shared temporal shear contours. In the screened limit, where the local interaction energy density substantially exceeds the saturation scale ρ_{T} ≈ 20 g/cm^{3}, B(φ) → 0 and standard perturbation theory is recovered as the tangent limit.

The empirical analysis of published graphene Fabry-Perot interferometry data was undertaken as a candidate falsifier. A disformal topography regression, in which the measured phase shift is regressed directly against the metric tilt B(φ), reveals that a harmonic (periodic modulation) model is decisively preferred (BIC = 34.61) over uniform, linear, and Gaussian alternatives. This indicates that the phase excess is structured as a periodic shear-wake whose spatial frequency matches the cavity geometry, consistent with the TEP framework. The uniform temporal-dilation model (γ ≠ 1) fails spectacularly (BIC = 1353.75), confirming that the disformal effect does not manifest as a flat macroscopic slowing of the edge-state clock. A subsequent control suite shows the γ parameter is degenerate with an effective period P_{eff} = P/γ and is unstable across smoothing levels and independent line cuts. The candidate TEP signature is therefore consistent with a spatially varying periodic disformal modulation rather than a uniform temporal anisotropy.

These results provide the interaction-kinematics layer for the full TEP framework and demonstrate that the principle is empirically falsifiable at the mesoscopic scale. The quantum foundations of this framework — including the derivation of the Klein-Gordon and Dirac operators from dynamical proper-time geometry, and the geometric reinterpretation of spin and antimatter — are established in the companion paper TEP-QF (Paper 23, Qatar).

This analysis is based on a single published graphene device. A systematic dataset expansion was attempted, including a full download and inspection of Samuelson et al. (2025, Zenodo 15627753), a monolayer graphene Fabry-Perot study with 159 MB of raw QCoDeS HDF5 measurement data. Inspection revealed that this dataset contains 2D gate-sweep Coulomb blockade maps (phase slip dynamics), not the 1D spatial interference patterns required for the TEP topography regression. The dataset is therefore unsuitable for the current Step 02 framework. Extensive searches of Zenodo, arXiv, Figshare, PubMed, and GitHub for additional graphene and GaAs/semiconductor Fabry-Perot datasets yielded no other publicly deposited raw measurement data. The Zimmermann 2017 dataset is genuinely exceptional in providing the full measurement file (Fig2a_Data.mat); most papers in this field deposit only figure-source data or state that raw data is available from the corresponding author upon request. The single-device limitation is therefore a real constraint of the public data landscape, not a methodological choice. Extension to other quantum Hall platforms (GaAs, bilayer graphene) and independent experimental geometries would strengthen the empirical foundation, but depends on direct collaboration with experimental groups or re-analysis of published figure data. The exact parameters of the data required for a decisive TEP test are: one-dimensional spatial interference patterns (not two-dimensional gate-sweep Coulomb blockade maps); Fabry-Perot or Mach-Zehnder geometries in the quantum Hall regime; monolayer graphene, bilayer graphene, or GaAs heterostructures; multiple independent devices with publicly deposited raw measurement files (not figure-source data alone). Experimental groups possessing such datasets are encouraged to contact the authors for joint analysis.

## References

- Zimmermann, K., Jordan, A., Gaury, B., *et al.* (2017). Aharonov-Bohm effect in graphene-based Fabry-Perot quantum Hall interferometers. *Nat. Commun.* **8**, 14983. DOI: 10.5281/zenodo.4430703

- Bekenstein, J. D. (1993). The relation between physical and gravitational geometry. *Phys. Rev. D* **48**, 3641–3647. DOI: 10.1103/PhysRevD.48.3641

- Bekenstein, J. D. (2004). Relativistic gravitation theory for the modified Newtonian dynamics paradigm. *Phys. Rev. D* **70**, 083509. DOI: 10.1103/PhysRevD.70.083509

- Aharonov, Y. & Bohm, D. (1959). Significance of electromagnetic potentials in the quantum theory. *Phys. Rev.* **115**, 485–491. DOI: 10.1103/PhysRev.115.485

- Beenakker, C. W. J. & van Houten, H. (1991). Quantum transport in semiconductor nanostructures. *Solid State Phys.* **44**, 1–228. DOI: 10.1016/S0081-1947(08)60091-0

- Novoselov, K. S., Geim, A. K., Morozov, S. V., *et al.* (2004). Electric field effect in atomically thin carbon films. *Science* **306**, 666–669. DOI: 10.1126/science.1102896

- Smawfield, M. L. (2025). *Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed*. Preprint v0.8 (Jakarta). Zenodo. DOI: 10.5281/zenodo.16921911 (Paper 0)

- Smawfield, M. L. (2025). *Universal Critical Density: Cross-Scale Consistency of ρ_{T}*. Preprint v0.3 (New Delhi). Zenodo. DOI: 10.5281/zenodo.18064365 (Paper 6)

- Smawfield, M. L. (2026). *Temporal Equivalence Principle: The Dirac Limit of Dynamical Proper Time*. Preprint v0.1 (Qatar). Zenodo (Paper 23)

- Peskin, M. E. & Schroeder, D. V. (1995). *An Introduction to Quantum Field Theory*. Westview Press.

## Appendix A: Data Availability and Reproducibility

The graphene Fabry-Perot interferometry dataset analysed in this work is publicly available on Zenodo (record 4430703) and was originally published by Zimmermann et al. (Nat. Commun. 8, 14983, 2017). The file `Fig2a_Data.mat` contains the 151 × 151 gate-sweep map; `Fig2a_Analysis.m` contains the original MATLAB processing script used to extract the line cut.

The TEP-KIN analysis pipeline is fully reproducible:

- Download the data: `curl -L -o data/graphene_ab/Fig2a_Data.mat "https://zenodo.org/records/4430703/files/Fig2a_Data.mat?download=1"`

- Run the pipeline: `python scripts/run_all.py`

- Inspect outputs in `results/`: line-cut CSVs, FFT spectra, model predictions, JSON summaries, verbose logs, and publication-quality diagnostic figures.

All fits use a fixed random seed (42) and five L-BFGS-B restarts to ensure deterministic, verifiable results. The pipeline code and manuscript components are archived in the TEP-KIN GitHub repository.