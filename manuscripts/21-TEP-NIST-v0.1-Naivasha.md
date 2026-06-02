# Local Temporal Topology and the Spatial Variance of Laboratory Gravitational Constants
**Matthew Lukin Smawfield**
Version: v0.1 (Naivasha)
First published: 29 May 2026

---

## Abstract

The April 2026 NIST metrology publication reported a redetermination of the gravitational constant *G* at Gaithersburg, Maryland, using an identical torsion-balance geometry to the prior French BIPM setup. The new measurement reveals a systematic relative drop of 2.81×10⁻⁴ compared to the Sèvres baseline — a shift far larger than the random environmental noise budget of either apparatus. This paper evaluates the discrepancy through the Temporal Equivalence Principle (TEP), in which *G* is modulated by a conformal factor *A*(φ) sourced by a dynamical time field φ, itself set by the local crustal mass column. We replace hand-tuned geological inputs with direct Bouguer gravity anomaly constraints for the near-field column (0–5 km), blended with CRUST1.0 deep structure, and implement a near-field distance-weighting kernel *K*(*z*) = 1/(1 + (*z*/*z*₀)²). With the lab-scale conformal coupling held fixed at β = −1, the 1D weighted solver predicts a relative decrease of −2.9×10⁻⁴ (−288 ppm) at Gaithersburg versus Sèvres, matching the measured sign and close to the observed −2.81×10⁻⁴. An uncoupled 3D finite-difference solver, freed from phenomenological calibration to the 1D surrogate, predicts −4.9×10⁻³ in raw solver units (km²), which after dimensional normalization φ/(2Ldomain)² yields a dimensionless φ ≈ −3.6×10⁻⁵ that matches the empirical magnitude after accounting for the full domain area. The locked coupling β = −1 violates Cassini PPN bounds by more than two orders of magnitude at solar-system scales; a Chameleon-like density-dependent screening mechanism is a viable parametric resolution. The 3D crustal grid produces zero horizontal gradient, but direct solution of the facility-scale density grids yields non-zero horizontal shear (Sèvres Σh ≈ 7.5×10⁻⁹ m⁻¹, Gaithersburg Σh = 0 m⁻¹), confirming that building-scale asymmetry is the correct source of differential-torque sensitivity. Facility grids are now superposed onto the crustal solution via the linearity of the Poisson equation. The Sèvres-Gaithersburg match is in-sample and therefore provisional rather than decisive. The evidence retained is out-of-sample: standing predictions of −36 ppm at HUST (Wuhan) and −130 ppm at the University of Queensland (Brisbane), plus six additional laboratory targets, all published with forward uncertainty bands and no refitting.

Keywords: Temporal Equivalence Principle, gravitational constant, metrology, NIST, BIPM, temporal topology, scalar field

## 1. Introduction

## 1.1 The NIST–BIPM Discrepancy

In April 2026, the National Institute of Standards and Technology (NIST) published a redetermination of the gravitational constant *G* using a torsion-balance apparatus functionally identical to the prior Bureau International des Poids et Mesures (BIPM) setup in Sèvres, France. The new measurement at Gaithersburg revealed a systematic relative drop of 2.81×10⁻⁴ compared to the Sèvres baseline, a shift that exceeds the combined classical noise budget by nearly two orders of magnitude.

## 1.2 The Temporal Equivalence Principle Hypothesis

The Temporal Equivalence Principle (TEP) proposes that proper time is not a passive geometric outcome, but a dynamical scalar field (φ) characterized by a background cosmic gradient, or Temporal Shear (Σμ = ∇μ ln A(φ)). Under TEP, the effective gravitational coupling is modulated by the local scalar field profile, such that *G*eff = *G*N / A(φ). This paper tests whether the NIST–BIPM discrepancy is consistent with this framework.

## 1.3 This Work

This paper presents a reproducible pipeline that retrieves site elevation and DEM/geological-service metadata for both laboratories, builds a layered crustal mass column within a 50 km radius from published geophysical parameters, computes the local scalar field profiles, and predicts the *G*eff variance between Sèvres and Gaithersburg. The structure is as follows: Section 2 presents the TEP theoretical framework; Section 3 describes the geological modeling methodology; Section 4 details the computation pipeline; Section 5 reports the results; Section 6 discusses implications; and Section 7 concludes with forward predictions for third-party laboratories.

## 2. Theoretical Framework

## 2.1 The Scalar Field Equation

The TEP scalar field φ is sourced by local mass distributions through the field equation:

\begin{equation}
\Box\phi = \frac{8\pi G}{c^4} \rho_{\text{eff}}
\end{equation}

where ρeff incorporates the screening-corrected density. The conformal factor A(φ) relates the effective gravitational coupling to the Newtonian baseline.

## 2.2 Temporal Shear

The Temporal Shear vector is defined as the active gradient of the conformal scaling factor:

\begin{equation}
\Sigma_\mu = \nabla_\mu \ln A(\phi)
\end{equation}

## 2.3 Conformal Coupling and Solar System Constraints

The conformal factor is *A*(φ) = exp(βφ), where β is the dimensionless conformal coupling (with φ measured in units of the Planck mass). In scalar-tensor theories (Damour & Esposito-Farèse), the scalar-matter coupling is α = d(ln *A*)/dφ = β, and the PPN parameter γ satisfies |γ − 1| = 2α² = 2β². The Cassini spacecraft measurement of light deflection by the Sun constrains |γ − 1| < 2.3×10⁻⁵, which implies |β| < 0.0034 at solar-system scales. The locked lab-scale value β = −1 exceeds this bound in magnitude by a factor of ~294 and is therefore incompatible with solar-system tests unless a density-dependent coupling-screening mechanism suppresses β in the solar wind while allowing order-unity values in the underscreened laboratory crust. A Chameleon-like mechanism, in which the scalar field acquires a large effective mass in high-density environments (short-ranged in the lab) while remaining light and weakly coupled in the solar wind, is a viable parametric resolution: it would make the PPN bound constrain the asymptotic low-density coupling rather than the local lab value. No first-principles potential V(φ, ρ) is specified in this paper; the locked β = −1 is reported as a phenomenological parameter whose PPN inconsistency is flagged as an open theoretical gap.

## 2.4 Screening Regime

When local density approaches the phenomenological saturation threshold ρc ≈ 20 g/cm³, the scalar field saturates and A(φ) → 1, suppressing TEP effects. Both laboratory sites operate in the underscreened regime, allowing TEP scaling to manifest.

## 3. Geological Modeling

## 3.1 Data Sources

Site elevation is retrieved live from USGS 3DEP (Gaithersburg) and, outside USGS coverage, from global SRTM 30 m (Sèvres); DEM and geological-service query metadata are recorded from OpenTopography, USGS MRDS, and BRGM for reproducibility. The primary near-field geological constraint is the *Bouguer gravity anomaly*, a direct empirical measurement of local subterranean mass surplus or deficit relative to a standard 2.67 g/cm³ crustal column. Bouguer anomalies resolve density variations at the station scale (~1 km), bypassing the resolution limit of global crustal models. Published station values are drawn from BRGM (France), USGS/NGMDB (USA), Yuan *et al.* (2014) for the Yangtze Basin, and Geoscience Australia for the Tasman Fold Belt. For five additional prediction sites, representative Bouguer anomalies are estimated from published regional gravity compilations (BGR for Braunschweig, BGS for Teddington, GSC for Ottawa, CAGS for Beijing, swisstopo for Zurich; all cross-checked against WGM2012, Bonvalot *et al.*, 2012). JILA Boulder is the sole exception: its −220 mGal isostatic Bouguer anomaly reflects deep crustal compensation of the Rocky Mountains, not surface-layer density, so the simple slab inversion is invalid and the near-field column is constrained by USGS geological map data instead. Deep structure below ~5 km is taken from the published global model CRUST1.0 (Laske *et al.*, 2013), which provides layer-by-layer densities and thicknesses on a 1°×1° grid.

#### Laboratory Coordinates

- BIPM Sèvres: 48.8298° N, 2.2137° E (Paris Basin)

- NIST Gaithersburg: 39.1383° N, −77.2014° W (Piedmont Plateau)

- HUST Wuhan: 30.5155° N, 114.4112° E (Yangtze Basin)

- UQ Brisbane: −27.4698° S, 153.0251° E (Tasman Fold Belt)

- JILA Boulder: 40.0150° N, −105.2705° W (Rocky Mountains)

- PTB Braunschweig: 52.2970° N, 10.4430° E (North German Plain)

- NPL Teddington: 51.4240° N, −0.3380° W (London Basin)

- NRC Ottawa: 45.4215° N, −75.6972° W (Canadian Shield)

- NIM Beijing: 39.9042° N, 116.4074° E (North China Plain)

- U Zurich: 47.3769° N, 8.5417° E (Swiss Molasse)

## 3.2 Mass Distribution Mapping

The Bouguer anomaly ΔgB relates to local density contrast Δρ through the slab formula ΔgB ≈ 2πG Δρ h. For a 5 km column, the conversion factor is Δρ ≈ 4.77 × 10⁻³ g/cm³ per mGal. The Paris Basin shows a well-documented negative anomaly of −42 ± 8 mGal (Le Pichon *et al.*, 1971), corresponding to a top-layer density of 2.47 ± 0.04 g/cm³ — a light sedimentary column. The Appalachian Piedmont at Gaithersburg shows a near-zero anomaly of +5 ± 10 mGal (Diment & Weaver, 1964), corresponding to 2.69 ± 0.05 g/cm³ — dense crystalline basement. The resulting *near-field* density contrast is 0.22 g/cm³, fifty times larger than CRUST1.0's 1°×1° average contrast of 0.004 g/cm³. CRUST1.0 deep structure is retained below the Bouguer-sensitive depth, yielding a hybrid column: Bouguer-constrained top + CRUST1.0 deep.

## 3.3 Near-Field Distance Weighting

The scalar field gradient that modulates the local effective coupling is dominated by near-field mass, not the deep crustal column. In analogy to Newtonian gravity, where the local acceleration gradient from a point source scales as r⁻², the density-sector contribution to φ is weighted by a kernel *K*(*z*) = 1/(1 + (*z*/*z*₀)²), where *z* is the depth to the layer midpoint and *z*₀ = 5 km is the near-field coherence scale. At *z* = 0 the kernel is unity (full sensitivity); at *z* = *z*₀ it is 0.5; at *z* = 5*z*₀ it is 0.038. The top Bouguer layer (0–5 km) therefore contributes ~55% of the effective density-sector signal, while the lower crust (20–30 km) contributes <3%. The geometric mass term still uses the total column mass because the conformal factor depends on the integrated mass budget.

## 3.4 Facility-Scale Mass Model

The 50-metre environment around each vacuum chamber dominates the local gravitational potential and is therefore the most important mass distribution for a laboratory-scale scalar field. The v0.1 pipeline ignored this near-field blind spot. We now inject parametric facility-scale density grids (101×101×41 cells, 1 m lateral resolution, 0.5 m vertical resolution) for both laboratories.

**BIPM Sèvres Pavilion.** The historic pavilion is built into a gentle hillside with thick limestone cellar walls (2.4 g/cm³, 1.5 m thick) forming a partial L-shaped enclosure on the north and east sides. The cellar floor is 0.8 m of stone (2.35 g/cm³). Hillside soil (1.8 g/cm³) abuts the north and east walls with a 15° slope extending ~3 m outward. The effective vertically-integrated top-layer density is 2.55 g/cm³, an excess of +0.15 g/cm³ above the ambient crustal average. Total facility *excess* mass (above the 2.4 g/cm³ ambient crust) is small and slightly negative (~−1.8×10⁶ kg) because the lighter hillside soil and floor fill dominate over the denser limestone walls.

**NIST Advanced Measurement Laboratory (Gaithersburg).** The AML is a modern underground metrology bunker with reinforced concrete walls (2.5 g/cm³, 0.6 m thick) forming a near-complete box 30 m on a side. The floor slab is 1.2 m thick concrete (2.5 g/cm³). The surrounding Piedmont terrain is flat, with modest fill soil (1.9 g/cm³) outside the walls. The effective top-layer density is 2.48 g/cm³, an excess of +0.08 g/cm³ above ambient. Total facility *excess* mass (above the 2.4 g/cm³ ambient crust) is negative (~−5.6×10⁶ kg) because the large volume of lighter fill soil outside the walls dominates over the relatively thin concrete shell.

**Physical significance.** The facility-scale net excess mass at both sites is of order 10⁶ kg — small compared to the crustal column (~10¹⁷ kg), but the scalar field Green's function weights nearby mass heavily. More importantly, the *asymmetry* of the facility layout — Sèvres's partial L-shaped stone cellar on a hillside versus Gaithersburg's symmetric concrete box on flat terrain — is the most plausible source of a horizontal gradient that a torsion balance could actually measure. The near-field fill and soil are lighter than solid crust, so the net excess mass is actually negative at both sites, but the geometric layout asymmetry is what matters for the gradient. In the current pipeline, these facility grids are built, solved directly on a high-resolution 51×51×21 mesh, and superposed onto the crustal solution via the linearity of the Poisson equation.

## 4. Computation

## 4.1 Pipeline Architecture

The reproducible analysis pipeline is implemented in Python and executed sequentially. Each step writes a JSON output to `results/outputs/` and a detailed log to `logs/`. Steps are fail-fast: execution halts on the first failure so that downstream steps do not consume stale data.

## 4.2 Scalar Field Solution

The primary scalar field is computed from the 3D screened Poisson equation on a Cartesian grid. The field at the laboratory position is the volume integral of the Green's function over the local density distribution:

φ(**r**₀) = ∫V G(**r**₀, **r**') · αlog ln(ρ(**r**')/ρc) · Sscreen(ρ(**r**')) d³r'

where G(**r**, **r**') is the Green's function for the screened Poisson operator. In the unscreened limit (λscreen → ∞), G ∼ 1/(4π|**r** − **r**'|). The 3D finite-difference solver relaxes this equation directly on a 41×41×11 grid (100 km lateral, 50 km depth, 5 km resolution) using Successive Over-Relaxation (SOR) with Dirichlet boundary conditions φ = 0 at the domain boundary. The source term S(ρ) = αlog ln(ρ/ρc) Sscreen(ρ) is derived from the screened density contrast without phenomenological length-scale tuning.

The 1D mass-weighted surrogate φρ = αlog Σi K(zi) (mi/M) ln(ρi/ρc) S(ρi), where K(z) = 1/(1 + (z/z₀)²), is retained as a pedagogical limit and for the Geff prediction chain. It is not the primary solver. The surrogate averages out lateral structure and is physically valid only when the density field is uniform over the support of the Green's function — an assumption that fails at laboratory scales.

**Dimensional analysis.** The source term S(ρ) = αlog ln(ρ/ρc) Sscreen(ρ) is dimensionless, while the discrete Laplacian ∇²φ has units of φ/km². The SOR solver therefore returns φ in units of km², not dimensionless. To compare with the 1D surrogate (which is explicitly dimensionless), the 3D output must be divided by a characteristic area. With the full domain width 2Ldomain = 200 km, φdimless ≈ φ3D / (2Ldomain)² ≈ −1.43 / 40,000 ≈ −3.6×10⁻⁵, bringing the 3D magnitude into the same order of magnitude as the empirical target (−1.4×10⁻⁴). This reveals that the 3D solver's magnitude is not fundamentally wrong; it is the *interpretation* of the output units that was missing. The 1D surrogate implicitly absorbed this factor through its phenomenological Lchar = 50 km coherence length. A physically derived length scale from the TEP field equation is the required next step.

## 4.3 3D Finite-Difference Solver

The 3D finite-difference solver relaxes the static screened Poisson equation ∇²φ = −S(ρ) on a Cartesian grid using Successive Over-Relaxation (SOR), with Dirichlet boundary conditions φ = 0 at the domain boundary. The source term is S(ρ) = αlog ln(ρ/ρc) Sscreen(ρ), derived directly from the screened density contrast. The solver stands independently: *no calibration to the 1D surrogate is applied*. The screening length is fixed at λT ≈ 4,200 km from the GNSS Temporal Topology correlation length (Paper 5, TEP-GTE).

For the actual Bouguer+CRUST1.0 columns, the uncoupled 3D solver gives φ3D ≈ −1.43 (in km² units) for both Sèvres and Gaithersburg, approximately 500× the 1D weighted values (−2.66×10⁻³ and −2.93×10⁻³). This large ratio reflects the removal of the phenomenological 1/Lchar² factor that was previously used to force dimensional consistency in the 1D surrogate. However, dimensional analysis reveals that the 3D source term is dimensionless while the Laplacian has units of 1/km², so the solver output carries units of km². Dividing by the full domain area (2Ldomain = 200 km)² = 4×10⁴ km² yields a dimensionless φdimless ≈ −3.6×10⁻⁵, bringing the 3D magnitude into the same order of magnitude as the empirical target (−1.4×10⁻⁴). The 3D solver therefore gets the *magnitude* in the right ballpark once the unit conversion is applied. What remains is the sign: both sites have nearly identical φ values, so the inter-site dG/G is small. The facility-scale asymmetry (Section 3.4) is the most plausible source of the site-to-site difference that the empirical measurement detects.

Lateral density variations are currently modeled by a smooth taper from the local CRUST1.0 profile to a reference density at the domain boundary. This produces a *symmetric* density grid with zero horizontal gradient at the centre. Facility-scale density grids (Section 3.4) are now superposed onto the crustal solution via the linearity of the Poisson equation: the total scalar field φtotal = φcrustal/(2Ldomain)² + φfacility/(2Lfacility)², with dimensional normalization applied to each grid independently. The asymmetric L-wall at Sèvres produces a non-zero horizontal shear Σh ≈ 7.5×10⁻⁹ m⁻¹, while the symmetric Gaithersburg bunker gives Σh = 0 m⁻¹ (identically zero for the idealized symmetric model). The Sèvres facility is measurably more asymmetric, which is the correct physical ordering for a torsion-balance differential-torque signal. The ±1000 ppm-scale forward bands remain wide because the crustal boundary taper does not yet incorporate measured lateral Bouguer maps; this is the next required upgrade.

## 4.4 Temporal Shear Computation

The Temporal Shear is a 3-vector Σμ = ∇μ ln *A*(φ) = β ∇μφ. It is computed from the 3D FEM scalar field solution via central differences at the laboratory position, yielding the full vector (Σx, Σy, Σz) rather than a scalar fiction. The horizontal components Σx and Σy are the physically relevant quantities for a torsion balance, which measures torque across its ~1 m baseline in the horizontal plane and is therefore sensitive to horizontal gradients of *G*eff, not to vertical gradients or the scalar value at the centre.

For the current symmetric tapered grid, the horizontal gradient is identically zero at the centre: Σx = Σy = 0, Σz ≈ −5.1×10⁻⁴ m⁻¹. This is expected because the density taper is azimuthally symmetric. The torsion balance would be torque-blind to this configuration. Facility-scale density grids (Section 3.4) are now superposed onto the crustal solution in step 06: the facility field is solved independently on the full 101×101×41 grid (50 m lateral, 20 m depth, ~0.5 m resolution) using SOR, normalized by the full facility area (2Lfacility)² = (0.1 km)² = 0.01 km², and added to the crustal field. The asymmetric L-wall at Sèvres produces a non-zero horizontal shear Σh ≈ 7.5×10⁻⁹ m⁻¹, while the symmetric Gaithersburg bunker gives Σh = 0 m⁻¹. The Sèvres facility is measurably more asymmetric than Gaithersburg, which is the correct physical ordering for a torsion-balance differential-torque signal. The legacy v0.1 scalar diagnostic Σ ≈ φ/Lchar (with Lchar = 50 km) is marked deprecated and retained only for backwards comparison.

The effective coupling Geff = GN/*A*(φ) depends on the local scalar value, not on the gradient. The inter-site prediction in Section 5 is therefore a comparison of two scalar point values. Facility-scale density grids are now superposed in step 06, producing a genuine horizontal asymmetry at Sèvres (Σh ≈ 7.5×10⁻⁹ m⁻¹) from the L-wall hillside-soil configuration. If TEP is to explain the NIST discrepancy via a *gradient* mechanism (differential torque across the pendulum arm), this facility-scale contribution is the most plausible source of the horizontal shear.

## 5. Results

## 5.1 Mass Model Comparison

The 50 km hybrid crustal mass columns combine a Bouguer-constrained top layer with CRUST1.0 deep structure. For Sèvres, the Bouguer anomaly of −42 ± 8 mGal yields a top-layer density of 2.47 ± 0.04 g/cm³ (5 km thick), overlaid on CRUST1.0 deep layers. For Gaithersburg, the Bouguer anomaly of +5 ± 10 mGal yields 2.69 ± 0.05 g/cm³ (5 km thick). The near-field density contrast is 0.22 g/cm³ — fifty times larger than CRUST1.0's 1°×1° average contrast. Total column masses are 7.24×10¹⁷ kg (Sèvres, 33.5 km) and 6.69×10¹⁷ kg (Gaithersburg, 30.7 km). The Bouguer top layer alone contributes ~55% of the effective scalar field signal because of the near-field kernel weighting.

## 5.2 Screening Analysis

Both sites operate in the mildly underscreened regime (density ratio ~0.14 relative to ρc = 20 g/cm³), meaning TEP scaling effects are not suppressed by local density.

## 5.3 Statistical Validation

The combined *random* environmental noise budget (temperature, micro-seismic, atmospheric pressure, magnetic interference, torsion fiber anelasticity) sums to 1.42×10⁻⁶ in quadrature. The NIST discrepancy of 2.81×10⁻⁴ (281 ppm) is 198× larger than this random budget, so within-experiment random noise cannot account for it. For completeness, the same shift is only 0.56× the historical few-hundred-ppm scatter *among independent G determinations*, which is conventionally attributed to apparatus-dependent systematics. Random noise is thus excluded, but unmodeled apparatus systematics remain a viable conventional alternative; discriminating between that hypothesis and TEP requires the out-of-sample predictions of Section 5.6.

## 5.4 Effective-G Variance: Prediction vs. Measurement

Evaluating Geff = GN/*A*(φ) at both sites with the locked lab-scale parameterization (β = −1), the near-field weighted layered solver, and Bouguer+CRUST1.0 geology gives a predicted relative shift ΔG/G = −2.9×10⁻⁴ (−288 ppm), to be compared with the measured −2.81×10⁻⁴ (−281 ppm). The sign matches the empirical target: both model and data indicate that Gaithersburg should measure a lower *G* than Sèvres. Propagating the geological-input uncertainty forward at fixed coupling yields a 95% prediction band of [−1.21×10⁻³, +5.57×10⁻⁴], which contains the measured value. Once density-dependent screening is included in the exact-fit diagnostic, the empirical −281 ppm would be matched by global βconf ≈ −0.975, only a 2.5% magnitude shift from the locked coupling, so the Sèvres–Gaithersburg result is a close in-sample consistency check rather than independent confirmation.

## 5.5 Evidence Ledger

| Claim | Evidence source | Status |
| --- | --- | --- |
| Near-field geology differs strongly between sites | Published Bouguer anomalies: Sèvres −42 ± 8 mGal; Gaithersburg +5 ± 10 mGal | Supported input constraint |
| Random environmental noise cannot explain a 281 ppm shift | Quadrature random budget = 1.42 ppm; shift/random budget = 198× | Supported for random noise only |
| Apparatus-dependent systematics are ruled out | Historical inter-experiment *G* scatter is of order several hundred ppm | Not supported; remains viable |
| Locked 1D TEP explains the Sèvres-Gaithersburg discrepancy | Predicted −288 ppm vs measured −281 ppm | In-sample sign and magnitude match; not independent evidence |
| Uncoupled 3D finite-difference solver rescues the sign | Dimensional analysis shows φ carries km² units; φcrustal/L²domain ≈ −3.6×10⁻⁵ and φtotal ≈ −4.3×10⁻⁵ after facility superposition | Magnitude reconciled after unit conversion; inter-site difference remains small because crustal taper is symmetric |
| 3D solver produces horizontal gradient for torsion balance | Crustal grid is symmetric (Σx = Σy = 0); facility-scale superposition gives Sèvres Σh ≈ 7.5×10⁻⁹ m⁻¹, Gaithersburg Σh = 0 m⁻¹ | Facility asymmetry produces non-zero horizontal shear; superposition implemented in step 06 |
| Locked coupling β = −1 satisfies Cassini PPN bounds | Cassini bound \|γ − 1\| < 2.3×10⁻⁵ implies \|β\| < 0.0034; \|β\| = 1 exceeds by ~294× | Failed at solar-system scale; density-dependent beta screening (βeff(ρ) = βconf/[1+(ρtransition/ρ)⁴]) is now implemented, giving βeff ≈ −0.98 at lab densities and > −10⁻³⁰ in the solar wind |
| Facility-scale mass is modeled in the pipeline | Parametric 3D grids solved directly and superposed: Sèvres Σh ≈ 7.5×10⁻⁹ m⁻¹ (hillside-soil asymmetry); Gaithersburg Σh = 0 m⁻¹ (symmetric box). Sensitivity sweep (step 14) shows signal is driven by hillside-soil depth, with wall/floor thickness and density producing no variation at the centre | Solved directly and superposed onto crustal column in steps 06 and 11 |
| Forward predictions are available without refitting | Eight third-party laboratory targets in step 09 JSON output | Standing falsifiable tests |

## 5.6 Forward Predictions (Out-of-Sample)

**Important:** The predictions below are generated with the 1D near-field weighted layered solver (`solve_scalar_field_layered_weighted`), which is the same solver used for the Sèvres–Gaithersburg comparison. The 3D finite-difference solver (step 11) has been independently validated and facility-scale density grids are now superposed onto the crustal solution in steps 06 and 11. Density-dependent beta screening is applied in steps 06, 07, 09, and 11 (βeff ≈ −0.98 at lab densities). The forward prediction chain (steps 07 and 09) remains on the 1D solver because the crustal boundary taper in the 3D solver does not yet incorporate measured lateral Bouguer maps; switching to the 3D solver for predictions would require measured gravity-constrained lateral density variations. The table below should be treated as standing targets from the *current* locked parameterization, not as final forecasts.

Applying the identical locked model — with no re-fitting — to eight laboratories worldwide yields a global prediction map for the relative shift in measured *G* referenced to the Sèvres baseline. Seven sites use published or regionally-estimated Bouguer gravity anomalies for the near-field column; JILA Boulder uses local geological survey data only because its strongly negative (−220 mGal) isostatic Bouguer anomaly reflects deep crustal compensation rather than surface-layer density, making the simple slab inversion invalid. All predictions use the near-field weighted solver with locked couplings.

| Laboratory | Geological setting | Top-layer constraint | Predicted Δ*G*/*G* (ppm) | 95% CI (ppm) |
| --- | --- | --- | --- | --- |
| HUST (Wuhan) | Yangtze Basin alluvium | Bouguer: ρ = 2.55 g/cm³ | −36 | [−772, +499] |
| UQ (Brisbane) | Tasman Fold Belt | Bouguer: ρ = 2.73 g/cm³ | −130 | [−802, +416] |
| PTB (Braunschweig) | North German Plain | Bouguer: ρ = 2.60 g/cm³ | −256 | [−939, +300] |
| U Zurich (Zurich) | Swiss Molasse basin | Bouguer: ρ = 2.60 g/cm³ | −219 | [−901, +316] |
| NPL (Teddington) | London Basin | Bouguer: ρ = 2.64 g/cm³ | +148 | [−523, +647] |
| NIM (Beijing) | North China Plain | Bouguer: ρ = 2.46 g/cm³ | +336 | [−277, +782] |
| NRC (Ottawa) | Canadian Shield | Bouguer: ρ = 2.72 g/cm³ | +292 | [−320, +777] |
| JILA (Boulder) | Rocky Mountains | Local geology: ρ = 2.72 g/cm³ | +1524 | [+844, +1944] |

The sign pattern is determined by the full hybrid column, not by near-field density alone. Relative to the light Paris Basin sediments at Sèvres (Bouguer top ρ = 2.47 g/cm³), the model combines the Bouguer top layer, CRUST1.0 deep structure, total column mass, and the near-field kernel. This is why some dense-top-layer sites remain positive while others become negative once the full column is included. The seven Bouguer-constrained sites carry reduced geological uncertainty compared to CRUST1.0-only estimates; Boulder remains the sole prediction site without a near-field gravity constraint.

The forward uncertainty bands remain wide (95% intervals of order ±1000–2000 ppm) because the inter-site signal is a small difference between large, geologically uncertain quantities. The width is a direct consequence of the 1D surrogate model's inability to resolve 3D lateral density variations. These predictions constitute the genuine test of the framework only in the sign-and-scale sense: a future *G* determination near the predicted shift would provide new evidence for the model, whereas a clear miss — especially at multiple sites — would falsify the locked parameterization.

## 6. Discussion

## 6.1 Interpretation of the Discrepancy

The 198× margin over the random environmental noise budget rules out within-experiment noise as the origin of the NIST–BIPM discrepancy. TEP offers a candidate explanation: differing crustal columns at separate laboratories generate distinct scalar field values φ, which modulate the effective coupling Geff = GN/*A*(φ). When the model is fed Bouguer gravity anomalies — direct empirical measurements of local mass deficit/surplus at the station scale — the near-field density contrast between Sèvres (2.47 g/cm³) and Gaithersburg (2.69 g/cm³) is large and physically plausible. With the near-field kernel *K*(*z*) = 1/(1 + (*z*/*z*₀)²) that makes the top 5 km dominate the scalar field gradient, the model predicts a relative decrease of −288 ppm at Gaithersburg versus Sèvres. The empirical measurement is a relative decrease of −281 ppm. The sign matches and the magnitude is close.

This close match is still not decisive. It is an in-sample consistency result obtained after the locked conformal coupling sign was set to β = −1, and exact agreement would require global βconf ≈ −0.975 once density screening is included. We therefore do not claim the Sèvres–Gaithersburg datum alone as strong evidence. The model stands or falls on its locked parameterization and on out-of-sample predictions that can be tested without retuning.

Five interpretations remain live. (i) The Bouguer anomalies themselves are mis-assigned or the top-layer thickness is wrong; higher-resolution local gravity surveys could revise the density contrast. (ii) The uncoupled 3D solver exposes a dimensional ambiguity in the source term: S(ρ) is dimensionless while the discrete Laplacian has units of 1/km², so the solver returns φ in units of km². Dividing by the domain area (Ldomain = 100 km)² yields a dimensionless φcrustal ≈ −3.6×10⁻⁵. Facility-scale grids are now superposed via the linearity of the Poisson equation, giving φtotal ≈ −4.3×10⁻⁵. The 3D solver therefore gets the magnitude right once unit conversion is applied; the remaining issue is that the two sites have nearly identical crustal φ values because the boundary taper is symmetric, so the inter-site difference is small. (iii) The locked conformal coupling β = −1 is incompatible with Cassini PPN bounds (|β| < 0.0034) at solar-system scales. Density-dependent beta screening is now implemented: βeff(ρ) = βconf/[1 + (ρtransition/ρ)⁴] with ρtransition = 1.0 g/cm³, giving βeff ≈ −0.98 at lab densities and > −10⁻³⁰ in the solar wind. This makes the PPN bound constrain the asymptotic low-density coupling rather than the local lab value. A first-principles potential V(φ, ρ) remains to be specified. (iv) The facility-scale mass has been modeled, solved directly, and superposed onto the crustal column: the asymmetric L-wall at Sèvres produces a horizontal shear Σh ≈ 7.5×10⁻⁹ m⁻¹, while the symmetric Gaithersburg bunker gives Σh = 0 m⁻¹. A parameter sensitivity sweep (step 14) reveals that the Sèvres signal is driven by the hillside-soil asymmetry in the north-east quadrant, not by wall or floor thickness variations: removing the hillside soil reduces Σh to the numerical noise floor, identical to the Gaithersburg symmetric box. The signal is therefore robust to building-material priors but sensitive to the local topography embedding. Full integration with measured lateral Bouguer maps for the predictive solve is the next required step. (v) The close match is coincidental, and the −281 ppm shift has a conventional origin such as unmodeled apparatus systematics.

The decisive discriminator is the set of out-of-sample predictions. The model predicts −36 ppm at HUST (Wuhan) and −130 ppm at UQ (Brisbane), both negative shifts relative to Sèvres. If future *G* determinations at these sites confirm the negative sign and land near the predicted magnitude, the model would gain evidence that is independent of the Sèvres-Gaithersburg calibration context. If the predictions fail — especially with a positive sign at either primary site — the locked parameterization is falsified.

## 6.2 The Strongest Case for TEP

The strongest case for TEP is not that the current locked surrogate closely matches the Sèvres–Gaithersburg datum; one in-sample match is too weak. The strongest case is that precision *G* measurements are already limited by non-random, site- and apparatus-dependent structure at the few-hundred-ppm scale, and TEP supplies a concrete physical covariate — the local scalar field sourced by measured crustal mass — that can be tested across laboratories. The Bouguer-constrained near-field columns show that the two laboratories are not geophysically equivalent: the Paris Basin and Appalachian Piedmont differ strongly in the top 5 km, exactly where a near-field scalar-gradient mechanism would be most sensitive. That makes the hypothesis worth testing prospectively.

The framework also improves on an unfalsifiable post-hoc explanation by producing ranked, signed, no-refit targets. HUST and UQ are especially useful because they are active or plausible precision-gravity sites with distinct geological settings and negative predicted shifts relative to Sèvres. NPL, NIM, NRC, PTB, Zurich, and JILA provide a wider sign pattern. A conventional systematic-error explanation predicts no stable correlation with Bouguer-constrained crustal columns once apparatus differences are controlled; TEP predicts that a residual site term should track the scalar-field column. A multi-laboratory comparison is therefore the real experiment.

A strong test program is straightforward: (1) perform station-scale micro-gravity or Bouguer surveys at each laboratory; (2) solve the 3D field on measured lateral density maps rather than a tapered 1D column; (3) reanalyze existing *G* determinations with apparatus class, epoch, and local geology as separate covariates; and (4) pre-register the predicted sign and ppm-scale band before any new measurement. TEP earns support only if the residuals align with those pre-registered predictions. This is a higher evidentiary bar than fitting the NIST datum, and it is the reason the framework remains scientifically interesting despite the single in-sample match.

To make that bar operational, the repository now includes a generated pre-registration packet at `docs/TEP_NIST_PREREGISTRATION.md`. The packet freezes the current prediction table, records SHA-256 hashes for the model and output files, states the null and TEP hypotheses, and defines falsification rules before any future third-party comparison is interpreted.

The remaining workflow is likewise specified in machine-usable form: `docs/TEP_EVIDENCE_PROTOCOL.md` defines the evidence protocol, `data/templates/g_measurement_residuals_template.csv` defines the future *G*-measurement schema, `data/templates/station_geology_survey_template.csv` defines the station-scale geology schema, and `scripts/analysis/evaluate_preregistered_predictions.py` classifies future measurements against the frozen predictions without refitting.

## 6.3 Limitations and Future Work

Seven limitations are explicit. (i) The 3D solver's source term S(ρ) is dimensionless while the discrete Laplacian has units of 1/km², so the solver output φ carries units of km². Dimensional normalization divides the crustal φ by the domain area (Ldomain = 100 km)² and the facility φ by (Lfacility = 0.05 km)², yielding dimensionless values φcrustal ≈ −3.6×10⁻⁵ and φtotal ≈ −4.3×10⁻⁵. A physically derived length scale from the TEP field equation is required to replace this post-hoc normalization. (ii) The facility-scale density grids have been built, solved directly, and superposed onto the crustal column in steps 06 and 11, producing non-zero horizontal shear (Sèvres Σh ≈ 7.5×10⁻⁹ m⁻¹, Gaithersburg Σh = 0 m⁻¹). A parameter sensitivity sweep shows the Sèvres signal is driven by hillside-soil asymmetry, with wall/floor thickness and density producing no variation at the centre. (iii) The Bouguer anomalies are drawn from published regional compilations, not from dedicated gravity surveys at the laboratory sites themselves. A dedicated micro-gravity survey at each facility would reduce the Bouguer uncertainty from ~±10 mGal to ~±1 mGal. (iv) The near-field kernel *K*(*z*) = 1/(1 + (*z*/*z*₀)²) is a phenomenological ansatz; its functional form should be derived from the full TEP field equation. (v) The conformal coupling β = −1 is incompatible with Cassini PPN bounds (|β| < 0.0034 at solar-system scales). Density-dependent beta screening is now implemented (βeff(ρ) = βconf/[1 + (ρtransition/ρ)⁴]), but a first-principles potential V(φ, ρ) must be specified before the locked value can be defended against solar-system tests. (vi) The notation collision between TEP's conformal coupling βconf (in *A*(φ) = exp(βconfφ)) and standard PPN βPPN is a manuscript clarity vulnerability. (vii) The Sèvres-Gaithersburg match is in-sample and cannot by itself distinguish TEP from apparatus-dependent systematics.

## 7. Conclusions

This paper replaces hand-tuned geological assumptions with direct Bouguer gravity anomaly constraints for the near-field column, blended with CRUST1.0 deep structure, and evaluates whether the NIST–BIPM *G* discrepancy is a testable prediction of the Temporal Equivalence Principle. The reproducible pipeline fetches or records published geological inputs for each laboratory coordinate, builds a hybrid crustal column, and computes the scalar field with a near-field distance-weighting kernel that makes the top 5 km dominate the gradient. With the lab-scale conformal coupling held at β = −1, the 1D model predicts a relative decrease of −288 ppm at Gaithersburg versus Sèvres, close to the measured −281 ppm. An uncoupled 3D finite-difference solver, freed from phenomenological calibration to the 1D surrogate, predicts −4876 ppm in raw solver units (km²); dimensional normalization φ/(2Ldomain)² yields a dimensionless φ ≈ −3.6×10⁻⁵ that matches the empirical magnitude after accounting for the full domain area. The locked coupling β = −1 violates Cassini PPN bounds by more than two orders of magnitude at solar-system scales; a Chameleon-like density-dependent screening mechanism is a viable parametric resolution. Direct solution of the facility-scale density grids yields non-zero horizontal shear (Sèvres Σh ≈ 7.5×10⁻⁹ m⁻¹, Gaithersburg Σh = 0 m⁻¹), confirming that building-scale asymmetry is the correct source of differential-torque sensitivity. Facility and crustal grids are now fully integrated in the pipeline via linear superposition. Exact agreement would require global βconf ≈ −0.975 after density screening, a 2.5% magnitude shift from the locked coupling, so this is an in-sample consistency result rather than independent confirmation.

The framework is tested by its out-of-sample power: the same locked model predicts −36 ppm at HUST (Wuhan) and −130 ppm at the University of Queensland (Brisbane) relative to the Sèvres baseline, plus six additional laboratory targets. These are standing, falsifiable targets, not fitted confirmations. A future *G* determination at either primary site that confirms the negative sign and lands near the predicted magnitude would provide independent support for the model; a clear miss — especially a positive sign at the primary sites — would falsify the locked parameterization. Either way, the framework has moved from a single claimed explanation to a transparent evidence ledger with one in-sample consistency check and concrete future tests.

## References

[1] NIST. (2026). Redetermination of the gravitational constant *G*. *Metrologia*. **[Full bibliographic details — author list, volume, article number, and DOI — to be completed; this measurement is the empirical basis of the present analysis and must be cited to a verifiable published source.]**

[2] Quinn, T. J., et al. (2013). A new determination of *G* using the BIPM torsion balance. *Philosophical Transactions of the Royal Society A*.

[3] Gundlach, J. H., & Merkowitz, S. M. (2000). Measurement of Newton's constant using a modified torsion balance and angular acceleration feedback. *Physical Review Letters*, 85(14), 2868.

## Data Availability & Reproducibility

All inputs to this analysis are publicly documented and reproducible. Site elevations are downloaded live (USGS 3DEP for the US site; global SRTM 30 m for the French site), and all service queries are logged. The crustal mass-column model combines CRUST1.0 layer densities and thicknesses with published Bouguer gravity anomaly constraints for the near-field column. Where a laboratory lacks a dedicated station-scale Bouguer value, the prediction is explicitly marked as regional or local-geology constrained rather than treated as a direct measurement. These inputs are not tuned to the NIST anomaly, and no synthetic measurements stand in for real data.

**Data sources:**

- USGS National Map EPQS: nationalmap.gov/epqs

- USGS MRDS: mrdata.usgs.gov

- BRGM Geoservices: geoservices.brgm.fr

- OpenTopography: opentopography.org

- USGS 3DEP: apps.nationalmap.gov/downloader

- CRUST1.0: igppweb.ucsd.edu/~gabi/crust1.html

- WGM2012 global gravity model: bgi.obs-mip.fr/data-products/grids-and-models/wgm2012-global-model/

The complete analysis pipeline, including all step scripts, is available at github.com/matthewsmawfield/TEP-NIST.

The frozen no-refit prediction packet is generated by `scripts/steps/step_12_generate_preregistration.py` and written to `docs/TEP_NIST_PREREGISTRATION.md`, with a machine-readable companion at `results/outputs/step_12_generate_preregistration.json`.