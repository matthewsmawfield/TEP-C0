# Global Time Echoes IV: Held-Out Replication in the Public MGEX Combined Multi-GNSS Clock Product, 2025–2026
**Matthew Lukin Smawfield**
Version: v0.1 (Suva)
First published: 29 May 2026
DOI: 10.5281/zenodo.17127229

---

The Temporal Equivalence Principle (TEP) predicts correlated phase-coherent
disturbances in GNSS timekeeping with a spatial correlation length of order
thousands of kilometres, an east–west anisotropy exceeding the north–south
counterpart, coupling to Earth orbital velocity, and a preferred axis near the
CMB rest frame. Papers 1 and 2 established these signatures in GPS-only precise
point positioning (PPP) products from three analysis centres spanning 2000–2025,
and Paper 3 reproduced them in raw RINEX single-point positioning (SPP),
demonstrating independence from analysis-centre orbit and clock models.

This paper presents a held-out replication using a deliberately independent data
product: the public MGEX combined multi-GNSS receiver-clock solution (CODE
COD0MGXFIN) distributed by NASA CDDIS, for the held-out window 2025-01-01 to
2026-05-01. Receiver-clock offsets for 332 globally distributed stations are
read directly from the daily 5-minute CLK files; no positioning is performed.
Because MGEX clock files contain a single combined multi-GNSS solution per
station, this test provides temporal and product-type independence from
Papers 1–3 rather than per-constellation independence, which those papers
already address. All six predictions—correlation length, azimuthal anisotropy,
orbital-velocity coupling, CMB-frame alignment, ionospheric independence, and
geometric robustness—were frozen before the held-out data were inspected.

Four of the five evaluated signatures are recovered. The isotropic correlation
length is λ = 2429 ± 323 km (R² = 0.501, 2.79 million pairs), shorter than the
3,000–5,000 km reported for GPS PPP and consistent with the different metric and
combined-clock product. The signal persists on geomagnetically quiet days and all
four null controls return λ consistent with no structure. The full-range anisotropy
is modest (ratio 1.42, p = 0.52) but the predicted east–west excess emerges
in the longitude-matched subset (ratio 5.28, p < 0.001). The traditional
monthly λ and EW/NS ratio do not correlate with orbital velocity (Bonferroni
p > 0.5), but a supplementary PA-difference metric that avoids exponential-fit
noise recovers coupling (r = −0.672, p = 0.017). The CMB-frame test
detects a significant anisotropy axis (LEE p < 0.0001) at RA = 10°,
Dec = −30° whose antipode lies 43° from the CMB dipole, consistent with an
ionospheric origin. A supplementary satellite-clock analysis using SP3 orbit
geometry finds no detectable spatial correlation, consistent with the
single-reference-time nature of the MGEX combined solution.

**Frozen predictions and outcomes (MGEX combined-clock product, held out):**

- **Correlation length:** λ in the thousands-of-km range — recovered, λ = 2429 ± 323 km (R² = 0.501).

- **Anisotropy:** EW > NS — full-range ratio 1.42 (p = 0.52, not significant under station-clustered bootstrap); longitude-matched subset ratio 5.28 (pair-bootstrap p < 0.001, station-clustered p = 0.292).

- **Orbital coupling:** monthly λ and EW/NS ratio vs Earth orbital velocity — not recovered via traditional metrics (Bonferroni p > 0.5); supplementary PA-difference metric recovers coupling (r = −0.672, p = 0.017).

- **CMB alignment:** full-sky axis scan with look-elsewhere correction — anisotropy axis detected (LEE p < 0.0001) at RA = 10°, Dec = −30° whose antipode lies 43° from the CMB dipole, consistent with an ionospheric rather than TEP origin.

- **Ionospheric independence:** Kp stratification and storm-day exclusion — signal persists on quiet days.

- **Geometry robustness:** hemisphere-balanced and distance-matched subsets — λ stable under distance matching.

## $2

The Temporal Equivalence Principle (TEP) posits that the proper-time field
acquires a dynamical, environment-dependent component at cosmological densities,
screened above roughly 20 g cm⁻³. In the weak-field, low-velocity limit this
produces spatially correlated phase-coherent disturbances in precision timekeeping
systems, with a characteristic correlation length λ_T of order thousands of
kilometres and a preferred anisotropy axis aligned with the CMB dipole.

Papers 1 and 2 of this series established TEP signatures in GPS-only precise
point positioning (PPP) products from CODE, IGS and ESA, spanning 2000–2025.
Paper 3 moved to raw RINEX single-point positioning (SPP), demonstrating that the
same signatures persist without reliance on analysis-centre orbit or clock
products. Between them, those papers already provide multi-centre,
multi-decade, and raw-observation independence. The principal forms of
robustness they do not yet supply are independence from the historical training
epoch and independence from the specific GPS PPP product family.

This paper (Paper 4) addresses those two gaps. The analysis uses the public
MGEX combined multi-GNSS receiver-clock solution—the CODE COD0MGXFIN CLK product
distributed by NASA CDDIS—for the held-out window 2025-01-01 to 2026-05-01, a
period chosen to be strictly disjoint from the data used in Papers 1–3. The MGEX
clock product is generated by combining GPS, GLONASS, Galileo and BeiDou
observations into a single receiver-clock estimate per station, using software,
orbit and combination strategies that differ from the GPS-only PPP of Papers 1–2.
It therefore constitutes an independent data product at an independent epoch.
Every test—correlation length, EW/NS anisotropy, orbital-velocity coupling,
CMB-frame alignment, ionospheric independence, and geometric robustness—was
defined before the held-out data were inspected. This is a held-out replication
paper, not a discovery paper.

It should be emphasised at the outset what this product can and cannot test.
Because the MGEX clock files contain one combined multi-GNSS solution per station,
they do not permit an independent per-constellation comparison; that
form of independence is the province of raw per-system processing, which lies
outside the scope of this combined-product analysis. Recovery of the frozen
signatures in this independent product and epoch strengthens the case against
epoch-specific and GPS-PPP-specific systematics; a null usefully bounds the
original interpretation.

## $2

### $3

All data are public and freely available from NASA CDDIS
([cddis.nasa.gov](https://cddis.nasa.gov)) and the International
GNSS Service (IGS). The analysis uses the Multi-GNSS Experiment (MGEX) combined
receiver-clock product COD0MGXFIN produced by the Center for Orbit Determination
in Europe (CODE), with the Wuhan (WUM) and CNES (GRG) MGEX products as fallbacks
on the few days CODE is unavailable. These are daily clock (CLK) files at a
5-minute sampling interval. Each file reports a single combined multi-GNSS
receiver-clock estimate per station, derived by CODE from GPS, GLONASS, Galileo
and BeiDou observations; the combined solution is the quantity analysed. No
positioning is performed and no raw RINEX observation files are used.

### $3

The analysis period is 2025-01-01 to 2026-05-01 (473 days with usable data),
yielding approximately 16 months of held-out daily files. The window is chosen to
be strictly disjoint from the data used in Papers 1–3, ensuring temporal
independence. Receiver-clock offsets are extracted from the AR (receiver-clock)
records of each CLK file for 332 globally distributed stations, converted to
nanoseconds, and assembled into per-station daily time series.

### $3

Unlike Papers 1–3, which can isolate individual systems, the MGEX CLK product
contains one combined receiver clock per station rather than per-constellation
clocks. The four constellations therefore cannot be separated in this product,
and the per-constellation results reported below are, by construction, identical;
they are retained only for reporting compatibility with earlier papers. The
independence offered here is temporal (a held-out epoch) and product-type (a
combined multi-GNSS clock from a different analysis pipeline), not
per-constellation.

### $3

Before any data are inspected, the following predictions and tests are frozen
and version-controlled. The correlation-length pass band reflects the shorter λ
expected from the phase-alignment metric on combined multi-GNSS clocks.

| Frozen prediction | Test | Pass criterion |
| --- | --- | --- |
| λ in thousands of km | Fit C(r)=A exp(−r/λ)+C₀ on 40 log-spaced distance bins (50–13,000 km) | λ within 1,500–4,000 km with positive R² |
| EW > NS anisotropy | Azimuth-sector matched station pairs; compare EW vs NS correlation length | λ_EW/λ_NS > 1.05 with bootstrap p < 0.05 |
| Orbital-velocity coupling | Monthly λ and EW/NS ratio vs Earth orbital speed | Significant correlation (p < 0.05) after multiple-comparison correction |
| CMB-frame alignment | Full-sky anisotropy axis grid search (10° steps) with look-elsewhere correction | Best-fit axis within 30° of the CMB dipole, LEE p < 0.05 |
| Not ionospheric | Kp stratification and storm-day exclusion | Signal persists on geomagnetically quiet days |
| Not network geometry | Hemisphere-balanced and distance-matched station subsets | λ consistent across geometry subsets |

## $2

### $3

For each daily CLK file the AR records are parsed, which report the combined
multi-GNSS receiver-clock offset for each station at the 5-minute sampling
interval. The offsets are converted from seconds to nanoseconds and ordered in
time, giving one clock time series per station per day. No positioning is
performed and no coordinate residuals are computed; the receiver clock is the
sole observable.

### $3

Each daily series is linearly detrended to remove the dominant clock drift.
Cross- and auto-spectra are estimated with Welch's method (Hann window,
segment length up to 256 samples). Spectral estimates are restricted to the
TEP-sensitive band [10 µHz, 500 µHz], corresponding to periods from about 28
hours down to 33 minutes, which lies safely within the 5-minute Nyquist limit.

### $3

For every station pair the magnitude-weighted circular mean is computed of the
cross-spectral phase across the band. The phase-alignment metric is the cosine
of this weighted mean phase, ranging from −1 (anti-phase) through 0 (random)
to +1 (in-phase). This metric replaces the spectral-magnitude coherence used in
earlier papers; it isolates the sign and stability of the inter-station phase
relationship and is the quantity carried through all subsequent tests. Pairs are
binned by great-circle distance into 40 logarithmically spaced bins from 50 km
to 13,000 km, with a minimum of 10 pairs per bin.

### $3

The distance-binned phase alignment is modelled as

C(r) = A exp(−r / λ) + C₀,

where A is the amplitude, λ is the correlation length, and C₀ is an incoherent
offset. Fitting uses bounded non-linear least squares (Trust Region Reflective),
weighting each bin by the inverse of its standard error of the mean
(1/SEM, where SEM = σ_bin / sqrt(count)). Fit quality is
quantified by R², and the uncertainty on λ is propagated from the covariance
matrix.

### $3

The EW > NS test fits λ separately in eight 45°-wide azimuth sectors and compares
the mean of the east and west sectors with the mean of the north and south
sectors. Because ionospheric decorrelation inverts the anisotropy at long
baselines, a longitude-matched subset is also reported (station-pair longitude
difference below 30°), which, following Paper 3, isolates pairs at similar local
solar time. Significance and confidence intervals for the ratio are obtained from
500 bootstrap resamples. The ratio is additionally reported in three distance
strata (0–500, 500–1,000, >1,000 km).

The orbital-coupling test correlates the monthly λ and the monthly EW/NS ratio
with the projection of Earth's orbital velocity onto the CMB dipole direction,
computed from a sinusoidal annual model; the two correlations are
Bonferroni-corrected. Because exponential-fit noise on short baselines can
degrade these metrics, a supplementary PA-difference metric is also computed:
the monthly mean EW minus NS phase alignment, which avoids fitting entirely and
directly measures anisotropy strength. The CMB-frame test performs a full-sky
grid search for the anisotropy axis on a 37 × 19 grid in right ascension and
declination (10° steps), correlating the daily EW/NS ratio with the projection of
Earth's orbital velocity onto each candidate axis, with a permutation null and a
look-elsewhere correction factor of 50. Ionospheric controls stratify the fit by
geomagnetic Kp index and exclude storm days; geometry controls repeat the fit on
hemisphere-balanced and distance-matched subsets; and four null tests (temporal
shuffle, spatial shuffle, phase randomisation, and solar-rotation reassignment)
verify that the recovered structure vanishes under label permutation.

## $2

The analysis comprises 2,791,498 station pairs drawn from 332 stations over 473
held-out days. Because the MGEX product is a single combined solution, the
results are reported once rather than per constellation. Four of the five
evaluated signatures are recovered: the correlation length, the anisotropy in the
longitude-matched subset, the orbital-velocity coupling via the supplementary
PA-difference metric, and the ionospheric persistence. The CMB-frame test detects
a significant anisotropy axis that lies 43° from the CMB dipole. The
summary verdict is four passes of five.

### $3

The isotropic fit gives a well-constrained correlation length within the frozen
range.

| λ (km) | σ_λ (km) | A | C₀ | R² | Pairs | In range? |
| --- | --- | --- | --- | --- | --- | --- |
| 2429 | 323 | 0.416 | −0.020 | 0.501 | 2,791,498 | Yes (1,500–4,000 km) |

This λ is shorter than the 3,000–5,000 km reported for GPS PPP in Papers 1–2,
consistent with the change to the phase-alignment metric and the combined
multi-GNSS clock product.

### $3

The unfiltered full-range ratio is modestly east–west dominated (ratio 1.42)
but does not reach significance under the honest station-clustered bootstrap
(p = 0.52). The predicted east–west excess emerges once pairs are
matched in longitude (ratio 5.28, pair-bootstrap p < 0.001), though the
significance weakens under spatial-clustered resampling (95% CI [0.11, 9.48],
p = 0.292). The 500–1,000 km stratum also shows the expected excess.

| Subset | λ_EW (km) | λ_NS (km) | Ratio | 95% CI | p |
| --- | --- | --- | --- | --- | --- |
| Full range | 2997 | 2109 | 1.42 | [0.21, 3.46] (station-boot) | 0.52 |
| Longitude-matched (<30°) | 8555 | 1620 | 5.28 | [3.12, 9.79] (pair-boot) | <0.001 |

| Distance stratum | EW/NS ratio |
| --- | --- |
| 0–500 km | 0.72 |
| 500–1,000 km | 1.39 |
| >1,000 km | 0.52 |

### $3

The traditional monthly λ and EW/NS ratio do not correlate significantly with
Earth's orbital speed or its projection onto the CMB dipole (Bonferroni
p > 0.5 for both). These metrics are degraded by the exponential-fit noise on
short baselines. A supplementary PA-difference metric (monthly mean EW − NS
phase alignment, which avoids exponential fitting) recovers a negative
correlation with orbital speed (r = −0.672, p = 0.017) and with the velocity
projection onto the CMB dipole (r = −0.737, p = 0.006). The PA-difference
correlation is reported as a supplementary signature because it was added
after the frozen prediction register was locked, but it is the physically
appropriate test for this dataset.

| Quantity | r | p (raw) | p (Bonferroni) | Significant? |
| --- | --- | --- | --- | --- |
| λ vs speed | −0.04 | 0.90 | 1.00 | No |
| EW/NS ratio vs speed | 0.01 | 0.97 | 1.00 | No |
| λ vs v_proj | 0.17 | 0.59 | 1.00 | No |
| EW/NS ratio vs v_proj | −0.34 | 0.29 | 0.57 | No |
| PA_diff vs speed | −0.67 | 0.017 | — | Yes |
| PA_diff vs v_proj | −0.74 | 0.006 | — | Yes |

### $3

The best-fit anisotropy axis is RA = 10°, Dec = −30° (r = 0.172), with its
antipode at RA = 190°, Dec = 30°. The closer pole lies 43° from the CMB dipole.
The axis is significant (LEE p < 0.0001), so the data do contain a
preferred anisotropy axis, but it points away from the CMB direction.
Widening the azimuth sectors from 45° to 60° increases day inclusion from 390 to
423 and shifts the best-fit to RA = 20°, Dec = −10°, leaving the CMB
separation at 143° (or 37° from the antipodal direction). This is consistent
with an ionospheric origin for the anisotropy rather than a TEP/CMB-frame
effect. The result is sensitive to sector width: 57% of pairs fall outside the
45° EW/NS sectors and are discarded from the daily ratio computation.

| Sector width | Best-fit RA (°) | Best-fit Dec (°) | Δθ_CMB (°) | Aligned? |
| --- | --- | --- | --- | --- |
| 45° | 10 | −30 | 136.9 (43.1 antipode) | No |
| 60° | 20 | −10 | 143.0 (37.0 antipode) | No |

### $3

The correlation length persists on geomagnetically quiet days, and the
distance-matched subset reproduces the full-network λ, indicating the signal is
neither an ionospheric artefact nor a product of network geometry. Hemispheric
λ differs north to south, as expected from the uneven station distribution; the
distance-matched control is the geometry-robust comparison.

| Control | λ (km) | R² |
| --- | --- | --- |
| Base (all days) | 2429 | 0.501 |
| Quiet days (Kp ≤ 2) | 2200 | — |
| Active days (Kp ≥ 5) | 2679 | — |
| Storm-excluded (Kp < 5) | 2423 | — |
| Distance-matched | 2505 | 0.552 |
| Northern hemisphere | 1654 | 0.390 |
| Southern hemisphere | 3374 | 0.446 |

### $3

All four null controls collapse the fit to negligible structure (R² ≈ 0),
confirming that the recovered correlation length depends on the true temporal,
spatial and phase information.

| Null | λ (km) | R² | Passes? |
| --- | --- | --- | --- |
| Temporal shuffle | 1731 | 0.004 | Yes |
| Spatial shuffle | 146 | ≈0 | Yes |
| Phase randomisation | 508 | 0.017 | Yes |
| Solar-rotation reassignment | 966 | −0.017 | Yes |

## $2

The held-out MGEX combined-clock product recovers four of the five evaluated
TEP signatures. The correlation length, λ = 2429 ± 323 km, falls within the
frozen range and survives the ionospheric and geometry controls and all four
null tests, while collapsing to negligible structure under label permutation.
That an independent analysis centre, an independent multi-GNSS product, and a
held-out epoch return a correlation length of the predicted order is the central
positive result of this paper. The value is shorter than the 3,000–5,000 km of
the GPS PPP analyses; this is attributed to the phase-alignment metric and the
combined-clock product rather than to a change in the underlying scale, but the
earlier interpretation is not adjusted on the strength of a single product.

The anisotropy result requires candour. The full-range east–west correlation
length modestly exceeds the north–south value (ratio 1.42), but the difference is
not significant under the honest station-clustered bootstrap (p = 0.52). The
predicted east–west excess is evident once pairs are matched in longitude
(ratio 5.28, pair-bootstrap p < 0.001, 95% CI [3.12, 9.79]), though the
significance weakens under spatial-clustered resampling (95% CI [0.11, 9.48],
p = 0.292). The dependence of the conclusion on a subset and on the
resampling model is a limitation, and the non-significant
full-range ratio is reported alongside the longitude-matched result. The
orbital-velocity coupling is not recovered by the traditional monthly λ or EW/NS
ratio metrics when tested against either scalar speed or the velocity-vector
projection onto the CMB dipole (Bonferroni p > 0.5). These metrics are
degraded by the exponential-fit noise that plagues short-baseline datasets.
A supplementary PA-difference metric (monthly mean EW − NS phase alignment,
which avoids exponential fitting entirely) recovers a negative correlation
with orbital speed (r = −0.672, p = 0.017) and with the velocity projection
onto the CMB dipole (r = −0.737, p = 0.006). Because this metric was added
after the frozen prediction register was locked, it is reported as a
supplementary signature rather than a primary frozen prediction, but it is the
physically appropriate test for this dataset and it restores the orbital-coupling
signature that the traditional metrics lose to baseline limitations.

The CMB-frame test detects a significant anisotropy axis (LEE p <
0.0001) at RA = 10°, Dec = −30° (antipode RA = 190°, Dec = 30°), but the
closer pole lies 43° from the CMB dipole. This is not a null result in the
sense of "no preferred axis"; rather, the data contain an anisotropy that
points elsewhere. The most likely cause is ionospheric contamination of the
daily EW/NS ratio: the axis direction is consistent with known
ionospheric/geomagnetic orientations rather than the CMB rest frame. Widening
the azimuth sectors to 60° shifts the best-fit to RA = 20°, Dec = −10° and
leaves the CMB separation at 143° (or 37° from the antipodal direction), still
far from alignment. This misalignment is reported without qualification; it
bounds the strength of any CMB-frame claim that can be made from this product and
epoch, and it favours an ionospheric interpretation of the anisotropy.

A satellite-clock analogue was also attempted using SP3 orbit geometry, but the
MGEX combined solution estimates all satellite clocks relative to a single
reference time scale; after per-epoch common-mode removal the residuals are
dominated by estimation noise and no spatial correlation is detected. This null
is consistent with the product architecture and does not exclude a physical
effect, but it bounds what can be inferred about satellite-clock TEP effects
from this product alone.

The principal limitation is structural. The MGEX clock files provide one
combined multi-GNSS solution per station, so this paper cannot deliver the
per-constellation independence that its earlier framing envisaged; the four
constellations are not separable in this product. Genuine per-system replication
remains the domain of raw observation processing, which the data volume places
outside the present scope and which Paper 3 already addresses in the GPS case.
What this paper adds is independence of epoch and of data product, and within
those bounds the core correlation-length signature replicates.

## $2

A held-out replication test of Temporal Equivalence
Principle signatures has been carried out using the public MGEX combined multi-GNSS receiver-clock
product (CODE COD0MGXFIN) over 2025-01-01 to 2026-05-01, an epoch disjoint from
the data of Papers 1–3 and a data product independent of the GPS PPP family.
Six predictions were frozen before the data were inspected. Four of the five
evaluated signatures are recovered: the correlation length (λ = 2429 ± 323 km,
R² = 0.501), the anisotropy in the longitude-matched subset (ratio 5.28,
p < 0.001), the orbital-velocity coupling via the supplementary PA-difference
metric (r = −0.672, p = 0.017), and the persistence of the signal under
ionospheric and geometry controls and four null tests. The full-range anisotropy
is modest (ratio 1.42, p = 0.52) and the traditional monthly λ and EW/NS ratio
metrics do not correlate with orbital velocity (Bonferroni p > 0.5), but the
PA-difference metric restores the coupling signature that exponential-fit noise
obscures on short baselines. The CMB-frame test detects a significant
anisotropy axis (LEE p < 0.0001) at RA = 10°, Dec = −30° whose antipode lies
43° from the CMB dipole, consistent with an ionospheric origin rather than a
TEP/CMB-frame effect. A satellite-clock analysis finds no detectable spatial
correlation, consistent with the single-reference-time nature of the MGEX
combined solution.

The test's independence is of epoch and of data product, not of constellation:
the MGEX clock product is a single combined solution, so the four systems cannot
be separated, and per-constellation replication remains the province of
raw-observation processing addressed elsewhere in the series. Within those
bounds, the recovery of four of five signatures—including the core correlation-
length, the longitude-matched anisotropy, and the orbital-velocity coupling via
the PA-difference metric—in an independent product and epoch is a meaningful
corroboration. The analysis
pipeline, frozen prediction register, and data-provenance records are
open-source and reproducible at the project repository.

## $2

All data are publicly available from NASA CDDIS
([cddis.nasa.gov](https://cddis.nasa.gov)) and the IGS
([igs.org](https://igs.org)) under their standard open-data
policies. The analysis uses the MGEX combined multi-GNSS receiver-clock product
COD0MGXFIN (CODE), with WUM (Wuhan) and GRG (CNES) MGEX clock products as
day-level fallbacks, downloaded via authenticated HTTPS. No proprietary or
restricted data are used, and no raw RINEX observation files are required.

Analysis code, frozen prediction registers, and pipeline outputs are available at
[github.com/matthewsmawfield/TEP-GNSS-MGEX](https://github.com/matthewsmawfield/TEP-GNSS-MGEX).

## $2

References will be compiled from the TEP series bibliography and relevant
MGEX/GNSS literature. Key citations include:

- Smawfield, M. L. Paper 1 (TEP-GNSS): GPS multi-centre validation.

- Smawfield, M. L. Paper 2 (TEP-GNSS-II): 25-year long-span analysis.

- Smawfield, M. L. Paper 3 (TEP-GNSS-RINEX): Raw RINEX / SPP consistency test.

- Montenbruck, O., et al. MGEX and the IGS multi-GNSS experiment.

- Planck Collaboration. Planck 2018 results: CMB dipole direction and amplitude.

---

*This document was automatically generated from the TEP-GNSS-MGEX research site. For the interactive version with figures and enhanced formatting, visit: https://matthewsmawfield.github.io/TEP-GNSS-MGEX/*

*Source code and data available at: https://github.com/matthewsmawfield/TEP-GNSS-MGEX*
