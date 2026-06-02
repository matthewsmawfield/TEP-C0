# Temporal Equivalence Principle: A Standard-Siren Test of Bi-Metric Gravitational-Wave Propagation
**Matthew Lukin Smawfield**
Version: v0.1 (Harare)
First published: 29 May 2026

---

## Abstract

Standard ΛCDM assumes that gravitational waves and electromagnetic radiation propagate through the same effective distance-redshift relation. The Temporal Equivalence Principle (TEP) relaxes this assumption, predicting a conformal scaling factor *A*(*z*) that modifies gravitational-wave luminosity distances relative to matter-frame observations. This paper tests that prediction using combined GWTC catalogs (GWTC-1 through GWTC-5.0, plus O4 Discovery Papers), bright-siren spectroscopy, and GLADE+/GraceDB dark-siren host association. Rather than treating the analysis as a standalone precision measurement of *H*0, we ask whether the GW data prefer TEP's redshift-dependent conformal scaling over the standard ΛCDM relation. The locked lab-scale model uses *A*(*z*) = exp(*β**φ*0(1+*z*)*n*); with *β* = -1 and *φ*0 < 0, *A*(*z*) rises above unity while the magnitude of the conformal departure grows with redshift. With corrected event-level uncertainties and the Hubble fit restricted to 15 candidate-quality host associations from 77 independent-redshift events, the single-host point-estimate fit gives ΛCDM *H*0 = 81.3 km/s/Mpc and lab-fixed TEP *H*0 = 82.5 km/s/Mpc. The lab-fixed comparison gives Δχ² = +0.031 and ΔBIC = +0.031 relative to ΛCDM, which is statistically inconclusive. Bootstrap resampling gives TEP lower χ² in 90.6% of resamples but closer to Planck in only 0.5%; leave-one-out tests give lower χ² in 100% of omissions but closer to Planck in 0%. A host-marginalized dark-siren diagnostic gives Δ(-2 ln *L*) = +0.037 with the sky prior, weakly favoring TEP, and all six host-prior ablations weakly favor TEP by likelihood. The redshift split gives small positive Δχ² in both halves (low *z* = +0.019, high *z* = +0.006), so the present sample does not yet resolve a decisive redshift-growth signature. The result is best read as a reproducible sensitivity diagnostic and a falsifiable target for deeper host marginalization rather than a finalized detection or Hubble-tension resolution.

Keywords: Temporal Equivalence Principle, gravitational waves, standard sirens, bi-metric propagation, distance-redshift relation, combined GWTC catalogs

## 1. Introduction

## 1.1 The Hubble Tension and Its Interpretation

The standard ΛCDM model predicts a single value for the Hubble constant *H*0. Measurements from the cosmic microwave background (Planck, *H*0 ≈ 67.4 km/s/Mpc) and the local distance ladder (SH0ES, *H*0 ≈ 73 km/s/Mpc) differ by approximately 5σ. An extensive literature has attempted to resolve this discrepancy by measuring *H*0 with gravitational-wave standard sirens (Abbott *et al.* 2017, 2019, 2023), by proposing new physics in the early universe, or by invoking systematic errors in local calibrations. This paper takes a fundamentally different approach.

Rather than adding another *H*0 measurement to an already crowded field, we test a specific, falsifiable prediction of the Temporal Equivalence Principle (TEP): that gravitational waves propagate through a conformally scaled metric, producing a *redshift-dependent* modification to the distance-redshift relation. This modification is not degenerate with *H*0 or with dark energy; it predicts a particular functional form for how GW luminosity distances deviate from the ΛCDM expectation as a function of redshift. If confirmed, it would imply that the Hubble tension is not a discrepancy in cosmological parameters, but a signature of bi-metric propagation.

## 1.2 The TEP Bi-Metric Prediction

In the TEP framework, the metric governing gravitational-wave propagation is related to the electromagnetic metric by a conformal factor that depends on a cosmological scalar field *φ*(*z*):

\begin{equation}
d_L^{(\text{TEP})}(z) = A(z) \, d_L^{\Lambda\text{CDM}}(z; H_0)
\end{equation}

where the redshift-dependent conformal factor is

\begin{equation}
A(z) = \exp\!\bigl[\beta\,\phi(z)\bigr], \qquad \phi(z) = \phi_0 (1+z)^n .
\end{equation}

In the lab-fixed test, the parameters *φ*0 and *n* are not free: *φ*0 follows from the locked 2025 lab-scale convention (the NIST/BIPM *G* discrepancy, Paper 21), and *n* ≈ 1 follows from the matter-density scaling of the scalar field. The dimensionless conformal coupling used in this LVK implementation is *β* = -1. With *φ*0 < 0, *A*(*z*) rises above unity while the magnitude of the departure from unity grows with redshift, so the TEP distance bias is stronger at higher *z*. This produces a redshift-dependent Hubble-diagram distortion that is tested against, but not yet decisively distinguished from, ΛCDM in the current public sample.

## 1.3 This Work

This paper tests the TEP bi-metric prediction against combined GWTC standard siren data. We fit both ΛCDM and TEP distance-redshift relations to the independent-redshift sample and compare their χ2, AIC, and BIC values. The primary lab-fixed comparison has the same number of fitted parameters as ΛCDM (only *H*0 is fitted in each case). A secondary joint-fit diagnostic lets *H*0, *φ*0, and *n* vary, with an explicit information-criterion penalty for the two additional parameters.

The structure is as follows: Section 2 derives the cosmological scalar field profile from lab-scale TEP parameters; Section 3 describes the combined GWTC catalog data selection and independent-redshift methodology; Section 4 details the computation pipeline; Section 5 reports the model comparison results; Section 6 discusses the implications for the Hubble tension; and Section 7 concludes.

## 2. Theoretical Framework

## 2.1 The Conformal Scaling Factor

Under TEP, the conformal factor *A*(*z*) relates the effective metric to the background metric. For gravitational waves propagating through regions of varying scalar field strength, the effective luminosity distance is rescaled by *A*(*z*) along the propagation path. The full redshift-dependent conformal factor is

\begin{equation}
A(z) = \exp\!\bigl[\beta\,\phi_0\,(1+z)^n\bigr] .
\end{equation}

In the small-field limit relevant to the present sample, this may be linearised as

\begin{equation}
A(z) \approx 1 + \beta\,\phi_0\,(1+z)^n + \mathcal{O}(\phi^2) ,
\end{equation}

where *β* = -1 is the dimensionless conformal coupling used in this LVK implementation, *φ*0 = −0.013 is the locked lab-scale convention, and *n* = 1 follows from matter-density scaling. The linearisation is accurate to better than 0.1% over the redshift range of the current GW sample.

## 2.2 Bi-Metric Propagation

The bi-metric framework distinguishes between the metric governing electromagnetic radiation (the observed metric) and the metric governing gravitational waves (the effective metric). This distinction arises naturally from the density-dependent coupling of the scalar field.

## 2.3 Hubble Diagram Prediction

The TEP-adjusted Hubble diagram overlays the standard ΛCDM prediction, the raw combined-GWTC data points with independent redshifts, and the TEP-scaled distance-redshift relation. In the locked lab-fixed convention currently implemented, the conformal scaling shifts the candidate-quality best-fit *H*0 upward from the ΛCDM value. The decisive test is therefore not whether the shift moves toward a particular external Hubble baseline, but whether the redshift-dependent scaling improves the distance fit after using the correct per-event uncertainties and model-comparison penalty.

## 3. Data Selection

## 3.1 Combined GWTC Catalogs

We combine all publicly available LVK event catalogs queried by the pipeline from the Gravitational-Wave Open Science Center (GWOSC): GWTC-1-confident, GWTC-2, GWTC-2.1-confident, GWTC-3-confident, GWTC-4.0, GWTC-4.1, O4 Discovery Papers, and GWTC-5.0. Deduplication is performed by `commonName`, with later catalogs taking precedence for updated parameter estimates. Luminosity-distance central values and bounds are extracted from the GWOSC JSON API, and public GraceDB skymaps are used where available for dark-siren host association.

## 3.2 Precision Filtering

Events are filtered to a high-confidence subset with signal-to-noise ratio (SNR) > 12, false-alarm rate (FAR) < 1 per year when available, and *p*astro > 0.9 when available. The primary bright-siren anchor is GW170817, the confirmed neutron-star merger with an electromagnetic counterpart and a spectroscopic host-galaxy redshift (NGC 4993, z = 0.0092).

## 3.3 Independent Redshifts — Circularity Avoidance

To avoid the circularity problem that invalidates cosmological tests using GWOSC-derived redshifts (which are computed from luminosity distances assuming &Lambda;CDM), we obtain redshifts from two independent sources:

**Bright sirens.** Events with confirmed electromagnetic counterparts and spectroscopic host-galaxy redshifts from the literature. Only GW170817 satisfies this criterion in the current sample.

**Dark sirens.** For events without electromagnetic counterparts, we perform candidate host-galaxy association using public GraceDB HEALPix skymaps where available and a merged redshift-bearing galaxy list. The baseline catalog is GLADE+ from VizieR VII/291; a deeper NED cone query around the highest-probability sky region adds literature redshift objects where available. Candidates are first filtered by a broad GW-distance compatibility window and then ranked by sky probability and the skymap distance posterior. Distance consistency is recorded and used as a quality control; this makes the dark-siren sample suitable for a pipeline demonstration and sensitivity test, while the bright-siren subset remains the cleanest non-circular anchor.

GWOSC redshift fields are *never* used in the fit, avoiding direct reuse of redshifts derived from GW luminosity distances under a fiducial ΛCDM cosmology.

## 4. Computation

## 4.1 Pipeline Architecture

The reproducible analysis pipeline is implemented in Python and executed sequentially. Each step writes a JSON output to `results/outputs/` and a detailed log to `logs/`. Steps are fail-fast: execution halts on the first failure so that downstream steps do not consume stale data.

## 4.2 GR Distance Extraction

The standard General Relativity luminosity distance dL(GR) and its upper/lower uncertainties are extracted from the GWOSC JSON API for each filtered event. Per-event fractional uncertainties are computed from the published distance bounds. Independent redshifts are taken from step 02 (bright-siren spectroscopy + merged GLADE+/NED dark-siren Bayesian host association); GWOSC redshift fields are *never* used.

## 4.3 TEP Distance Transformation

The locked TEP conformal scaling factor A(φ) is used in the forward GW-distance model:

\begin{equation}
d_L^{(\text{GW})}(z) = A(z) \, d_L^{\Lambda\text{CDM}}(z; H_0)
\end{equation}

For observed GW-inferred distances, the corresponding matter-frame corrected distance is dLGR / A(z). Downstream fits use the forward model above, the propagated per-event distance uncertainties from the GWOSC bounds, and redshift uncertainty in distance space. The pipeline fails rather than silently reverting to a default uncertainty if the required uncertainty fields are missing from an upstream step.

## 5. Results

## 5.1 Hubble Diagram

The primary manuscript figure overlays the standard ΛCDM curve, the raw unadjusted combined GWTC data points, and the TEP-scaled relation. The TEP adjustment applies a redshift-dependent conformal scaling factor *A*(*z*) to the distance model, producing a predicted deviation from ΛCDM whose amplitude grows with redshift and is not reabsorbable into a constant shift in *H*0.

![Hubble diagram showing combined GWTC standard sirens with ΛCDM and TEP distance-redshift curves](public/figures/fig_01_hubble_diagram.png)

Figure 1. Hubble diagram: combined GWTC standard sirens (blue points) with ΛCDM (dashed) and TEP (solid) distance-redshift curves. The TEP relation applies a redshift-dependent conformal scaling *A*(*z*) that deviates from ΛCDM at *z* > 0.1.

## 5.2 H₀ Reconciliation

The Hubble constant is computed from the independent-redshift standard siren sample by fitting ΛCDM and TEP distance-redshift relations. Three models are compared: (1) ΛCDM with *H*0 free; (2) TEP with locked lab-scale convention *φ*0 and *n* fixed, *H*0 free; (3) TEP joint fit with *H*0, *φ*0, and *n* all free. Results are reported relative to the early-universe CMB baseline (Planck: *H*0 ≈ 67.4 km/s/Mpc) and the local distance ladder (SH0ES: *H*0 ≈ 73 km/s/Mpc).

With the corrected per-event distance uncertainties, 77 events have independent redshifts and 15 candidate-quality host associations enter the primary Hubble fit. The single-host point-estimate fit gives ΛCDM *H*0 = 81.3 km/s/Mpc and lab-fixed TEP *H*0 = 82.5 km/s/Mpc. Thus the locked TEP scaling shifts the inferred scale upward by 1.2 km/s/Mpc under the current *β* = -1 convention, moving it farther from both Planck and SH0ES while keeping the same number of fitted cosmological parameters as ΛCDM. The TEP joint fit gives *H*0 = 88.1 km/s/Mpc with best-fit *φ*0 = -0.050 and *n* = 3.0, at the fit bounds, so it is treated as a diagnostic rather than a physical parameter recovery.

![Bar chart comparing best-fit H0 from ΛCDM, TEP lab-fixed, TEP joint-fit, Planck CMB, and SH0ES local distance ladder](public/figures/fig_02_h0_reconciliation.png)

Figure 2. Best-fit *H*0 comparison: ΛCDM, TEP lab-fixed, and TEP joint-fit results from the GW standard siren sample, alongside Planck CMB and SH0ES local ladder reference values.

## 5.3 Model Comparison

Frequentist model comparison (&chi;&sup2;, AIC, BIC) is performed between ΛCDM and TEP bi-metric as competing hypotheses. The joint TEP fit incurs a BIC penalty of 2 ln(*N*) for its two additional free parameters (*φ*0, *n*) and must improve &chi;&sup2; by more than this to be preferred. A full Bayesian analysis with posterior samples from emcee MCMC provides credible intervals on (*H*0, *φ*0, *n*) and tests whether the GW data independently reproduce the locked lab-scale convention TEP parameters.

The lab-fixed comparison gives Δ&chi;&sup2; = +0.031 and ΔBIC = +0.031 relative to ΛCDM, which is statistically inconclusive (|ΔBIC| < 2). A generic one-parameter linear redshift-bias model with *H*0 fixed to the ΛCDM best-fit (same k = 1 as lab-fixed TEP) also does not achieve decisive evidence, confirming that the current sample is too small to distinguish the specific TEP functional form from a generic redshift-dependent distance rescaling. The joint TEP fit incurs a BIC penalty of 2 ln(*N* = 15) = 5.42 chi2 units for its two extra parameters and does not overcome this penalty (ΔBIC = -5.01). A host-marginalized dark-siren diagnostic, summing over all plausible merged galaxy candidates with skymap sky probability as the host prior, gives Δ(-2 ln *L*) = +0.037, weakly favoring lab-fixed TEP. The evidence therefore supports only a weak model-sensitivity signal, not a robust model-selection detection or a Hubble-tension resolution.

![Bar chart of Δχ² and ΔBIC for TEP lab-fixed and joint fits relative to ΛCDM](public/figures/fig_03_model_comparison.png)

Figure 3. Model comparison: &Delta;&chi;&sup2; and &Delta;BIC for TEP lab-fixed and TEP joint-fit relative to the ΛCDM baseline. Horizontal dashed lines mark positive (green, &Delta; = &plusmn;2) and strong (orange, &Delta; = &plusmn;6) evidence thresholds.

## 5.4 Robustness Diagnostics

The strongest support for the lab-fixed TEP interpretation comes from fit stability rather than from the small absolute Δ&chi;&sup2;. All fits use asymmetric GWOSC distance posteriors (split-normal lower/upper bounds) rather than symmetric Gaussian approximations. Bootstrap resampling on the candidate-quality subset gives TEP lower in &chi;&sup2; in 90.6% of resamples, but closer to Planck in only 0.5% of resamples. Leave-one-out tests are stable by fit statistic: TEP is lower in &chi;&sup2; in 100% of event omissions, but closer to Planck in 0%. The redshift split gives small positive Δ&chi;&sup2; in both halves (low *z* = +0.019, high *z* = +0.006), so the present sample does not yet resolve a decisive growth of |*A*(*z*) &minus; 1| with redshift. A synthetic injection test confirms that the current sample size is below the threshold needed for decisive recovery of the locked lab-scale TEP amplitude; as the event count grows, the false-positive rate under the null and the recovery rate for the true signal provide calibrated sensitivity limits.

![Four-panel robustness diagnostic showing Planck alignment gain, resampling support, redshift split, and fixed-reference H0 residuals for lab-fixed TEP](public/figures/fig_06_tep_robustness.png)

Figure 6. Robustness diagnostics for the lab-fixed TEP signal. Positive Δ&chi;&sup2; values favor lab-fixed TEP. The diagnostics show a stable tiny fit-statistic preference, while the Hubble-scale shift moves away from external baselines under the current sign convention.

## 5.5 Adversarial Controls

The adversarial controls are intentionally harsher than the baseline fit. The locked sign of *φ*0 gives a tiny fit improvement (Δ&chi;&sup2; = +0.031) but moves the inferred scale farther from Planck by 1.2 km/s/Mpc; the wrong sign gives the opposite Hubble-scale direction but does not improve the fit. Zero coupling returns exactly to ΛCDM. A generic one-parameter linear redshift-bias model with *H*0 fixed to the ΛCDM best-fit (same degrees of freedom as lab-fixed TEP) does not outperform TEP by likelihood. Redshift shuffling destroys the event-distance pairing, while ΛCDM mock catalogs can reproduce a Δ&chi;&sup2; at least as large as observed with *p* = 0.154 and the observed Planck-alignment statistic with *p* = 0.611. Chronological splitting gives small positive Δ&chi;&sup2; in both halves (early +0.003, late +0.031), so the current sample does not yet reach discovery-level significance.

![Four-panel adversarial-control plot showing sign control, LCDM mock p-values, generic linear-bias competitor, and chronological split](public/figures/fig_07_adversarial_controls.png)

Figure 7. Adversarial controls. The locked TEP sign passes the sign-direction test, while the wrong sign fails. Mock-calibrated p-values and chronological splitting show that the present sample is suggestive but not decisive.

## 5.6 Host-Prior Ablation

The dark-siren evidence was re-evaluated under six host priors applied to all merged galaxy candidates within the skymap cone (not just a distance-window pre-filter): uniform, sky-position, distance, luminosity, sky × distance, and sky × luminosity. Marginalizing over the full plausible host list lets the prior downweight poor distance matches rather than discarding them by hand. Across these priors, the lab-fixed TEP model consistently raises the best-fit *H*0 by about 0.7–1.4 km/s/Mpc relative to ΛCDM, and all six priors weakly favor TEP by host-marginalized likelihood. The median Δ(-2 ln *L*) across priors is +0.023, so the robust conclusion is a weak likelihood preference paired with an upward Hubble-scale shift, not a Hubble-tension resolution.

![Two-panel host-prior ablation showing likelihood preference and H0 shift across host-prior choices](public/figures/fig_08_host_prior_ablation.png)

Figure 8. Host-prior ablation. The Hubble-scale shift is upward across plausible host-prior choices, while the likelihood preference remains weakly on the TEP side for all tested priors.

## 5.7 Conformal Scaling

The TEP conformal scaling factor *A*(*z*; *φ*0, *n*) quantifies the predicted deviation from GR propagation as a function of redshift. For the locked lab-scale convention parameters (*φ*0 = &minus;0.013, *n* = 1.0), *A*(*z*) departs from unity at the percent level by *z* &sim; 0.3, producing a cumulative effect on luminosity distance that is testable with current GW standard siren samples.

![Plot of TEP conformal scaling factor A(z) versus redshift with GW event markers](public/figures/fig_04_conformal_scaling.png)

Figure 4. TEP redshift-dependent conformal scaling *A*(*z*) for locked lab-scale convention parameters (*φ*0 = &minus;0.013, *n* = 1.0, red curve). Grey dashed line marks the GR limit *A* = 1. Blue points show the inferred *A*(*z*) for individual GW events.

## 5.8 Posterior Constraints

The joint MCMC fit to (*H*0, *φ*0, *n*) yields a 3-D posterior that tests whether the GW data independently reproduce the locked lab-scale convention TEP parameters. Red lines mark the locked lab-scale convention values. In the current sample the posterior is broad and the chain is diagnostic rather than converged enough for a precision parameter claim; the direct optimizer supplies the best-fit &chi;&sup2; used in model comparison.

![Corner plot of MCMC posterior samples for H0, phi0, and n from TEP joint fit](public/figures/fig_05_corner_posterior.png)

Figure 5. TEP joint-fit posterior *P*(*H*0, *φ*0, *n* | GW data) from emcee MCMC. Red lines mark locked lab-scale convention parameter values. Marginal distributions show 16th, 50th, and 84th percentiles.

## 6. Discussion

## 6.1 Implications for Cosmology

The corrected analysis separates two claims that were previously too easy to conflate. First, the locked lab-fixed TEP scaling gives a tiny and fairly stable fit-statistic improvement over ΛCDM without adding fitted parameters. Second, under the current *β* = -1 convention, that same scaling raises the candidate-quality best-fit *H*0 from 81.3 to 82.5 km/s/Mpc, moving the inferred scale farther from both SH0ES and Planck. The strongest current case for TEP is therefore a weak model-sensitivity signal, not a Hubble-tension reconciliation: the fit preference is stable under event omission, bootstrap resampling, and host-prior ablation, but the absolute Δχ² is tiny and the redshift split does not yet show a decisive growth signature.

## 6.2 Limitations

The present dark-siren implementation uses public skymaps and GLADE+ host candidates, with optional NED and local DESI DR1 subsets, so its redshift sample is still limited by galaxy-catalog completeness, localization area, cone truncation, and host ranking. The pipeline now reports both single-host and host-marginalized diagnostics; the latter shows that the current model-comparison evidence is weakly TEP-favoring once retained host candidates are summed over, but by far less than conventional evidence thresholds. ΛCDM mock calibration also shows that the observed Δ&chi;&sup2; is not rare under the null. The locked 2025 parameters preserve the out-of-sample structure of the lab-fixed TEP test, but definitive cosmological inference requires full posterior-sample treatment and deeper galaxy-catalog marginalization over all plausible host galaxies.

## 7. Conclusions

This paper implements a falsifiable standard-siren test of the locked 2025 TEP parameterization using combined public GWTC catalogs. The corrected pipeline shows that the lab-fixed conformal scaling gives a tiny, reproducible fit-statistic improvement over ΛCDM without adding fitted cosmological degrees of freedom, but under the current *β* = -1 convention it raises the candidate-quality inferred Hubble scale rather than reconciling it with external baselines. The robustness tests show that this weak fit preference is usually preserved under event omission, bootstrap resampling, and host-prior ablation, but the redshift split, ΛCDM mock calibration, and synthetic-injection sensitivity test do not yet support a detection claim. The reproducible pipeline records each analysis step, propagates event-level uncertainties, and exposes the remaining host-marginalization and catalog-completeness limitations that must be tightened in future releases.

## References

[1] Abbott, B. P., et al. (LIGO/Virgo). (2019). GWTC-1: A gravitational-wave transient catalog of compact binary mergers observed by LIGO and Virgo during the first and second observing runs. *Physical Review X*, 9(3), 031040.

[2] Abbott, B. P., et al. (LIGO/Virgo). (2021). GWTC-2: Compact binary coalescences observed by LIGO and Virgo during the first half of the third observing run. *Physical Review X*, 11(2), 021053.

[3] Riess, A. G., et al. (2022). A comprehensive measurement of the local value of the Hubble constant with 1 km/s/Mpc uncertainty from the Hubble Space Telescope and the SH0ES team. *The Astrophysical Journal Letters*, 934(1), L7.

[4] Planck Collaboration. (2020). Planck 2018 results. VI. Cosmological parameters. *Astronomy & Astrophysics*, 641, A6.

## Data Availability & Reproducibility

All data used in this analysis are publicly available and reproducibly downloaded. No synthetic, fabricated, or simulated data is used in the main analysis.

Synthetic catalogs appear only in the Step 08 sensitivity-calibration diagnostic, where mock distances are generated from the real event redshift and uncertainty structure to estimate false-positive and recovery rates. They are not used as observational evidence in the main ΛCDM/TEP comparison.

**Data sources:**

- GWOSC combined catalogs: GWTC-1-confident, GWTC-2, GWTC-2.1-confident, GWTC-3-confident, GWTC-4.0, GWTC-4.1, O4 Discovery Papers, GWTC-5.0 — gwosc.org

- GraceDB public skymaps (bayestar.fits.gz): gracedb.ligo.org

- GLADE+ Galaxy Catalog (VizieR VII/291): glade.plus

- NASA/IPAC Extragalactic Database redshift-bearing objects: ned.ipac.caltech.edu

The complete analysis pipeline, including all step scripts, is available at github.com/matthewsmawfield/TEP-LVK.