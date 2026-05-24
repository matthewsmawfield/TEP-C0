Yes. The fix is not mainly in the paper wording now — it is in the **analysis architecture**. The pipeline currently proves that several modules run and recover standard limits, but it does not yet force the decisive comparison:

> **Can Temporal Shear replace the scale factor, rather than merely perturbing it?**

Below is what needs to change.

---

# 1. Separate “validation gates” from “claim gates”

Right now the results summary treats everything as “research grade” and says “all claim gates open.” But most current results are **validation gates**, not replacement-claim gates. For example, BBN matching (\Lambda)CDM proves the Jakarta branch preserves standard light-element yields, not that expansion is reconstructed from Temporal Shear.  Likewise, the CMB zero-limit test proves the code recovers (\Lambda)CDM when TEP parameters vanish, but the file explicitly says the current CLASS build does not yet include TEP-specific source modifications. 

Change the pipeline output categories to:

| Gate type                 | Meaning                                                     | Current status |
| ------------------------- | ----------------------------------------------------------- | -------------- |
| **Infrastructure gate**   | Code runs; data ingest works; baseline cosmology reproduced | mostly passed  |
| **Zero-limit gate**       | TEP (\to\Lambda)CDM when TEP parameters vanish              | CMB passed     |
| **Preservation gate**     | BBN/CMB/growth not broken                                   | partly passed  |
| **Model-comparison gate** | TEP beats/matches (\Lambda)CDM under same likelihood        | not yet shown  |
| **Replacement gate**      | pure Temporal Shear works without primitive expansion       | not yet shown  |

The results summary should only say “pipeline validation gates passed” until the model-comparison and replacement gates are passed.

---

# 2. Add matched (\Lambda)CDM baselines for every result

The static metric result currently reports:

[
\chi^2/{\rm dof}=1.906
]

for 1,694 Pantheon+ supernovae.  But this number is not interpretable unless the exact same pipeline also reports the (\Lambda)CDM baseline.

Add this mandatory table to every output:

| Model                 | (\chi^2) | dof | (\chi^2/{\rm dof}) | AIC | BIC | (\Delta)BIC vs (\Lambda)CDM |
| --------------------- | -------: | --: | -----------------: | --: | --: | --------------------------: |
| (\Lambda)CDM baseline |          |     |                    |     |     |                           0 |
| Static metric         |          |     |              1.906 |     |     |                             |
| Mixed Temporal Shear  |          |     |                    |     |     |                             |
| Pure Temporal Shear   |          |     |                    |     |     |                             |

Without this, the analysis cannot claim model preference.

Also: (\chi^2/{\rm dof}=1.906) should be flagged as **not competitive by default** unless the matched (\Lambda)CDM result is similar or worse.

---

# 3. Add a pure Temporal Shear model

The current analysis appears to test a mixed form:

[
\ln(1+z_{\rm obs})
==================

(1-f_T)\ln(1+z_{\rm FLRW})
+
f_T\ln(1+z_T).
]

But a mixed model cannot prove that expansion is an emergent reconstruction. It only proves that a Temporal Shear component may improve or alter a standard expansion model.

You need three explicit models:

## M0: standard expansion

[
1+z = \frac{a_0}{a_{\rm em}}.
]

## M1: mixed expansion + shear

[
\ln(1+z_{\rm obs})
==================

(1-f_T)\ln(1+z_{\rm FLRW})
+
f_T\ln(1+z_T).
]

## M2: pure Temporal Shear reconstruction

[
\ln(1+z_{\rm obs})
==================

\int_{\gamma}
\Sigma_\parallel^{\rm eff} d\ell,
]

with

[
a_{\rm eff}
===========

\exp\left[
-\int_{\gamma}
\Sigma_\parallel^{\rm eff} d\ell
\right].
]

The Level 3/radical claim only starts to become viable if **M2 is competitive** with M0 across SNe + BAO + CMB acoustic anchors.

---

# 4. Fix the meaning of (f_T)

The current results include a Temporal Shear fraction:

[
f_T = 0.113 \pm 0.008.
]

If this is really a fraction, it says the effect is only about 11%, not a full reconstruction. That conflicts with the radical claim.

Fix this in one of two ways.

## Option A: rename it as an amplitude

Use:

[
\epsilon_T
]

instead of (f_T), and define it as a coupling amplitude:

[
\Sigma_\parallel^{\rm eff}
==========================

\epsilon_T , \Sigma_\parallel^{\rm model}
+
\mathcal C_T.
]

Then it no longer implies “only 11% of expansion.”

## Option B: make it redshift-dependent

Use:

[
f_T(z)
======

1-\exp[-(z/z_T)^n].
]

Then the model can say:

* local universe: near-(\Lambda)CDM / weak shear reconstruction;
* high redshift: Temporal Shear dominates;
* (f_T\to1) is the replacement limit.

The pipeline should report:

[
f_T(0),\quad f_T(0.5),\quad f_T(1),\quad f_T(2),\quad f_T(1100).
]

Not just one number.

---

# 5. Diagnose the supernova fit properly

The static model result is currently treated too positively. A (\chi^2/{\rm dof}=1.906) for Pantheon+ is not an “excellent” result without context. 

Add a diagnostic step:

## SN likelihood audit

Check whether the pipeline uses:

* full Pantheon+ covariance matrix;
* systematic covariance;
* intrinsic scatter treatment;
* nuisance parameters / absolute magnitude marginalization;
* peculiar velocity covariance;
* redshift covariance;
* consistent sample cuts;
* identical treatment for (\Lambda)CDM and TEP.

Then output:

```text
SN_AUDIT_PASS = true/false
uses_full_covariance = true/false
nuisance_marginalization = true/false
matched_LCDM_baseline = true/false
```

Until this passes, the static metric result should be labeled:

> proof-of-concept SN fit, not model preference.

---

# 6. Treat the 31% sound-horizon shift as a critical stress test

The TEP perturbation module reports:

[
r_s^{\rm TEP}=103.9 {\rm \ Mpc},
\qquad
r_s^{\Lambda{\rm CDM}}=150.9 {\rm \ Mpc},
]

a (-31.1%) shift. 

This is the most urgent thing to test. A 31% shift is huge. It can only work if the reconstructed angular-diameter distance changes consistently:

[
\theta_*=\frac{r_s}{D_A(z_*)}.
]

Add a dedicated acoustic consistency step:

```text
step_acoustic_consistency.py
```

It should compute:

[
r_s,\quad D_A(z_*),\quad \theta_*,\quad \ell_A,\quad R.
]

And compare:

| Quantity            | Planck/(\Lambda)CDM |       TEP | residual |
| ------------------- | ------------------: | --------: | -------: |
| (r_s)               |           150.9 Mpc | 103.9 Mpc |   -31.1% |
| (D_A(z_*))          |                     |           |          |
| (\theta_*)          |                     |           |          |
| (\ell_A)            |                     |           |          |
| shift parameter (R) |                     |           |          |

If (\theta_*) fails, the model fails CMB-I. If it passes, this becomes a major result.

---

# 7. Add BAO ratio tests immediately

Because (r_s) shifts by 31%, BAO cannot be treated casually. BAO surveys constrain ratios, not just distances:

[
D_M(z)/r_d,
\qquad
D_H(z)/r_d,
\qquad
D_V(z)/r_d.
]

Add:

```text
step_bao_ratio_consistency.py
```

For each BAO point, compute:

[
\Delta_i
========

\frac{O_i^{\rm TEP}-O_i^{\rm obs}}{\sigma_i}.
]

Output:

| BAO observable | redshift | observed | (\Lambda)CDM | TEP | residual |
| -------------- | -------: | -------: | -----------: | --: | -------: |

This is mandatory before claiming the sound-horizon shift helps.

---

# 8. Add supernova time-dilation and Tolman gates

A non-expansion model must not look like tired light. The same shear field that fits redshift must also predict time dilation:

[
\Delta t_{\rm obs}
==================

(1+z_T)\Delta t_{\rm em}.
]

Add:

```text
step_sn_time_dilation.py
```

Then Tolman dimming:

[
S_{\rm obs}
===========

S_{\rm em}
(1+z_T)^{-4}
\Xi_T(z,\mathcal E).
]

Add:

```text
step_tolman_surface_brightness.py
```

The pipeline should output:

| Test             | Pass criterion                                   |
| ---------------- | ------------------------------------------------ |
| SN time dilation | observed stretch compatible with (1+z_T)         |
| Tolman dimming   | (\Xi_T\approx1), or predicted deviation detected |
| Distance duality | (\eta(z)) residual matches TEP prediction        |

These are not optional. They are the first things a critic will ask.

---

# 9. Make distance duality the flagship discriminator

The current results mention a predicted distance-duality deviation:

[
\Delta\eta \approx 3.5% \quad \text{at } z=2.
]

That should become a formal output, not prose.

Add:

```text
step_distance_duality.py
```

Compute:

[
\eta(z)
=======

\frac{D_L}{D_A(1+z)^2}.
]

Then fit:

[
\eta(z)
=======

1+\eta_1 z+\eta_2 z^2
]

and compare with the TEP prediction:

[
\eta_{\rm TEP}(z,\mathcal E)
============================

1+\Delta\eta_T(z,\mathcal E).
]

Output:

| Model        | (\eta(z=1)) | (\eta(z=2)) | compatible with data? |
| ------------ | ----------: | ----------: | --------------------- |
| (\Lambda)CDM |           1 |           1 |                       |
| TEP          |             |     (1.035) |                       |

This is a clean, publishable discriminator.

---

# 10. Add null-injection tests

Before any detection claim, prove the pipeline does not hallucinate TEP.

Create mock data:

## Mock A: (\Lambda)CDM universe

Generate synthetic SNe + BAO + CMB anchors from (\Lambda)CDM.

Expected recovery:

[
\epsilon_T \approx 0
]

or:

[
f_T \approx 0.
]

## Mock B: injected Temporal Shear universe

Generate synthetic data with known (\Sigma_\parallel^{\rm eff}).

Expected recovery:

[
\hat{\epsilon}_T \approx \epsilon_T^{\rm injected}.
]

Output:

| Mock              | injected TEP | recovered TEP | pass/fail |
| ----------------- | -----------: | ------------: | --------- |
| (\Lambda)CDM mock |            0 |            ~0 |           |
| TEP mock          |        known |     recovered |           |

This is essential. Otherwise, a flexible pipeline may detect “Temporal Shear” in anything.

---

# 11. Add held-out prediction

Split the data:

## Training

* Pantheon+ subset;
* older BAO;
* Planck compressed anchors.

## Held-out

* DES-SN / Union3 / DESI BAO if available in the pipeline;
* cosmic chronometers;
* strong-lens time delays;
* lensed SNe;
* standard sirens.

Report:

[
\chi^2_{\rm train},
\quad
\chi^2_{\rm test},
\quad
{\rm LOO-CV},
\quad
{\rm WAIC}.
]

A radical model that only fits the training data is not convincing. A radical model that predicts held-out residuals is serious.

---

# 12. Implement CMB in three stages, not one

Given the current file says the CLASS build does not yet contain TEP-specific modifications, the pipeline should explicitly stage CMB validation. 

## CMB-I: compressed acoustic geometry

Fit:

[
\theta_*,
\quad
\ell_A,
\quad
R.
]

## CMB-II: compressed peak structure

Fit:

[
\omega_b,\quad \omega_c,\quad n_s,\quad A_s,\quad \tau,
]

and peak-height proxies.

## CMB-III: full spectra

Modify CLASS/CAMB source equations and compute:

[
C_\ell^{TT},\quad C_\ell^{TE},\quad C_\ell^{EE},\quad C_\ell^{\phi\phi}.
]

Your current pipeline has passed CMB-zero-limit, but not CMB-III.

---

# 13. Add matter-frame BBN guardrail

BBN currently matches (\Lambda)CDM, which is good.  But if later TEP parameters change early-time evolution, the pipeline needs a guardrail.

Add:

```text
step_bbn_guardrail.py
```

Output:

[
Y_p,\quad {\rm D/H},\quad ^3{\rm He/H},\quad ^7{\rm Li/H}
]

and compare to observational tolerances.

Pass criterion:

[
|\Delta Y_p|<1\sigma\text{ or stated threshold},
\qquad
|\Delta{\rm D/H}|<1\sigma\text{ or stated threshold}.
]

Label this:

> Matter-frame thermal-history preservation.

---

# 14. Proposed revised step list

Your pipeline should look like this:

```text
00_config_and_priors
01_ingest_pantheon_plus
02_ingest_bao
03_ingest_planck_compressed
04_ingest_growth_data
05_ingest_bbn_constraints

10_lcdm_baseline_fit
11_static_metric_fit
12_mixed_temporal_shear_fit
13_pure_temporal_shear_fit

20_sn_likelihood_audit
21_sn_time_dilation_gate
22_tolman_surface_brightness_gate
23_distance_duality_gate

30_sound_horizon_calculation
31_acoustic_consistency_gate
32_bao_ratio_gate
33_cmb_zero_limit_gate
34_cmb_modified_class_gate

40_bbn_preservation_gate
41_growth_sigma8_gate
42_pk_and_fs8_gate

50_null_injection_lcdm
51_null_injection_tep
52_train_test_split
53_loo_cv_waic

60_evidence_matrix
61_claim_gate_classifier
```

The final claim-gate classifier should output something like:

```json
{
  "infrastructure_gates": "passed",
  "zero_limit_gates": "passed",
  "preservation_gates": "partial",
  "model_comparison_gates": "pending",
  "replacement_gate": "pending",
  "allowed_claim": "Temporal Shear Cosmology pipeline validated; replacement claim not yet established"
}
```

Later, if pure shear passes:

```json
{
  "allowed_claim": "Pure Temporal Shear reconstruction is competitive with primitive expansion across SNe+BAO+CMB-I"
}
```

---

# 15. What needs to change in the current “results” wording

Replace:

> “All claim gates OPEN. Zero blockers. Pipeline complete.” 

with:

```markdown
All foundational infrastructure gates are open. The pipeline now reproduces standard BBN yields, recovers \(\Lambda\)CDM in the CMB zero limit, validates baseline growth calculations, and produces a nontrivial static-distance fit to Pantheon+. The decisive claim gates now move to model comparison: pure Temporal Shear reconstruction must be tested against matched \(\Lambda\)CDM baselines across SNe, BAO ratios, CMB acoustic geometry, supernova time dilation, Tolman dimming, distance duality, BBN preservation, and growth.
```

That is the correct status.

---

# Final answer

To fix the analysis, do these five things first:

1. **Add matched (\Lambda)CDM baselines everywhere.**
2. **Add a pure Temporal Shear reconstruction model, not just mixed FLRW+TEP.**
3. **Run acoustic consistency and BAO ratio tests for the 31% sound-horizon shift.**
4. **Audit the Pantheon+ likelihood, because (\chi^2/{\rm dof}=1.91) is not “excellent.”**
5. **Add null-injection and held-out prediction tests so the pipeline cannot falsely detect Temporal Shear.**

Only after those pass should the results say the expansion-replacement claim is supported. Right now, the honest interpretation is:

> **The pipeline is technically viable and has passed several foundational validation gates. The decisive replacement analysis still needs matched model comparison, pure Temporal Shear fitting, and acoustic/BAO consistency tests.**
