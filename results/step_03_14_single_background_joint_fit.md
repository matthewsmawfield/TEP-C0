# Step 03.14: Single-background joint fit (Pantheon+ × DESI BAO × CMB θ*)

Script: `scripts/steps/step_03_14_single_background_joint_fit.py`. Data: Pantheon+SH0ES (N = 1580; z_HD > 0.01, non-calibrators, full covariance, M profiled); official DESI DR2 and DR1 Gaussian BAO (`data/raw/desi_bao/`); Planck 2018 100θ* = 1.04110 ± 0.00031.

Every probe is evaluated on **one expansion history per model**. Only the ratio r_s*/r_d = 0.9819 is fixed, which TEP-HC preserves. The absolute ruler α = c/(H0 r_d) is free and shared by BAO and the CMB, so this is a pure test of the distance-redshift shape from z = 0.01 to z = 1090.

## Results (DESI DR2; DR1 in brackets)

| model | params | SN+BAO Δχ² vs ΛCDM | SN+BAO+CMB Δχ² | ΔAIC (all three) |
|---|---|---|---|---|
| ΛCDM | Ωm = 0.297 | 0 | 0 (χ² = 1400.9) | 0 |
| **Native (Paper 18 branch)** | Ωm = 0.296, ε_T = 0.019, z_T = 0.50 | −7.2 [−4.9] | **−8.1 [−6.3]** | **−4.1 [−2.3]** |
| M1 zero-Λ kernel, n = 2 | ε = 0.66, z_T = 18.7 (forced) | −5.0 [−4.5] | +112.7 [+65.2] | +114.7 [+67.2] |
| M1, n free | ε = 1.59, z_T = 3.5, n = 0.30 | −6.5 [−4.6] | +6.0 [+3.4] | +10.0 [+7.4] |

- **Corpus Step 03.05 point** (M1, ε = 0.883, z_T = 4.24) on its own background: with the BAO-fitted ruler it predicts 100θ* = 1.279, against 1.0411 observed. The comoving distance to last scattering comes out about 19% short. Step 03.05 passes the CMB only because it evaluates θ* on a separate hi_class ΛCDM-like background.
- **Native branch.** Profile in ε_T (Ωm, z_T, α re-optimised), Δχ²: 0 → −3.6 (0.005) → −6.2 (0.01) → −7.6 (0.015) → −8.0 (0.02) → −7.4 (0.025) → −5.7 (0.03) → +0.7 (0.04). The preferred transition is late, z_T ≈ 0.5. At the Paper 18 values (ε_T = 0.0043, z_T = 5) the fit is ΛCDM-equivalent, Δχ² = +0.8.
- **Native apparent w(z)** (with Ωm = 0.296): −0.95 (z = 0.1), −0.90 (0.5), −1.10 (1.0), ≈ −1.0 for z ≥ 1.5. The expansion modifier M(z) = H/H_ΛCDM has M → 1 for z ≳ 1.5, so the early universe and θ* are untouched.

## Conclusions
1. Only one TEP background is consistent with all three probes on a single history: the native H = M(z)·H_ΛCDM branch. It is modestly preferred over ΛCDM (Δχ² ≈ −8, about 2.4σ for 2 extra parameters).
2. The zero-Λ M1 kernel cannot fit SN, BAO and θ* together. Its SN-only and SN+BAO successes do not survive the CMB distance.
3. Recommended corpus change: Papers 26 and 27 should take their background from the native branch (or quote M1 as SN-only). The Paper 18 / TEP-HC fit should explore the late-transition region z_T ≈ 0.5 with the full Boltzmann likelihood.
