# Freeze of RGD v6 (primary) and v5-calibrated (comparator) for blind tests — 2026-10-25

**Status:** parameters below are frozen. Do not change them after the tag `v6-freeze-20261025`. Any change is a new version with a new tag.

## Models (code: `v5/frozen_v6.py`, constants in `FROZEN`)
- Common kernel: K = 1 + λ e^(−x/s) / (1 − e^(−x/s))^n, λ = 0.999, n = 0.5 + 0.4068 ξ, g† = cH(z)/2π (H0 = 70, Ωm = 0.3), Θ(ξ) = [1 + (0.115/ξ)^4]^−1.
- **v6 (primary):** s = (g/g500)^(qΘ), **q = 1.26** (equivalently x_eff = x (ξ/y)^(qΘ)).
- **v5-calibrated (comparator):** s = (Δ/Δ500)^(pΘ), **p = 0.61**.
- Stars: Kravtsov+2018, φ = 0.7 inside r2500, 1.0 inside r500 and r200 (archived prescription).
- Amplitude systematics attached to every prediction: sample-level ln-offset ~ N(0, 0.20²). Intrinsic scatter τ (ln f): v6 0.151, v5 0.198.

## Calibration provenance
`v5/ref58.py` / `ref58_results.txt`: 331 points, 6 samples (Vikhlinin, Ettori, Allen, Groups-L, Groups-S, XXL-C1), per-sample offset marginalised with σδ = 0.20.
Best fits: p = 0.610 ± 0.015, q = 1.259 ± 0.026 (statistical; LOSO choices spanned p 0.57–0.62, q 1.20–1.30, so use p ± 0.03, q ± 0.05).
Not blind: the radial family (g/g500 vs Δ/Δ500) was selected after seeing these data (`ref52`, `ref52b`). The calibration data are therefore *not* eligible for the blind test.
The code reproduces the fit code to machine precision when H0 = 70; samples with H0 = 71 (Vikhlinin) or 73 (Sun) differ by 0.6–2%.

## Blind-test protocol (to apply BEFORE inspecting a new sample's gas fractions)
1. Compute f_bar,pred for every system with both `predict_v6` and `predict_v5` (inputs: M(<r), r, M500, R500, z, from the sample's own mass profile).
2. Observed baryon fraction = gas fraction + stars (frozen φ).
3. Primary statistic: marginal log-likelihood with one free amplitude offset for the sample (prior N(0, 0.20²)) and τ fixed at the values above; report ΔlnL = lnL(v6) − lnL(v5).
   - v6 is favoured if ΔlnL > +3, v5 if < −3, otherwise inconclusive.
4. Secondary (offset-free): for systems measured at two radii, ρ = ln[f_obs(r)/f_obs(r500)] − ln[f_pred(r)/f_pred(r500)]; report mean ± s.e.m. and rms. A model is *falsified* by |mean ρ| > 0.15 at > 3σ in a sample of ≥ 15 systems.
5. Report posterior offset of the new sample; |δ| > 0.40 means the sample lies outside the calibration range and is reported separately.
6. Eligible new data: HSC-XXL per-system gas at r500, eFEDS bins (tabulated), X-COP resolved profiles (no constant-f assumption), CLASH, any sample not listed in the provenance above.
