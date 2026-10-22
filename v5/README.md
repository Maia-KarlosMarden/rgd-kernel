# RGD v5 — extended tests and radial-term update

Scripts and transcribed tables behind the second paper
*"Relational Gravitational Dynamics: um kernel dependente de H(z) para sistemas ligados e seus limites cosmológicos"*.
Everything is plain Python (numpy, scipy). Run scripts from inside this folder
(some use `exec(open('other.py'))`, so the working directory matters; do not use `python -I`).

## What changed relative to v4

| | v4 | v5 |
|---|---|---|
| Radial term | exponent modulation `n = A(ξ) w^(-α(ξ))`, 2 extra parameters | scale shift `s = [Δ(r)/Δ500]^(p·Θ(ξ-ξc))`, `Δ/Δ500 = ρ̄(<r)/(500 ρc)`, 1 extra parameter (p ≈ 0.45–0.5) |
| X-COP 39 points (13 clusters × 0.5/1.0/1.5 Mpc) | 0.037 dex (reproduced; paper: 0.040) | 0.016–0.021 dex |
| Vikhlinin r2500 (w = 0.39–0.48) | n > 1, kernel diverges (1 of 10 no solution) | obs/pred 0.84, f_bar(r2500)/f_bar(r500) = 1.00 (obs 0.90) |
| Samples | SPARC, 23 groups, X-COP, CLASH | + 43 groups, Allen, Mantz, Vikhlinin, X-COP f_gas(r200), CLASH lensing masses, KMOS3D BTFR |

## Key results (see the paper for caveats)

- g† ∝ E(z)^α: clusters prefer α ≳ 1 (Allen χ² 69 at α=1 vs 126 at α=0; Mantz 117 vs 165); KMOS3D galaxies prefer α ≈ 0.5–0.75 (within ~10% systematics).
- Step switch in ξ (ξc ≈ 0.09–0.14) joins galaxies to groups/clusters; no data between 0.09 and 0.136.
- Tensions: f_bar–mass trend at r2500 (-0.31 per dex with p = 0.5), ~20% normalisation offset between Vikhlinin and X-COP, p(ξ) rejected (wrong sign r500→r200).
- Cosmology: the kernel applied to linear perturbations boosts gravity ~10²–10³×; a Δ-damper with plateau K_c = Ωm/Ωb reproduces ΛCDM linear growth by construction. No CMB claim.

## Files

| Topic | Scripts |
|---|---|
| Kernel, X-COP machinery, validation | `kernel_scenarios.py`, `kernel_v5_validation.py` |
| v4 vs v5 radial law | `v4_vs_v5.py` (X-COP 39 pts), `v4_indep.py` (Vikhlinin r2500, X-COP r200) |
| Clusters at r2500 / r500 | `vik_test.py`, `vik_test_p0.py`, `allen_test.py`, `mantz_z_test.py` (data: `allen_data.py`, `mantz_data.py`) |
| X-COP r500→r200 (Ettori+2019 Table 2) | `xcop200.py`, `xcop200_p0.py` |
| Groups | `groups_test.py` (data `groups_data.py`, `data/grupos_transcritos_Lovisari2015_Sun2009.csv`) |
| Stars, p(ξ), hydrostatic bias, LODO | `stars_test.py`, `p_xi_stars.py`, `hse_bias.py`, `lodo.py` |
| Lensing | `lens_test.py` (CLASH, Umetsu+2016) |
| g†(z) exponent and high-z BTFR | `alpha_scan.py`, `highz_tf.py` (KMOS3D, Übler+2017), `z_predict.py`, `switch_predictions.py` |
| ξ switch | `joint_switch.py`, `joint_switch_delta.py`, `joint_groups.py` |
| SPARC | `sparc_test.py`, `sparc_mbar.py`, `sparc_bridge.py` (need SPARC Table 2: set `SPARC_TABLE2=/path/to/SPARC_Table2.txt`, from astroweb.cwru.edu/SPARC) |
| ΛCDM comparison | `lcdm_compare.py` |
| Cosmological extension | `cmb_regime.py` (regime map and toy growth), `damper_test.py` |

## Data notes

Tables in `allen_data.py`, `mantz_data.py`, `groups_data.py` and the CSV were transcribed by hand from the published papers
(Allen+2008, Mantz+2014, Vikhlinin+2006 inside `vik_test.py`, Lovisari+2015, Sun+2009, Ettori+2019, Umetsu+2016) and
should be checked against the originals before reuse. Mantz f_gas are shell values, Allen and Vikhlinin cumulative.
The X-COP rows in `kernel_scenarios.py` use an f_gas scaling relation (the v4 "paper mode"), not measured profiles;
the measured-profile test is `xcop200.py`. Stellar masses follow Kravtsov+2018.

## Manuscript
`paper/` holds the English manuscript (PDF, LaTeX, plain text), the figures and `make_figs.py`.

## Referee-type checks (added with the final manuscript)
- `ref31.py` – cluster bootstrap, frozen-p predictions on Vikhlinin / Ettori r200 / groups / Allen, best-fit p per sample.
- `ref32.py` – sensitivity to hydrostatic bias and lensing calibration.
- `ref34.py`, `ref34b.py` – width and location of the xi switch (logistic m, xi_c, p) and the power-law alternative.
- `predict_frozen.py` – frozen prescription (p=0.5, xi_c=0.115, m=4, phi=0.7) to be applied to new samples before inspecting them.
