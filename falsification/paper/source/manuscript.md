---
title: "Mock controls show that host-dependent gravitational kernels reduce to the baryon fraction"
subtitle: "A falsification of the RGD kernel from galaxies to clusters"
author: "Karlos Maia — Independent Researcher"
date: ""
---

```{=latex}
\begin{center}\Large\bfseries Abstract\end{center}\vspace{-0.5em}
```

We subject the RGD (v4) effective gravitational kernel to non-circular falsification tests. None of the tests finds predictive content beyond the baryon bookkeeping expected in ΛCDM. The results that remain are methodological and exclusionary.

- **Circularity.** Mocks with no new physics, containing only the baryon fraction $f_{\rm bar}$(M), reproduce the RGD (v4) exponent law n(ξ) = ½ + 0.4068ξ: r = 0.978 and A = 0.411. Any kernel inversion that uses a total-mass host variable inherits this problem.
- **Baryonic host variable $\xi_b$.** For 23 groups (Sun et al. 2009) and 12 X-COP clusters, the law of this paper, with a purely baryonic $\xi_b$, predicts groups from clusters with an RMS of 0.084 dex under the stellar fraction of Gonzalez et al. (2013). The same gain appears in ΛCDM mocks with no new physics (this paper's law wins in 79% of replicas; p = 0.45). Under the stellar fraction of Budzynski et al. (2014), it loses (0.199 vs 0.076 dex).
- **Galaxies.** With M/L marginalized, λ = 0.93–0.95 rather than 0.999, and the RAR shape is preferred in test χ² (Δχ² = +166). Galaxies require A ≈ −9 ± 1, against A ≈ +2.5 in clusters: the single linear law is falsified (Δχ² = +73).
- **Exclusion of local laws.** At the same ($g_b$, $\alpha_b$), CLASH BCGs have $g_{\rm obs}$ 2.4 times higher than SPARC spirals (11.5σ). No local law in baryonic acceleration and potential depth describes both populations.
- **Theoretical bound.** The kernel admits a QUMOND-type action only for n ≤ 1. Explaining the BCG excess would require n ≈ 2.3–2.5.
- **H₀.** Replacing 70 with 67.4 km/s/Mpc changes the results by less than 2%.

**Conclusion.** In galaxies, this paper's kernel reproduces the RAR with one extra parameter and no gain (with M/L marginalized the RAR shape is preferred). In groups and clusters, it predicts nothing beyond what the baryon fraction already predicts, and in BCGs it fails quantitatively (0.27 dex predicted against 0.62 observed).

Code, data and logs: <https://github.com/Maia-KarlosMarden/rgd-kernel> (folder `falsification/`).

# Introduction

This work tests, and largely falsifies, the effective gravitational kernel of RGD (v4). The RGD (v4) kernel described disks, groups, hydrostatic clusters and lensing clusters with two global parameters and an exponent n(ξ) that depends on the host system.

**Context in the literature.** It is already known that the relation between observed and baryonic acceleration in clusters has a higher normalization than in galaxies (Chan & Del Popolo 2020; Tian et al. 2020; Pradyumna et al. 2021). In X-COP it has a complex shape that departs from the RAR of spirals (Eckert et al. 2022). For isolated galaxies, the early- versus late-type difference in the lensing RAR is disputed: Brouwer et al. (2021) find at least 6σ, while Mistele et al. (2024) find none under a strict isolation criterion.

**Contributions of this work:**

1. A mock control showing that exponent laws obtained by inverting a kernel with total mass reproduce the baryon fraction, not a gravitational law (Section 4).
2. A functional-form-independent test that excludes every local law $g_{\rm obs}$ = F($g_b$, $\Phi_b$) across BCGs and spirals (Section 1.1).
3. A purely baryonic host variable, $\xi_b$, tested against null models and ΛCDM mocks (Section 5).
4. The condition n ≤ 1 for the kernel to admit an action (Section 7).

## Exclusion of local laws

At the same pair ($g_b$, $\alpha_b$), where $\alpha_b$ is the internal photon deflection computed from baryons only, the 20 CLASH BCGs have $g_{\rm obs}$ higher than SPARC galaxies by +0.37 ± 0.03 dex (factor 2.4; 11.5σ). Including the potential depth of the cluster gas raises the difference to +0.47–0.49 dex. Any explanation needs a property of the host system, or unseen mass.

RGD (v4) used ξ = g₅₀₀/g†(z), but g₅₀₀ = GM₅₀₀/R₅₀₀² is $g_{\rm obs}$ itself at R₅₀₀. This paper asks what remains when the host variable is built without total mass. With these data, the answer is that nothing remains beyond what ΛCDM already predicts.

# The kernel and the cosmological scale

The RGD (v4) kernel is kept unchanged; what changes in this paper is the variable that sets the exponent n.

$$
g_{\mathrm{obs}} = K_{\mathrm{eff}}\, g_b, \qquad K_{\mathrm{eff}} = 1 + \lambda\, \frac{q}{(1-q)^{n}}, \qquad q = e^{-g_b/g^\dagger(z)}, \qquad g^\dagger(z) = \frac{c\,H(z)}{2\pi}
$$

In the deep regime, $K_{\rm eff}$ ≈ λ $(g^\dagger/$g_b$)^n$, and n = ½ yields flat rotation curves and the BTFR with a₀ = λ² g†.

**Adopted H₀.** RGD (v4) labelled H₀ = 70 km/s/Mpc as "Planck 2018"; the Planck value is 67.4. This paper adopts 67.4 and reports 70 as a sensitivity test. The two give g†₀ = 1.042 and 1.082 × 10⁻¹⁰ m/s².

| Quantity | H₀ = 70 | H₀ = 67.4 |
|:--|--:|--:|
| λ (fixed M/L) | 0.984 | 1.007 |
| λ (marginalized M/L) | 0.933 | 0.952 |
| $A_b$ (X-COP only) | 2.26 | 2.24 |
| This paper, LOO RMS, groups + X-COP (dex) | 0.074 | 0.074 |
| Predicted BCG log($g_{\rm obs}$/$g_b$) | 0.27 | 0.26 |

Only λ shifts noticeably, by about 2%, as expected: it compensates the change in g† to keep a₀ = λ²g†. No conclusion depends on the choice of H₀.

# Disk galaxies and the λ test

Item 5 of the RGD (v4) falsification agenda predicted that λ = 0.999 would hold with independent M/L and that the scatter would decrease. Without DESI BGS data, the test was run on SPARC with M/L marginalized per galaxy. The scatter decreases as predicted, but λ does not hold.

Method: 159 galaxies, 3253 points, χ² with $\sigma_{\rm int}$ = 0.057 dex and velocity uncertainties. Per-galaxy nuisance parameters are $\Upsilon_d$, $\Upsilon_b$, D and i, with 0.1 dex priors on Υ and Gaussian priors on D and i, as in Li et al. (2018). Validation is 5-fold cross-validation split by galaxy, with H₀ = 70.

| Model | M/L | Test χ²/N | Δχ² vs RAR | RMS/galaxy (dex) | λ across folds |
|:--|:--|--:|--:|--:|:--|
| RAR (free a₀) | fixed | 3.635 | 0 | 0.2131 | — |
| RGD, free λ | fixed | 3.902 | +869 | 0.2137 | 0.960–0.997 |
| RAR (free a₀) | marginalized | 1.010 | 0 | 0.1334 | — |
| RGD, free λ | marginalized | 1.061 | +166 | 0.1329 | 0.921–0.951 |
| RGD, λ = 0.999 | marginalized | 1.088 | +255 | 0.1343 | — |

Three conclusions:

- λ depends on the M/L treatment: 0.984 with fixed M/L and 0.933 ± 0.011 marginalized. Every fold falls below 0.999.
- With the same number of parameters, the RAR shape is preferred over the kernel shape. The RMS per galaxy is nearly equal, but χ² weights well-measured points, where the RAR fits the transition better.
- The marginalized value matches the "theoretical" λ of 0.937 in Table 1 of RGD (v4). This suggests declaring λ as derived from a₀ and g†, not as a parameter close to 1.

The BTFR remains a consequence of the kernel with n = ½. The RGD (v4) radial law n(w, ξ), however, broke it: with α(0.03) = 0.156, n varies across the disk and v² ∝ $R^{2n-1}$ is no longer flat. This paper drops that law.

# Circularity control: what the RGD (v4) law n(ξ) measures

The RGD (v4) law n(ξ) = ½ + 0.4068ξ is reproduced by simulations that contain no new physics. It encodes the baryon fraction of the systems, not a gravitational law.

In the mocks, 5000 systems have M₅₀₀ log-uniform between 2×10¹³ and 2×10¹⁵ M$_\odot$ and z between 0.03 and 0.7. They are assigned the baryon fraction of Sun et al. (2009), $f_{\rm bar}$ = 1.09 × 0.0708 $(M/10^{13})^{0.135}$, with 0.05 dex scatter. The kernel is then inverted with $K_{\rm obs}$ = 1/$f_{\rm bar}$, exactly as in RGD (v4).

| Test | r($n_{\rm emp}$, ξ) | Fitted A |
|:--|--:|--:|
| Mock with $f_{\rm bar}$(M) of Sun et al. | 0.978 | 0.411 |
| Mock with constant $f_{\rm bar}$ = 0.13 | 0.983 | 0.397 |
| RGD (v4), real data | 0.97 | 0.4068 |

The correlation arises from the inversion itself. With $K_{\rm obs}$ ≈ 7 nearly fixed and x = $f_{\rm bar}$ ξ, n must grow with ξ to keep K constant. In the deep regime:

$$
n \approx \frac{\ln\!\left(1/f_{\mathrm{bar}} - 1\right)}{\ln\!\left(1/(f_{\mathrm{bar}}\,\xi)\right)}
$$

The RGD (v4) control ($r_{\rm sim}$ = 0.31) used NFW mocks without a realistic $f_{\rm bar}$ and therefore missed this effect. Three consequences for the RGD (v4) results:

- The LOO RMS of 0.006 dex on CLASH does not measure predictive power: the predictor ξ contains the observable g₅₀₀.
- The r = 0.999 correlation of the residual with log(g₅₀₀/g†₀) is largely algebraic, for the same reason.
- The evolution of g†(z) is degenerate with $f_{\rm bar}$(M, z) as long as mass and redshift are correlated in the sample. This paper treats it as a hypothesis, not a detection.

# Baryonic host variable $\xi_b$

With a purely baryonic $\xi_b$, this paper's law does not predict groups better than ΛCDM mocks with no new physics lead one to expect. The apparent gain over the null models is produced by the kernel algebra combined with the $f_{\rm bar}$(M) trend.

**Definition.** $R_b$ is the radius where the mean baryonic density equals 500 $f_{\rm cos}$ $\rho_c$(z), with $f_{\rm cos}$ = $\Omega_b/\Omega_m$ = 0.157:

$$
\frac{M_b(<R_b)}{\tfrac{4\pi}{3}R_b^3} = 500\, f_{\mathrm{cos}}\, \rho_c(z), \qquad \xi_b = \frac{G\,M_b(<R_b)}{R_b^{2}\, g^\dagger(z)} \propto M_b^{1/3}, \qquad n = \tfrac12 + A_b\, \xi_b
$$

$\xi_b$ uses no total mass and is defined in the same way for disks and clusters.

**Sample and baryon accounting.**

- **Groups:** the 23 Tier 1 and 2 groups of Sun et al. (2009), Table 3, with measured M₅₀₀, r₅₀₀, r₂₅₀₀ and $f_{\rm gas}$ at r₅₀₀ and r₂₅₀₀. The table was checked row by row, which corrected an earlier version with 13 groups and an error in the r₅₀₀ of 3C 442A (495 kpc, not 405).
- **Clusters:** 12 X-COP clusters, with $M_{\rm HSE,500}$ and $f_{\rm gas}$ at R₅₀₀ and R₂₀₀ from Eckert et al. (2019).
- **Baryons:** $M_b$ = $M_{\rm gas}$ + $f_*$(M₅₀₀) $M_{\rm tot}$(<r), with the same rule for groups and clusters and stars tracing total mass. Two scenarios bracket the literature: G13, $f_*$ = 0.032 $(M_{500}/10^{14})^{-0.48}$ (Gonzalez et al. 2013); and B14, $f_*$ ≈ 0.013, nearly constant (Budzynski et al. 2014).
- **Target:** $g_{\rm obs}$($R_b$), from the hydrostatic mass. The local slope of M(<r) comes from the measured r₂₅₀₀ in groups and from the typical NFW ratio in clusters.

**Mocks.** 200 ΛCDM replicas with the same masses, redshifts and gas profiles as the real samples, $f_{\rm gas}$ = 0.088 $(M_{500}/10^{14})^{0.26}$ with 0.05 dex scatter, and $g_{\rm obs}$ given by the total mass, with no new physics.

| Model | k | LOO, G13 | X→G, G13 | Mocks, G13 | X→G, B14 | Mocks, B14 |
|:--|--:|--:|--:|--:|--:|--:|
| This paper: n = ½ + $A_b$ $\xi_b$ | 1 | 0.074 | 0.084 | 0.060 | 0.199 | 0.188 |
| Null 1: $M_{\rm tot}$ = $M_b$/$f_u$ | 1 | 0.094 | 0.146 | 0.135 | 0.321 | 0.318 |
| Null 2: $f_{\rm bar}$ ∝ $M_b^\beta$ | 2 | 0.074 | 0.121 | 0.091 | 0.076 | 0.062 |

RMS in dex. k: number of parameters. LOO: leave-one-out over all 35 systems. X→G: calibrated on X-COP only, predicting the 23 groups. Values are for H₀ = 67.4; with 70 nothing changes beyond the third decimal. Mock values are medians over the 200 replicas.

**Reading:**

- Under G13, this paper's law beats the nulls in the cross-prediction by 0.037 dex. In the ΛCDM mocks it wins in 79% of replicas, and an equal or larger gain occurs in 45% of them. The real gain is typical of a universe with no new physics.
- Under B14, it loses to null 2 in the data (0.199 vs 0.076 dex) and in the mocks (0.188 vs 0.062).
- In no scenario does the real data behave differently from the mocks. The $\xi_b$ law, like the RGD (v4) ξ law, reproduces the baryon bookkeeping, and $A_b$, from 2.5 to 3.7 depending on $f_*$, carries no gravitational information.

An earlier version of this manuscript reported the gain in groups as a positive result. The mock control reclassifies it as null.

# Predictions without refitting: galaxies and BCGs

The n($\xi_b$) law calibrated on clusters fails at both ends: it has the wrong sign in galaxies, and it falls far below the observations in BCGs.

**SPARC galaxies.** In galaxies, $\xi_b$ ranges from 0.0006 to 0.0105. With $A_b$ = 2.26, n ranges from 0.501 to 0.524, and the RMS per galaxy barely changes (0.2157 vs 0.2148 dex). The prediction is a slight curvature of the BTFR at the massive end: v² ∝ $R^{+0.05}$ in the deep regime. Leaving A free in galaxies, with χ² and nuisances:

| Configuration | H₀ = 70 | H₀ = 67.4 |
|:--|--:|--:|
| A measured in galaxies (mean of 5 folds) | −9.1 ± 1.2 | −8.2 ± 1.3 |
| Test Δχ², free A vs n = ½ | −157 | −130 |
| Test Δχ², A = 2.26 vs n = ½ | +73 | +67 |

Galaxies strongly prefer a negative slope, and the cluster value is disfavored. Part of the negative A may come from systematics, such as a mass trend in M/L or the shape difference between the kernel and the RAR. None of these readings makes $A_{\rm gal}$ compatible with +2.26. A single linear law in $\xi_b$ from galaxy to cluster scales is falsified.

**CLASH BCGs.** The sample is 20 BCGs with masses from Tian et al. (2020). The host has $\xi_b$ from 0.066 to 0.120, giving n = 0.65–0.77. This paper predicts log($g_{\rm obs}$/$g_b$) = 0.27 in the BCG region, against 0.62 observed (RMS 0.38 dex). Closing the gap would require n ≈ 2.3–2.5, and Section 7 shows that n > 1 is not admissible. The excess in cluster interiors remains unexplained in this paper.

# Action formulation

With fixed n, the kernel admits a non-relativistic QUMOND-type action, but only for n ≤ 1. The dependence on $\xi_b$ requires a non-local action that has not yet been written.

**QUMOND with fixed n.** Following Milgrom (2010), $\Phi_N$ is the Newtonian potential of the baryons and the physical potential Φ solves:

$$
\nabla^2\Phi_N = 4\pi G\rho_b, \qquad \nabla^2\Phi = \nabla\cdot\left[\nu\!\left(\frac{|\nabla\Phi_N|}{g^\dagger}\right)\nabla\Phi_N\right], \qquad \nu(y) = 1 + \lambda\,\frac{e^{-y}}{\left(1-e^{-y}\right)^{n}}
$$

The pair of equations follows from the Lagrangian density below, with Q′(z) = ν(√z). The action guarantees momentum and energy conservation for stationary systems.

$$
\mathcal{L} = -\frac{1}{8\pi G}\left[\,2\,\nabla\Phi\cdot\nabla\Phi_N - (g^\dagger)^2\, Q\!\left(\frac{|\nabla\Phi_N|^2}{(g^\dagger)^2}\right)\right] - \rho_b\,\Phi
$$

**Well-posedness condition.** The force law is well defined only if $g_{\rm obs}$ increases monotonically with $g_b$, that is, if y·ν(y) is increasing. The numerical check for 10⁻⁴ ≤ y ≤ 10³:

| n | min d ln $g_{\rm obs}$ / d ln $g_b$ | Status |
|--:|--:|:--|
| 0.50 | +0.505 | well posed |
| 0.80 | +0.200 | well posed |
| 1.00 | 0.000 | limit: $g_{\rm obs}$ → λ g†, constant |
| 1.20 | −0.200 | not invertible |
| 2.50 | −1.500 | not invertible |

For n > 1, $g_{\rm obs}$ **decreases** as $g_b$ increases in the deep regime, and no potential generates such a law. This excludes two elements of RGD (v4): the radial law n(w, ξ), which reaches n = 2.7–4.4 at BCG scales, and any attempt to close the BCG excess within this family of kernels.

**Non-locality.** $\xi_b$ is a functional of the mass distribution of the whole system, $M_b$ within $R_b$, not a function of fields at a point. There are two ways out:

- Promote $\xi_b$ to a local auxiliary field, such as the depth of $\Phi_N$. But local baryonic potential depth has already been tested and fails: BCGs and spirals differ at 11.5σ at the same ($g_b$, $\alpha_b$).
- Write an explicitly non-local action. It must define $\xi_b$ for embedded systems (a galaxy inside a group inside a cluster) and preserve momentum conservation between subsystems with different $\xi_b$.

**Relativity.** QUMOND is not relativistic. Using lensing masses, such as those from CLASH, requires a relativistic completion with $\gamma_{\rm PPN}$ = 1 in the weak field, such as that of Skordis & Złośnik (2021) for MOND. The RGD (v4) claim that the PMG-σ field reproduces $\gamma_{\rm PPN}$ = 1 could not be verified here and has been removed from the list of results.

# Discussion

None of the non-circular tests finds, in current data, any content in the RGD kernel beyond what ΛCDM baryon bookkeeping already explains.

**What this work establishes:**

1. Exponent laws obtained by inverting a kernel with a total-mass host variable are reproduced by mocks with no new physics (r = 0.978). The same holds for the baryonic version $\xi_b$: the gain in groups appears in 79% of ΛCDM mocks. Before any such kernel is interpreted, the mandatory test is a comparison with mocks that share the same $f_{\rm bar}$(M).
2. No local law in baryonic acceleration and potential depth describes BCGs and spirals at the same time (11.5σ).
3. Kernels of the form 1 + λ $q/(1-q)^n$ admit a QUMOND-type action only for n ≤ 1. This rules out explaining cluster interiors within this family.
4. In galaxies, the kernel with n = ½ is a viable interpolating function, but inferior to the RAR in test χ² (Δχ² = +166). Its λ depends on the M/L treatment.

**What this work does not support:**

- The RGD (v4) claim that dark matter is unnecessary in bound systems.
- A single law from galaxies to clusters: galaxies require A ≈ −9, clusters A ≈ +2.5.
- A detection of g†(z) ∝ H(z): it is degenerate with $f_{\rm bar}$(M, z) as long as mass and redshift are correlated in the sample.
- Any positive result for the $\xi_b$ variable.

**Relation to ΛCDM.** In groups and clusters, the data behave exactly like mocks in which $g_{\rm obs}$ comes from the total mass and $f_{\rm bar}$ follows the observed scaling relations. This does not prove ΛCDM, but it shows that the kernel does not distinguish itself from it.

**Relation to MOND.** In galaxies, the kernel is a choice of interpolating function with a₀ motivated by cH₀/2π. The MOND deficit in clusters remains: the host dependence that RGD introduced to resolve it does not survive the mock control.

**Why publish.** The value of this work lies in the controls, which apply to any modified-gravity program that uses cluster masses, and in documenting a hypothesis that was tested and failed.

# Falsification agenda with pre-specified criteria

Any future version of the program must beat ΛCDM mocks with the same $f_{\rm bar}$(M), not only empirical null models. The criteria below are fixed before the data, and a failure must be recorded without redefining the function a posteriori.

| Test | Data | Requirement | Falsified if |
|:--------|:------------|:------------|:------------|
| Gain over mocks | Any sample of groups or clusters | Cross-prediction gain above that expected in ΛCDM mocks | $P_{\rm mock}$(gain ≥ real) > 0.05 |
| Measured stellar mass | Groups with per-group photometry (e.g. eROSITA with optical catalogues) | Same result without choosing an $f_*$ scenario | Conclusion changes between G13 and B14 |
| Independence from $f_{\rm bar}$ | Groups and clusters with the same $M_b$ and different $f_{\rm bar}$ | n follows $\xi_b$, not $f_{\rm bar}$ | Partial correlation of n with $f_{\rm bar}$ at fixed $\xi_b$ ≥ 0.5 |
| BTFR curvature | BIG-SPARC or DESI, spectroscopic M/L | Outer slope of log v² vs log R of +0.05 ± 0.02 for $M_b$ > 10¹¹ M$_\odot$ | Slope ≤ 0 with error < 0.02 |
| Isolated ETGs vs BCGs | SLUGGS, ATLAS-3D, weak lensing (KiDS, Euclid); same $M_*$, strict isolation | Difference in log($g_{\rm obs}$/$g_b$) ≥ 0.15 dex | Difference < 0.05 dex |
| g†(z) at fixed mass | Clusters with the same $M_b$ at different z | g† ∝ H(z) | Better fit with constant g† (Δχ² > 4) |
| Non-local action | Theoretical | Action that generates n($\xi_b$), conserves momentum and handles embedded systems | No consistent formulation |

**Already decided in this work:**

- The single linear law between galaxies and clusters: falsified (Δχ² = +73).
- Explaining the BCG excess within the kernel family: excluded, because it requires n > 1.
- The gain of $\xi_b$ over the nulls in groups: indistinguishable from the mocks (p = 0.45).
- BTFR curvature: the first result on SPARC is unfavorable, but degenerate with M/L.

The isolated-ETG test should be read in light of the disagreement between Brouwer et al. (2021) and Mistele et al. (2024) about the isolation criterion.

# Reproducibility and data

**Code and data availability.** All code, input data, logs and results are public in the GitHub repository [Maia-KarlosMarden/rgd-kernel](https://github.com/Maia-KarlosMarden/rgd-kernel), folder `falsification/` (branch `main`). The command `bash falsification/code/run_all.sh` regenerates every logged result. The folder `v5/` of the same repository holds a separate, earlier line of work (a radial-term extension of the kernel) and is not the subject of this paper.

Every number in this manuscript comes from re-runnable Python scripts (NumPy, SciPy, scikit-learn), with validation always split by system, never by radial point.

| Script | Section | What it produces |
|:--|:--|:--------------|
| groups_test2.py ⟨H₀⟩ | 5 | 23 groups + 12 X-COP, G13 and B14 scenarios, LOO, cross-prediction and 200 ΛCDM mocks |
| rgd_tests.py ⟨H₀⟩ | 3, 6, 7 | λ with fixed and marginalized M/L; A in galaxies; kernel monotonicity |
| rgd_v5.py ⟨H₀⟩ | 4, 6 | RGD (v4) circularity mocks; predictions for SPARC and BCGs |
| form_test.py | 1.1 | BCGs vs SPARC at the same ($g_b$, $\alpha_b$) |
| audit3.py | 3 | χ² likelihood with per-galaxy nuisances |

**Checks performed:**

- Table 3 of Sun et al. (2009) checked row by row against the arXiv HTML version: 23 Tier 1 and 2 groups with measured M₅₀₀ and $f_{\rm gas,500}$.
- Stellar fraction: the published relations of Gonzalez et al. (2013) and Budzynski et al. (2014) replace the typical 3% value used in an earlier version.
- The ΛCDM mocks were applied to this paper's law as well, not only to RGD (v4).

**Known limitations:**

- Stellar masses come from scaling relations, not per-group measurements. The two relations used lead to different conclusions about this paper's law versus null 2, but to the same conclusion about this paper's law versus the mocks.
- Stars were distributed like the total mass. A more concentrated distribution would change $R_b$ in groups.
- X-COP redshifts are catalogue values. Masses were not rescaled with h; the effect is a few per cent.
- In the BCG test, M₅₀₀ comes from kT via the relation of Arnaud et al. (2005), but it enters only as $M_b^{1/3}$.
- The λ test uses SPARC with marginalized M/L, not DESI spectroscopic M/L.

# References {-}

::: {.refs}
Arnaud, M., Pointecouteau, E. & Pratt, G. W. 2005, The structural and scaling properties of nearby galaxy clusters. II. The M–T relation, *Astronomy & Astrophysics*, 441, 893. doi:[10.1051/0004-6361:20052856](https://doi.org/10.1051/0004-6361:20052856)

Brouwer, M. M., Oman, K. A., Valentijn, E. A., et al. 2021, The weak lensing radial acceleration relation: Constraining modified gravity and cold dark matter theories with KiDS-1000, *Astronomy & Astrophysics*, 650, A113. doi:[10.1051/0004-6361/202040108](https://doi.org/10.1051/0004-6361/202040108); arXiv:2106.11677

Budzynski, J. M., Koposov, S. E., McCarthy, I. G., et al. 2014, The similarity of the stellar mass fractions of galaxy groups and clusters, *Monthly Notices of the Royal Astronomical Society*, 437, 1362. doi:[10.1093/mnras/stt1965](https://doi.org/10.1093/mnras/stt1965)

Chan, M. H. & Del Popolo, A. 2020, The radial acceleration relation in galaxy clusters, *Monthly Notices of the Royal Astronomical Society*, 492, 5865. doi:[10.1093/mnras/staa225](https://doi.org/10.1093/mnras/staa225); arXiv:2001.06141

Eckert, D., Ghirardini, V., Ettori, S., et al. 2019, Non-thermal pressure support in X-COP galaxy clusters, *Astronomy & Astrophysics*, 621, A40. doi:[10.1051/0004-6361/201833324](https://doi.org/10.1051/0004-6361/201833324); arXiv:1805.00034

Eckert, D., Ettori, S., Pointecouteau, E., van der Burg, R. F. J. & Loubser, S. I. 2022, The gravitational field of X-COP galaxy clusters, *Astronomy & Astrophysics*, 662, A123. doi:[10.1051/0004-6361/202142507](https://doi.org/10.1051/0004-6361/202142507); arXiv:2205.01110

Ettori, S., Ghirardini, V., Eckert, D., et al. 2019, Hydrostatic mass profiles in X-COP galaxy clusters, *Astronomy & Astrophysics*, 621, A39. doi:[10.1051/0004-6361/201833323](https://doi.org/10.1051/0004-6361/201833323)

Gonzalez, A. H., Sivanandam, S., Zabludoff, A. I. & Zaritsky, D. 2013, Galaxy cluster baryon fractions revisited, *The Astrophysical Journal*, 778, 14. doi:[10.1088/0004-637X/778/1/14](https://doi.org/10.1088/0004-637X/778/1/14); arXiv:1309.3565

Lelli, F., McGaugh, S. S. & Schombert, J. M. 2016, SPARC: Mass models for 175 disk galaxies with Spitzer photometry and accurate rotation curves, *The Astronomical Journal*, 152, 157. doi:[10.3847/0004-6256/152/6/157](https://doi.org/10.3847/0004-6256/152/6/157); arXiv:1606.09251

Li, P., Lelli, F., McGaugh, S. & Schombert, J. 2018, Fitting the radial acceleration relation to individual SPARC galaxies, *Astronomy & Astrophysics*, 615, A3. doi:[10.1051/0004-6361/201732547](https://doi.org/10.1051/0004-6361/201732547)

Maia, K. M. 2026, A Cosmological-Background Effective Gravitational Kernel: Unifying disk galaxies, galaxy groups, galaxy clusters, and gravitational lensing under a single effective law (RGD v4), Zenodo preprint, version 1, 29 April 2026. doi:[10.5281/zenodo.19868177](https://doi.org/10.5281/zenodo.19868177); code: [Maia-KarlosMarden/rgd-kernel](https://github.com/Maia-KarlosMarden/rgd-kernel)

McGaugh, S. S., Lelli, F. & Schombert, J. M. 2016, Radial acceleration relation in rotationally supported galaxies, *Physical Review Letters*, 117, 201101. doi:[10.1103/PhysRevLett.117.201101](https://doi.org/10.1103/PhysRevLett.117.201101)

Milgrom, M. 2010, Quasi-linear formulation of MOND, *Monthly Notices of the Royal Astronomical Society*, 403, 886. doi:[10.1111/j.1365-2966.2009.16184.x](https://doi.org/10.1111/j.1365-2966.2009.16184.x); arXiv:0911.5464

Mistele, T., McGaugh, S., Lelli, F., Schombert, J. & Li, P. 2024, Radial acceleration relation of galaxies with joint kinematic and weak-lensing data, *Journal of Cosmology and Astroparticle Physics*, 2024(04), 020. doi:[10.1088/1475-7516/2024/04/020](https://doi.org/10.1088/1475-7516/2024/04/020); arXiv:2310.15248

Planck Collaboration (Aghanim, N., et al.) 2020, Planck 2018 results. VI. Cosmological parameters, *Astronomy & Astrophysics*, 641, A6. doi:[10.1051/0004-6361/201833910](https://doi.org/10.1051/0004-6361/201833910); arXiv:1807.06209

Postman, M., Coe, D., Benítez, N., et al. 2012, The Cluster Lensing And Supernova survey with Hubble (CLASH): An overview, *The Astrophysical Journal Supplement Series*, 199, 25. doi:[10.1088/0067-0049/199/2/25](https://doi.org/10.1088/0067-0049/199/2/25); arXiv:1106.3328

Pradyumna, S., Gupta, S., Seeram, S. & Desai, S. 2021, Yet another test of radial acceleration relation for galaxy clusters, *Physics of the Dark Universe*, 100765. doi:[10.1016/j.dark.2020.100765](https://doi.org/10.1016/j.dark.2020.100765); arXiv:2011.06421

Skordis, C. & Złośnik, T. 2021, New relativistic theory for modified Newtonian dynamics, *Physical Review Letters*, 127, 161302. doi:[10.1103/PhysRevLett.127.161302](https://doi.org/10.1103/PhysRevLett.127.161302); arXiv:2007.00082

Sun, M., Voit, G. M., Donahue, M., Jones, C., Forman, W. & Vikhlinin, A. 2009, Chandra studies of the X-ray gas properties of galaxy groups, *The Astrophysical Journal*, 693, 1142. arXiv:0805.2320

Tian, Y., Umetsu, K., Ko, C.-M., Donahue, M. & Chiu, I-N. 2020, The radial acceleration relation in CLASH galaxy clusters, *The Astrophysical Journal*, 896, 70. doi:[10.3847/1538-4357/ab8e3d](https://doi.org/10.3847/1538-4357/ab8e3d); arXiv:2001.08340

Umetsu, K., Zitrin, A., Gruen, D., et al. 2016, CLASH: Joint analysis of strong-lensing, weak-lensing shear and magnification data for 20 galaxy clusters, *The Astrophysical Journal*, 821, 116. doi:[10.3847/0004-637X/821/2/116](https://doi.org/10.3847/0004-637X/821/2/116); arXiv:1507.04385
:::

*Author lists are abbreviated with "et al." where long.*
