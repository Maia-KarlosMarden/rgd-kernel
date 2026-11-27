# Addendum to the RGD falsification study: cluster profiles from 100 kpc outwards

Karlos Maia. Draft for author review, 2026-11-27.

This addendum extends the study in [`falsification/`](../../falsification/). It changes none of its
conclusions and reports no positive evidence for the RGD kernel or for physics beyond LambdaCDM.

## Data and conventions
- **X-COP** (hydrostatic masses, M_FORW), six clusters with measured stellar profiles: A85, A644,
  A1795, A2029, A2142, ZW1215.
- **CLASH** (lensing, spherical NFW fits of Umetsu et al. 2016, ApJ 821, 116), 15 X-ray-selected
  clusters; gas from ACCEPT2.0; stars bracketed with the X-COP stellar-to-gas ratios.
- K = g_obs/g_bar. g_dagger uses H0 = 67.4; cluster data are kept at H0 = 70.
- Using lensing masses as g_obs assumes that lensing and dynamical masses coincide
  (gamma_PPN = 1). The kernel has no relativistic completion that establishes this (paper, Sec. 7).

## 1. Relations calibrated on galaxies fall far below clusters
Bias = mean of log10 K_predicted - log10 K_observed.

| Sample | Radii (kpc) | RGD, n = 1/2, SPARC | RAR, SPARC |
|---|---|---|---|
| X-COP, hydrostatic | 100-1000 | -0.485 | -0.490 |
| CLASH, lensing | 100-600 | -0.570 | -0.571 |

Neither relation was refitted on clusters. This agrees with the long-known cluster deficit of
MOND-like relations and is not a new finding.

## 2. A host-baryonic-mass branch is not supported
H-host used the H5 transition form with M_gas(<500 kpc) as host variable, calibrated on SPARC and
X-COP; its CLASH predictions were hashed and pushed (commit d7e21f0) before the lensing masses of
that test were opened.

- Primary criterion: **inconclusive** (RMS 0.173 dex, bias -0.122; local control 0.182).
- LambdaCDM control: **much better** (RMS 0.081 dex at f_gas = 0.12; 0.08-0.10 over 0.10-0.14).
- Predicted K rises steeply with radius; observed K is nearly flat.
- The host variable saturates for all clusters, so mass dependence among clusters was not tested.
- The H5 form had already failed its calibration criteria before this test.

H-host is a phenomenological branch, not the paper's xi_b law, and has no action. See
[`README.md`](README.md) for the full record.

## 3. The departure of I1 from zero replicates but does not discriminate
R0 = ln(K-1) + x - ln(lambda) - L/2, with x = g_bar/g_dagger and L = -ln(1 - exp(-x));
I1 = dR0/dln r. The pure n = 1/2 kernel gives I1 = 0.

- Earlier X-COP run, 100-1000 kpc: median I1 = -0.31 [-0.41, -0.24]. SPARC galaxies in the same x
  window (79): +0.01. A radial hydrostatic bias rising from 0 to 30 % would reduce the cluster
  value to -0.14.
- Common radii, 100-600 kpc:

| Sample | Median I1 | 95 % CI | Clusters with I1 < 0 |
|---|---|---|---|
| CLASH, lensing | -0.31 | [-0.36, -0.23] | 15 / 15 |
| X-COP, hydrostatic | -0.24 | [-0.32, -0.09] | 6 / 6 |

- By the criterion written before this computation, the result **replicates** in all three stellar
  scenarios.
- It is **reproduced by LambdaCDM mocks**: the observed median lies inside their central 95 %
  interval for f_gas = 0.10, 0.12 and 0.14, in both samples. For CLASH at 0.12: mock -0.29,
  observed -0.31.
- In CLASH the baryonic terms of R0 supply almost all of the slope (median -0.39), while ln(K-1)
  contributes +0.05: K stays nearly flat while g_bar falls.

## Relation to the falsification agenda (paper, Sec. 9)
Nothing here satisfies the requirement of a gain over LambdaCDM mocks, and no non-local action is
offered. The mocks used here fix f_gas at three values with NFW halos; they are **not** the
f_bar(M) mocks that Sec. 9 asks for, so this addendum should not be read as a Sec. 9 test.

## Limitations
- Neither cluster analysis is blind. The X-COP and I1-replication protocols have no external
  timestamp; only the H-host hashes do.
- CLASH data are NFW fits and the mocks use NFW halos. This favours agreement between them and
  imposes part of the radial shape. Non-parametric masses would be needed to decide.
- Small samples (6 and 15); bootstrap intervals with 6 clusters are coarse.
- M200/c200 uncertainties were not propagated. The Umetsu table was machine-extracted twice with
  identical results; a check by a person against the paper is still pending.
- Satellite stellar profiles are missing for CLASH.

## Conclusion
No further functional form should be fitted to these samples. A continuation needs a non-local
action written in advance, mocks with f_bar(M), and an unseen sample.

## Provenance
`runs/replica_I1/` (criterion, results, per-cluster slopes) and `clash/replica_I1.py`. The script
needs pipeline modules and X-COP files that are not included here; it is provided for provenance.
The X-COP run behind Sec. 1 and the first bullet of Sec. 3 is in [`../xcop-2026-11-26/`](../xcop-2026-11-26/).
