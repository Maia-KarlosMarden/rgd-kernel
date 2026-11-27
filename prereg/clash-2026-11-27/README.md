# Negative record: host baryonic mass prediction for CLASH clusters (2026-11-27)

**Status: H-host is not supported in this form.** The primary criterion is INCONCLUSIVE and the
LambdaCDM control reproduces the data far better, so there is no evidence beyond LambdaCDM.

## Hypothesis and freeze
H-host: the gravitational boost depends on the host system's baryonic mass, here M_gas(<500 kpc),
using the existing H5 functional form. Calibrated on SPARC and X-COP only. The prediction used CLASH
baryons only. Its hashes (`CONGELADO.sha256`, commit d7e21f0) were pushed before the lensing masses
of this test were opened. Criteria: `runs/clash_pre/PROTOCOLO.md` (Portuguese).

Not blind: CLASH BCG data (14-30 kpc) had been seen in an earlier run. Radii of 100-600 kpc had not
been used. The H5 form had already failed the calibration criteria before this test.

## Result (15 X-ray-selected clusters, median stellar scenario)
e = log10 K_predicted - log10 K_observed, K = g_obs/g_bar. Brackets: 95% cluster bootstrap.

| Model | RMS (dex) | Bias (dex) |
|---|---|---|
| LambdaCDM mock, f_gas = 0.12 | 0.081 [0.060, 0.103] | +0.015 |
| H-host | 0.173 [0.142, 0.201] | -0.122 |
| Local control, Phi = g_b r | 0.182 [0.155, 0.203] | -0.155 |
| RAR calibrated on SPARC+X-COP | 0.341 | -0.330 |
| RAR calibrated on SPARC only | 0.577 | -0.571 |

- **T1 (transport): INCONCLUSIVE** in all three stellar scenarios. RMS lies between the pass (0.15)
  and fail (0.20) thresholds; the gain over the local control is 0.009 dex (0.05 required).
- **T2 (radial shape, descriptive):** observed K is nearly flat (log10 K = 0.91-0.98); H-host rises
  from 0.57 to 0.97 between 100 and 600 kpc (error -0.28 dex at 100 kpc, +0.01 at 600 kpc).
- **T3 (mass trend at 300 kpc):** Theil-Sen slope 0.28, 95% CI [-0.04, 0.84]; not detected, weak test.
- **T4 (LambdaCDM):** mocks with f_gas = 0.10-0.14 give RMS 0.08-0.10 dex, against 0.17 for H-host.

## Limits on interpretation
- The host variable saturates: every X-COP and CLASH cluster has t >= 0.998. The test is therefore
  "does the cluster branch learned on X-COP transfer to CLASH", not a test of mass dependence.
- The data are NFW fits and the mocks use NFW halos, which favours the mocks in radial shape.
- In the 4 high-magnification (merging) clusters the order reverses (H-host 0.13, f_gas=0.12 mock
  0.20). This is outside the criteria, small, and not counted in favour of H-host.
- The SPARC-only RAR falls below lensing by a factor of about 3.7. This is a statement about that
  relation, not support for H-host.
- Deviations from the protocol and open checks: `runs/clash_eval/DESVIOS.md`. In short: primary
  sample 15, not 16; the Umetsu et al. table was machine-extracted (see the update in that file)
  and should be checked by a person against the paper; M200/c200 uncertainties were not propagated.

**Decision:** no new functional form will be fitted to CLASH, which has now been seen. A mass-continuous
form would need a non-local action written in advance and a new sample (preferably galaxy groups).

## Files and reproduction
- `runs/clash_pre/`: protocol, frozen models and predictions. Verify with
  `cd runs/clash_pre && sha256sum -c CONGELADO.sha256`.
- `clash/avaliar.py`: evaluation; run `python3 clash/avaliar.py` from this folder
  (numpy, pandas, scipy, astropy).
- `clash/congelar.py`, `clash/cobertura.py`: provenance only. `congelar.py` imports pipeline modules
  that are not included here, and `cobertura.py` needs the ACCEPT2.0 and X-COP archives.

## Data sources
- Gas profiles: ACCEPT2.0 (Donahue et al.), https://zenodo.org/records/17740132, CC BY 4.0;
  gas masses here are derived from those electron density profiles.
- Lensing NFW parameters: Umetsu et al. 2016, ApJ 821, 116 (arXiv:1507.04385), Table 2.
- Calibration (not redistributed here): SPARC; X-COP.

## Addendum
`ADDENDUM.md` summarizes this record together with the I1 replication (`runs/replica_I1/`).
