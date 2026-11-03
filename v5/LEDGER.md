# RGD improvement ledger

Every attempt to improve RGD gets one status. A status changes only by citing a new test and its result file.

## Statuses

- **Fruitful**: improves on data not used in the fit (leave-one-sample-out or blind), without being circular, and survives free per-sample amplitude offsets.
- **Fruitful (provisional)**: passes the above on calibration data only. The blind test in `FREEZE.md` is still owed.
- **Failed**: rejected by a test that could have passed it.
- **Inconclusive**: depends on data not available, or the only test is circular.
- **Invalid**: the evidence was circular and is withdrawn.

## Criteria (in order)

1. New data: does it improve on samples excluded from the fit?
2. Parameters: does it reduce free parameters, or justify them?
3. Circularity: could the test have failed?
4. Systematics: does it survive per-sample offsets (prior sd 0.20 in ln f)?

## Ledger

| # | Attempt | Status | Evidence | Files |
|---|---|---|---|---|
| 1 | Virial construction Delta_RGD ~ 300 and two other earlier approaches | Failed | Falsified in the M_X00 bootstrap (earlier program records) | n/a |
| 2 | v5 with p = 0.5 fixed ("parameter-free") | Failed as pure prediction | Fit prefers p = 0.610 +/- 0.015; valid as a calibration only | `ref58_results.txt`, `frozen_v6.py` |
| 3 | Empirical equation for p from covariates (254 systems) | Failed | No validated equation; best covariates (z, instrument) are sample-driven | `ref51.py`, `ref51_results.txt` |
| 4 | p as a dimensionless formula | Failed | Same outcome as 3 | `ref52.py`, `ref52b.py`, `ref52_results.txt` |
| 5 | v6 radial law, s = (g/g500)^(q Theta), q = 1.26 | Fruitful (provisional) | dlnL = +47.8, dAIC = -95.5; LOSO +24.5 vs -30.5; tau 0.151 vs 0.198; better or tied in 6/6 samples. Family chosen after seeing the data | `ref56_results.txt`, `ref57_results.txt`, `ref58_results.txt` |
| 6 | v6 with q = 1 ("threshold follows local dynamical frequency") | Failed as identity | q = 1 rejected, dlnL = 43; survives only as approximation | `ref58_results.txt` |
| 7 | Marginalizing per-sample offsets (Sherman-Morrison) | Fruitful | p, q stable for prior sd 0.1 to 0.3 | `ref58.py` |
| 8 | Forward selection model | Partially fruitful | Explains part of the sample disagreements, not all | `blind_selection_model.py`, `blind_selection_results.txt` |
| 9 | X-COP stress test with q = 1 (rms 0.002 dex) | Invalid | Circular: constant f_bar input per radius. Not to be used as evidence | `ref53.py` |
| 10 | eFEDS tension | Inconclusive | Bins read by eye from a figure; eFEDS also disagrees with HSC-XXL | `blind_xxl_efeds_results.txt` |
| 11 | Offsets between samples and XMM vs Chandra (up to +/-0.3 in ln f) | Inconclusive | Pattern present, cause unidentified (calibration, stars, selection) | `ref57_results.txt`, `ref58_results.txt` |

## Status triggers

- Row 5 becomes **Fruitful** if the `FREEZE.md` blind test gives dlnL > +3 for v6 over v5 on eligible data.
- Row 5 becomes **Failed** if v5 wins by more than 3, or if |mean ratio| > 0.15 at more than 3 sigma in a sample with at least 15 systems.
- Row 9 would be replaced by a valid test if the resolved X-COP profiles are ingested with f_bar measured per radius.
- Rows 10 and 11 need tabulated eFEDS bins and per-system HSC-XXL gas at r500.

## Scope note

All gains in rows 5 to 8 concern the empirical gas law of clusters. None of them tests the galactic side (SPARC, lambda, g-dagger), so "fruitful" here means a better cluster fit, not a stronger derivation.

## Freeze reference

Commit `d9d5461` on `main`. The tag `v6-freeze-20261025` exists locally only; the push returned HTTP 403.
