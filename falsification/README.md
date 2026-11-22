# Falsification study — Mock controls show that host-dependent gravitational kernels reduce to the baryon fraction

*A falsification of the RGD kernel from galaxies to clusters.*
Karlos Maia — Independent Researcher

This folder holds the manuscript, code, input data and logged results of a falsification study of
the RGD (v4) kernel: **it does not confirm RGD (v4)**. It is independent of the `v5/` folder of this
repository, which holds a separate line of work (radial-term extension and the v6 freeze). The surviving results are
methodological and exclusionary (see the manuscript abstract).

## Headline results

| Result | Value |
|:--|:--|
| Circularity: mocks with only f_bar(M) reproduce the v4 law n(ξ)=½+0.4068ξ | r = 0.978, A = 0.411 |
| This paper (baryonic ξ_b) X-COP→groups RMS, f* = G13 | 0.084 dex vs nulls 0.146 / 0.121 |
| Same gain in ΛCDM mocks with no new physics | wins 79 % of replicas, p = 0.45 |
| This paper, f* = B14 | 0.199 dex vs nulls 0.321 / 0.076 (loses) |
| Galaxies vs clusters (single linear law) | A_gal ≈ −9 ± 1 vs A_cl ≈ +2.5; Δχ² = +73 |
| λ with M/L marginalized | 0.933 ± 0.011; RAR preferred, Δχ² = +166 |
| BCGs vs SPARC at same (g_b, α_b) | g_obs 2.4× higher (11.5σ) |
| QUMOND-type action | requires n ≤ 1; BCG excess would need n ≈ 2.3–2.5 |
| H0 = 70 → 67.4 | changes results by < 2 % |

## Layout

```
falsification/
├── paper/                 manuscript (PDF + source)
├── code/
│   ├── run_all.sh         reproduces every logged result
│   ├── data/              SPARC (Lelli+2016) tables and precomputed intermediates
│   ├── results/           logs (.txt) and results (.json) behind the manuscript numbers
│   └── *.py               analysis scripts (see below)
└── docs/                  data sources and provenance
```

## Reproduce

```bash
pip install -r requirements.txt       # numpy, scipy, scikit-learn
bash falsification/code/run_all.sh               # main results (~ 10 min)
bash falsification/code/run_all.sh --exploratory # also the earlier radiative-induction / lensing scripts
```

Scripts change to their own directory, so they can be run from anywhere. Results are written to
`falsification/code/results/`. Re-running reproduces the committed JSON files to 1e-6 relative tolerance.

## Scripts ↔ manuscript

| Script | Section | Produces |
|:--|:--|:--|
| `groups_test2.py <H0>` | 5 | 23 groups + 12 X-COP, G13/B14, LOO, cross-prediction, 200 ΛCDM mocks |
| `rgd_tests.py <H0>` | 3, 6, 7 | λ with fixed/marginalized M/L; A in galaxies; kernel monotonicity |
| `rgd_v5.py <H0>` | 4, 6 | RGD (v4) circularity mocks; predictions for SPARC and BCGs |
| `form_test.py` | 1.1 | BCGs vs SPARC at the same (g_b, α_b) |
| `audit3.py <σ_int>` | 3 | χ² likelihood with per-galaxy nuisances (Υd, Υb, D, i) |
| `audit2.py` | 3 | SPARC parser (159 galaxies / 3253 points after cuts) and CV evaluation |
| `clusters_test.py` | 1.1 | 20 CLASH BCGs (Tian+2020 / Postman+2012) |
| `groups_test.py` | — | earlier group test with a fixed 3 % stellar fraction (superseded; kept for history) |
| `analyze_radiative.py`, `fit_signed_exponents.py`, `audit.py`, `lens_test.py`, `lens_check.py`, `robust1.py` | — | exploratory radiative-induction / lensing-geometry audits that preceded this study |

## Data and provenance

See `docs/DATA_SOURCES.md`. Cluster and group tables were transcribed from the cited papers
(Sun+2009, Eckert+2019, Tian+2020, Postman+2012) and are embedded in the scripts.

## Known limitations

Stellar masses come from scaling relations (G13, B14), not per-group photometry; stars are
distributed like the total mass; X-COP redshifts are catalogue values. See the manuscript,
Reproducibility and data.

## How to cite

See `CITATION.cff` in this folder. Zenodo DOI to be added once the archived version is published.
