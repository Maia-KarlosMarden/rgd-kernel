# Data sources

| Data | Source | Used in |
|:--|:--|:--|
| SPARC rotation curves and photometry (175 galaxies) | Lelli, McGaugh & Schombert 2016, AJ 152, 157 (`data/SPARC_*.txt/.mrt`) | galaxies, RAR, BTFR |
| Nuisance priors (Υ, distance, inclination) | Li, Lelli, McGaugh & Schombert 2018, A&A 615, A3 | `audit3.py`, `rgd_tests.py` |
| 23 groups (Tier 1+2), M500, f_gas,500 | Sun et al. 2009, Table 3 | `groups_test2.py` |
| 12 X-COP clusters | Eckert et al. 2019, A&A 621, A40; Ettori et al. 2019 | `groups_test2.py` |
| CLASH BCG accelerations | Tian et al. 2020; Postman et al. 2012 (L_bol) | `clusters_test.py` |
| M500 from kT | Arnaud, Pointecouteau & Pratt 2005 | `clusters_test.py` |
| Stellar fraction f*(M500) | Gonzalez et al. 2013 (G13); Budzynski et al. 2014 (B14) | `groups_test2.py` |

`SPARC_tripla_test_results.json` and `SPARC_expoentes_livres_resultados.json` are intermediate
outputs of `analyze_radiative.py` / `fit_signed_exponents.py`, shipped so the later scripts run
without recomputing them. Cluster/group tables are transcribed by hand; the transcription of
Sun+2009 Table 3 was checked row by row against the arXiv HTML version.
