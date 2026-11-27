# X-COP multi-radius relational test (2026-11-26)

Provenance for Sec. 1 (X-COP row) and the first bullet of Sec. 3 of
[`../clash-2026-11-27/ADDENDUM.md`](../clash-2026-11-27/ADDENDUM.md).
Protocol and report are in Portuguese: `runs/xcop_multi/PROTOCOLO.md`, `RELATORIO_XCOP.md`.

## Status of the freeze
Predictions and mocks were hashed before the total masses were opened (`E1_congelado.sha256`), but the
record was internal to the working session: **no external timestamp**, and **not blind**. These files
were published here after the results were known.

## What is included, and what is not
Protocol, report, hashes, scripts and per-cluster or aggregate results. Point-level files derived
from the X-COP archive (`E1_barions.csv`, `E1_previsoes_logK.csv`, `E2_pontos.csv`) and
`E1_mocks_lcdm.csv` are not included, so `E1_congelado.sha256` cannot be fully verified from this
folder. `PROTOCOLO.sha256` and the hashes of the listed files present here can be:

    cd runs/xcop_multi && sha256sum -c --ignore-missing E1_congelado.sha256 PROTOCOLO.sha256

The scripts need the pipeline modules and the X-COP archive (https://dominiqueeckert.wixsite.com/xcop/data),
which are not included; local paths in them are those of the original run.
