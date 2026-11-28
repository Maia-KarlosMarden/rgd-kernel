# Eficiência de isolamento de lentes: medidas em mocks

Pergunta (preprint "isolation efficiency"): que fração epsilon do excesso de massa ambiental o corte de
isolamento remove num universo CDM? Janela registrada: 0,43-0,68, a ser julgada com Delta-Sigma em mock com
cone de luz. **O teste registrado ainda não foi feito.**

| Pasta | Situação em 2026-11-28 |
|---|---|
| `mice_not_run/` | Protocolo e código; o catálogo MICE nunca foi obtido. |
| `smdpl/` | Rodado (caixa, um instante). Preliminar: 0,448±0,061 como especificado; 0,557±0,022 em análise posterior. |
| `cosmodc2/` | Protocolo e código publicados; resultado ainda não registrado aqui. |

## Ordem dos fatos (o que a data deste commit prova e o que não prova)
- SMDPL: os arquivos só foram publicados **depois** das medições. Os hashes foram gerados na sessão de trabalho,
  sem testemunho externo.
- cosmoDC2: protocolo e código escritos antes de qualquer análise. Este commit foi feito **enquanto a medida
  rodava** no computador do autor, antes de o autor relatar qualquer número. Não é, portanto, um pré-registro
  público estritamente anterior à execução.
- As escolhas do cosmoDC2 (incluindo a regra das satélites) foram feitas já conhecendo os números do SMDPL.

Os catálogos não são redistribuídos. cosmoDC2: portal.nersc.gov/project/lsst/cosmoDC2 (LSST DESC).
SMDPL UniverseMachine: ver `smdpl/extrair_umachine.py`. Verificação: `sha256sum -c HASHES*.sha256`
(nas pastas, para os arquivos presentes).
