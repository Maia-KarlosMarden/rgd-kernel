# Teste relacional multirraio no X-COP (26/11/2026)

Protocolo: `PROTOCOLO.md` (hash em `PROTOCOLO.sha256`). Previsões e mocks congelados antes de abrir as massas totais (`E1_congelado.sha256`). O registro é interno à sessão (sem carimbo de tempo externo) e **não é cego**: massas do X-COP em dois raios e o desvio publicado em relação à RAR já eram conhecidos.

## Dados
- `allfiles.tar.gz`, SHA-256 0edf5038...0d94. 8 raios fixos, 100–1000 kpc. x mediano de 0,27 a 0,087.
- Primária: A85, A644, A1795, A2029, A2142, ZW1215 (estrelas medidas). A2319 como estresse. Os 12 como secundária (estrelas = 0,09 do gás onde não medidas).
- M_estrelas/M_gás mediana (primária): 0,575 a 100 kpc, 0,067 a 1 Mpc.
- Comparação da massa NFW do arquivo com a Tabela 1 de Ettori (valores transcritos no run01): 9 clusters diferem <0,15%; A644 (+3,5/+4,4%), A2029 (−1,4/−1,9%) e A2319 (+3,5/+4,2%) diferem. Coincide com o relato público de terceiros (extragravitation/atlas, issue 19).

## Resultados pelos critérios registrados
| Critério | Veredicto |
|---|---|
| P1: contraste de inclinação clusters − galáxias | **PASSA** (primária, 5 massas, com e sem A2319) |
| P2: alguma lei fixada no SPARC explica os clusters | **FALHA** (todas) |
| P3/P4: contraste exige algo além do ΛCDM | **INCONCLUSIVO** |

### P1
Primária, M_FORW: I1 dos clusters −0,31 [−0,41; −0,24]; D = −0,33 [−0,49; −0,21]. R1: D = −0,92 [−1,12; −0,66]. Galáxias (79, mesma janela de x): I1 = +0,01. IC por bootstrap com 6 clusters: grosseiro. Na amostra secundária de 12, o R1 não é estável (IC inclui zero em M_FORW e M_NFW).

### P2: leis transferidas (primária, M_FORW)
| Lei | RMS log10 K | viés (previsto − observado) | I1 previsto |
|---|---|---|---|
| RGD | 0,493 | −0,485 | 0,000 |
| RAR | 0,498 | −0,490 | −0,002 |
| Spline | 0,476 | −0,468 | −0,046 |
| F3 | 0,474 | −0,467 | −0,068 |
| F3+G | 0,489 | −0,481 | −0,029 |
| G1 | 0,509 | −0,500 | +0,035 |
| N2 | 0,527 | −0,515 | +0,126 |
| F3+G+N | 0,506 | −0,496 | +0,062 |

F3 e spline ganham ~0,02 dex sobre RGD (cerca de 4%, abaixo dos 5% exigidos) e continuam errando por um fator ~3. Nenhuma mediana prevista de I1 cai no IC observado. G e N com coeficientes do SPARC pioram. Não houve reajuste em clusters, portanto isto falsifica a **transferência**, não toda forma possível de geometria ou não localidade.

### P3: ΛCDM
Previsão registrada antes das massas (f=0,12): I1 = −0,38 (IC 95% do mock −0,51 a −0,26). Observado com M_FORW: −0,31. Nível: log10 K observado 0,92, mock 0,92–1,01 conforme f.

p bicaudal de I1, primária, f = 0,10 / 0,12 / 0,14:
| Massa | I1 obs | p |
|---|---|---|
| FORW | −0,31 | 0,66 / 0,30 / 0,12 |
| NFW | −0,38 | 0,47 / 0,91 / 0,65 |
| EIN | −0,53 | 0,003 / 0,018 / 0,059 |
| ISO | −0,27 | 0,26 / 0,10 / 0,028 |
| BUR | −0,32 | 0,79 / 0,38 / 0,16 |

Com a massa primária o ΛCDM é compatível; com Einasto (sobretudo somando A2319) não. Pela regra P4, inconclusivo.

## Sensibilidades pós-hoc (não registradas)
| Variante | I1 | D | log10 K | excesso sobre RAR |
|---|---|---|---|---|
| registrada | −0,31 | −0,33 | 0,92 | 0,52 |
| sem estrelas | −0,37 | −0,39 | 1,00 | 0,55 |
| estrelas ×2 | −0,28 | −0,30 | 0,87 | 0,47 |
| estrelas ×3 | −0,26 | −0,28 | 0,83 | 0,43 |
| viés hidrostático 5%→15% | −0,26 | −0,27 | 0,98 | 0,55 |
| viés hidrostático 0→30% | −0,14 | −0,15 | 1,01 | 0,57 |

- Um viés hidrostático que cresce com o raio pode reduzir o contraste pela metade. O sinal e o nível do excesso se mantêm.
- Com ruído de mock igual ao erro mediano real (3,2% em vez de 10%), a leitura não muda (`PH_mock_ruido_real.csv`).
- I1 recalculado por um segundo código, a partir dos FITS: −0,3137 nos dois.

## Limitações
Seis clusters; equilíbrio hidrostático; sem luz intra-aglomerado; sem covariância entre raios; mocks com halo NFW esférico e massa estimada pelo gás, sem hidrodinâmica; incertezas de bárions não propagadas. P+ depende do raio onde o perfil de gás termina.

## O que muda em relação ao run01
- A quebra aparece com 8 raios, mas com metade da amplitude (−0,33 contra −0,62 com dois raios).
- A não localidade em clusters deixa de estar bloqueada quanto à transferência: falha.

## Reproduzir
Extrair o arquivo em `/home/claude/xcop_in` (ou editar `DATA` em `xcop_multi/common.py`); dentro de `xcop_multi/`, rodar `stage1.py`, `stage2.py`, `posthoc.py`. Requer o pacote do run01 e astropy.
