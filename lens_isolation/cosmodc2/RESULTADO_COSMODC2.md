# cosmoDC2: resultado PRELIMINAR (28/11/2026)

Rodado pelo autor com medir_dc2.py (hash conferido por ele antes de rodar). Números transcritos da saída do terminal.
Não decide o critério registrado (sem Delta-Sigma confiável em 0,3-1 Mpc; inclinação de xi não medida).

Diagnóstico: 30 recortes, 102,2 deg2; r<20: 1075/deg2; r<22,5: 8877/deg2. Shear: sinal +1
(g_t=+0,00322; g_x=0,000001). Satélites com hospedeiro a menos de 3 Mpc: 99,9%.

## epsilon em 0,3-1 Mpc
| variante | isoladas / lentes | principal (com hospedeiro) | sem hospedeiro |
|---|---|---|---|
| base | 17431 / 86165 | 0,543 +- 0,050 | 0,721 +- 0,103 |
| pz018 | 17433 / 86287 | 0,575 +- 0,040 | 0,729 +- 0,140 |
| m0 | 17309 / 86165 | 0,532 +- 0,045 | 0,702 +- 0,086 |
| m020 | 17711 / 86165 | 0,523 +- 0,056 | 0,738 +- 0,137 |
| viz_r20 | 39005 / 86165 | 0,391 +- 0,062 | 0,678 +- 0,169 |

Por faixa de log Mb (base): 9,0-10,5: 0,695+-0,071; 10,5-10,8: 0,650+-0,081; 10,8-11,1: 0,595+-0,059;
11,1-11,5: 0,374+-0,051. Baixas menos altas: +0,188+-0,059.

Delta-Sigma (secundário, 1-3 Mpc, base): epsilon = 0,348+-0,033 (isoladas 1,482; controle 2,273 Msun/pc2).

## Ressalvas (identificadas DEPOIS de ver o resultado)
1. A variante principal usa vizinhos até r<22,5. Mistele et al. descrevem o limite m=20 da amostra KiDS-bright
   como causa de vizinhos não detectados, o que indica que o corte real usa vizinhos r<20, isto é, a variante
   viz_r20 (0,391). A escolha da principal foi um erro do protocolo. O número declarado não é trocado;
   os dois devem ser citados juntos.
2. Fração isolada: M24b informa ~16%. Aqui: 20% (base) e 45% (viz_r20), com denominadores não idênticos.
   A seleção do mock ainda NÃO reproduz a real. Hipótese não testada: erros de photo-z independentes e
   gaussianos desfazem pares físicos.
3. Faixa por faixa e tendência de massa da variante viz_r20: ainda não examinadas.
