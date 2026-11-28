# Protocolo: eficiência de isolamento no MICE (fixado antes de ver o catálogo)

Critério do preprint (inalterado): epsilon_mock fora de 0,43-0,68 desfavorece CDM com c(M) padrão
(fora de 0,30-0,75 se a inclinação de xi no mock diferir de 1,8).

ESTA MEDIDA É PRELIMINAR. O critério registrado pede o Delta-Sigma ambiental. Aqui epsilon(r) usa a massa de
halos vizinhos em cascas 3D. O veredito final só vem com o Delta-Sigma empilhado.

- Corte: rejeitar a lente com vizinho de M* > 0,1 M*_lente a menos de 4/h70 Mpc (3,836 Mpc), 3D, comóvel,
  com distância radial dada pelo photo-z. Geometria H0=73, Om=0,2793.
- Photo-z: z + N(0; 0,02(1+z)). Variante 0,018.
- Ruído em log M*: 0,12 dex (SUPOSTO). Variantes 0 e 0,20.
- Lentes: r<20, 0,1<z_phot<0,5. Vizinhos: todos com r<22,5. Margem de borda: 4 Mpc projetados.
- Excesso(<casca) = massa de halos centrais (exceto o da lente), posições verdadeiras, menos o mesmo em
  pontos aleatórios de mesmo redshift.
- Controle: todas as lentes, sem corte, reponderadas para a mesma distribuição de (M*, z) das isoladas.
- epsilon(r) = 1 - E_isoladas/E_controle. Número principal: razão dos excessos somados em 0,3-1 Mpc.
- Faixas: log Mb = 9,0 / 10,5 / 10,8 / 11,1 / 11,5; Mb = M*(1+f), f = 10^(-0,69 log M* + 6,63) (A CONFERIR).
- Erros: jackknife em 50 faixas de ascensão reta.

A conferir no portal antes de rodar: nomes das colunas, convenção de flag_central (suposto 0 = central),
unidades de lmstellar/lmhalo (suposto Msun/h^2 e Msun/h, h=0,7), e se o corte de M24b é comóvel ou físico.

Regra adicional: se o excesso do controle em 0,3-1 Mpc não for positivo a mais de 3 sigma, o resultado é INCONCLUSIVO.
