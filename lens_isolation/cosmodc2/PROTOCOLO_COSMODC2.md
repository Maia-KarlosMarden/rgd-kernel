# Protocolo cosmoDC2 (escrito antes de extrair os dados)

Já se conhecem os resultados do SMDPL (0,448 como especificado; 0,557 posterior com hospedeiro das satélites).
As escolhas abaixo foram feitas conhecendo esses números. Nenhum dado de lente do cosmoDC2 foi analisado.
O critério do preprint (janela 0,43-0,68) continua sendo julgado como PRELIMINAR: ele pede Delta-Sigma,
que aqui não é confiável em 0,3-1 Mpc.

## Área
30 recortes HEALPix (nside 32, ring), ~100 deg2, RA 54,8-67,5, Dec -40,2 a -30,0:
9685 9558 9814 9686 9557 9813 9684 9429 9941 9559 9430 9815 9942 9428 9556 9687 9940 9812 9302 9301
10070 10069 9431 9683 9303 9943 9560 9816 9173 10071

## Medida principal: contagem de massa de halos
- Lentes: r_SDSS (observada, com poeira) < 20 e 0,1 < z_obs < 0,5. Vizinhos: r < 22,5 (variante: r < 20).
- z_obs = `redshift` + N(0; 0,02(1+z)) (variante 0,018). log M* = log totalMassStellar + N(0; 0,12) (variantes 0; 0,20).
- Corte: rejeita se houver vizinho com M* > 0,1 M*_lente a menos de 4/h70 Mpc, 3D comóvel; H0=73, Om=0,2793.
- Massa ambiental: galáxias centrais com hostHaloMass > 1e11, posições verdadeiras (ra_true, redshiftHubble),
  cascas comóveis; número principal em 0,3-1 Mpc.
- SATÉLITES (decidido de antemão): principal = o halo hospedeiro da lente satélite ENTRA no ambiente.
  Secundária: excluído. Para lentes centrais, o próprio halo é sempre excluído.
- Controle: todas as lentes, reponderadas em (log M* observado, z_obs). Aleatórios dentro da área.
- Margem de borda de 4 Mpc projetados. Jackknife por recorte. Inconclusivo se S/N do excesso do controle < 3.
- Faixas log Mb: 9,0 / 10,5 / 10,8 / 11,1 / 11,5, agora com a abundância do próprio catálogo.

## Medida secundária: Delta-Sigma
- Fontes: 5% aleatórios das galáxias com 0,4 < z < 1,0, posições lenteadas (ra, dec); z_fonte > z_lente + 0,2.
- Convenção de sinal de shear2 fixada ANTES, em halos > 1e13: sinal tangencial positivo e componente cruzada nula.
- Só R entre 1 e 3 Mpc: eps = 1 - DS_isoladas/DS_controle. Aproximada (inclui resíduo do halo próprio).
  Não decide o critério.

## A conferir na extração, antes de analisar
Unidades de totalMassStellar e hostHaloMass (comparar a distribuição com o esperado em Msun).

## Detalhes de implementação (fixados com o código, antes de ver qualquer resultado)
- Área útil: células de 0,05 grau ocupadas por galáxias. 5 pontos aleatórios por lente.
- Margem de borda calculada com o menor entre z_obs e z verdadeiro.
- S/N do controle = excesso / erro jackknife (30 recortes). Delta-Sigma só na variante base, com z verdadeiros
  e sem subtração de sinal em pontos aleatórios.
- Tendência de massa: média das 2 faixas baixas menos média das 2 altas, com erro jackknife.
- Código testado só em catálogo artificial: sem corte dá 0; sinal do shear recuperado nas duas convenções.
