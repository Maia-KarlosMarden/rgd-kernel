# Pré-registro: teste relacional multirraio no X-COP (arquivo público allfiles.tar.gz)
SHA-256 do arquivo: 0edf5038b419b70d070b73b22f4801e27f318b0854db61eec52142c27c140d94

## O que já era conhecido ao escrever (NÃO é cego)
- Massas totais do X-COP em 2 raios e K nesses raios (run01); contraste I1 clusters-galáxias ~ -0,62.
- Resultado publicado de que a RAR dos X-COP se afasta da das espirais.
- Lidos até agora do arquivo: nomes, colunas, unidades, comentários e faixas das colunas de RAIO. Nenhum valor de massa, densidade, temperatura.

## Ordem de execução
E1: ler SÓ bárions (MGAS*, MSTAR*, raios, R500 do header apenas para converter R/R500 em kpc). Calcular e congelar previsões e mocks (hash).
E2: só então ler M_tot (hydro_mass).

## Amostras
- Primária: A85, A644, A1795, A2029, A2142, ZW1215 (estrelas medidas).
- Estresse: A2319 (fusão), reportado à parte e somado.
- Secundária: os 12, com Mstar = 0,09 Mgas nos sem perfil estelar (e, para comparação, nos 12).
- Estrelas: MSTAR_SMOOTHED. ICL não incluída.

## Variáveis
- 8 raios fixos, log-espaçados de 100 a 1000 kpc (independentes da massa total). Interpolação log-log.
- g_b = G(Mgas+Mstar)/r²; x = g_b/g†, g† com H0=67,4 e z=0 (como no run01).
- K = M_tot/M_b. Primária: M_FORW. Sensibilidade: M_NFW, M_EIN, M_ISO, M_BUR.
- R0 = ln(K-1) + x - ln λ - ½L, L = -ln(1-e^-x). I1 = inclinação de R0 em ln r; R1 = inclinação em L (mínimos quadrados, 8 raios). λ do RGD ajustado no SPARC.
- Galáxias: mesmas I1, R1 no SPARC, pontos com x entre os percentis 5 e 95 do x dos clusters primários, >=4 pontos por galáxia.

## Leis fixadas só no SPARC (sem reajuste em clusters)
RGD, RAR, SPL, F3, G1, F3+G, N2, F3+G+N. Z = (1 - dln g_b/dln r)/3. P-, P+ = integrais de g_b do primeiro ao último raio do perfil de gás / (r g_b).

## Mocks ΛCDM (2000)
Bárions reais + halo NFW: M500 obtido SÓ da massa de gás (raio onde Mgas = f * (4π/3) R³ 500 ρc), f = 0,12 (variantes 0,10 e 0,14), M200 = 1,45 M500, dispersão 0,15 dex, c(M) como no run01 (dispersão 0,15 dex); g = g_b + G (1-0,157) M_NFW/r². Ruído por raio: independente, 10% (a covariância não é conhecida).

## Critérios
- P1 (existência): D = mediana I1(clusters) - mediana I1(galáxias) < 0, com IC95% (bootstrap por sistema) excluindo zero, amostra primária, M_FORW. Idem R1.
- P2 (lei transferida explica): na primária, a lei supera RGD e RAR em RMS de log10 K por >=0,005 dex E >=5% (IC95% pareado por sistema, Holm sobre as leis) E sua mediana prevista de I1 cai no IC95% da observada. Falhar qualquer parte = FAIL.
- P3 (exclusividade): D observado e mediana I1 dos clusters "não explicados por ΛCDM" só se p bicaudal <= 0,05 em TODAS as variantes de f.
- P4 (robustez): mesma conclusão de P1-P3 com todas as cinco massas e com/sem A2319. Se mudar: INCONCLUSIVE.
