# Replicação de I1 no lenteamento do CLASH (registrado em 2026-11-27, antes do cálculo)

Exploratória e NÃO cega: as massas do CLASH já foram abertas e já se viu que d ln(K-1)/d ln r observado é ~0,05
(que não é I1). O valor anterior do X-COP (I1=-0,31, 100-1000 kpc) também é conhecido.

## Medida
- I1 = dR0/dln r por cluster, função src.diagnostics.slopes (a mesma do X-COP), lambda da RGD no SPARC; mediana entre clusters.
- Raios comuns: 100,150,200,300,400,500,600 kpc (CLASH: até onde há gás).
- CLASH: 15 primários, g_obs do ajuste NFW de Umetsu 2016, estrelas no cenário mediano (baixo e alto = sensibilidade).
- X-COP: 6 primários, M_FORW, mesmos raios, estrelas medidas.
- IC 95%: bootstrap por cluster (10000, semente 1).

## Critério
REPLICA: IC do CLASH inteiro abaixo de zero E com sobreposição ao IC do X-COP nos mesmos raios.
NAO REPLICA: IC do CLASH inclui zero ou está acima E não se sobrepõe ao do X-COP.
INCONCLUSIVO: demais casos; também se o IC do X-COP nos mesmos raios incluir zero, ou se o veredito mudar com o cenário estelar.
## Controle LCDM
Mesmos mocks (halo pelo gás, f=0,10/0,12/0,14, 2000 realizações, ruído 0,10/ln10 por ponto). Se a mediana observada
cair no intervalo central de 95% dos mocks, o valor não distingue LCDM. Reportar os três f.
R1 é apenas descritivo. Ressalva: o NFW impõe parte da forma radial; nem concordância nem discordância são decisivas.
