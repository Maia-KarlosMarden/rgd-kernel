# CLASH: previsão registrada antes de abrir as massas por lenteamento (2026-11-27)

## Hipótese
H-host: o reforço gravitacional depende da massa bariônica do sistema hospedeiro.

## Definições
- Variável do hospedeiro S = log10 M_gas(<500 kpc) (clusters). SPARC: log10 da massa bariônica no último ponto
  (supõe que não há bárions relevantes além dele; halo circungaláctico ignorado).
- Forma funcional: a H5 já existente (largura s=0,3 dex fixa; transição S0 ajustada). Nenhuma forma nova.
- Calibração: só SPARC e X-COP (12 clusters, M_FORW), catálogos com peso igual. BCG, grupos e CLASH fora.
- CLASH: gás ACCEPT2.0 (H0=71, Om=0,27) convertido para 70/0,30 (convenção do X-COP usada na calibração;
  g† continua com H0=67,4; diferença em g_b ~0,008 dex). mu_e=1,17. Sem correção do desvio ACCEPT2.0-X-COP (-0,023 dex).
- Raios: 100,150,200,300,400,500,600 kpc, sem extrapolar além do perfil.
- Estrelas: M*/M_gas(r) mínimo / mediana / máximo dos 7 clusters X-COP com perfil medido.
- Controles: CTRL-local (Phi=g_b r; classe já falseada), RAR recalibrada em SPARC+X-COP, RAR e RGD só SPARC,
  mocks LCDM (halo definido só pelo gás, f_gas=0,10/0,12/0,14).
- Gás nos mocks além do último raio: lei de potência ajustada ao terço externo, inclinação limitada a [0; 1,5].

## O que se sabe antes (importante)
- A H5 com logMb_out JÁ FALHOU na rodada de correções (C1 0,211 > 0,208; C2 RMS 0,147 > 0,10; não monotônica;
  grupos piores que a RAR em 0,054 dex). Isto é um teste de falsificação, não a validação de um candidato aprovado.
- Não é cego: BCGs do CLASH (14-30 kpc) já foram vistas no run01. Os raios >=100 kpc não foram usados.
- Calibração em massas hidrostáticas; se houver viés hidrostático, o lenteamento ficará acima (esperado 0,05-0,10 dex).
  Não será corrigido depois.

## Fonte observacional a abrir depois
Umetsu et al. 2016 (ApJ 821, 116), M200c e c200c NFW por cluster; g_obs = G M_NFW(<r)/r², convertido para 70/0,30.
Perfil paramétrico: a inclinação radial depende do modelo.
Primária: 16 selecionados por raios X. Estresse: MACSJ0416, J0647, J0717, J1149.

## Critérios (e = log10 K_prev - log10 K_obs; RMS por sistema)
T1 transporte (primária, estrelas na mediana):
  PASSA se RMS<=0,15 E |viés|<=0,10 E RMS menor que CTRL-local e que RAR(SPARC+X-COP) em >=0,05 dex.
  FALHA se RMS>0,20 OU |viés|>0,15. Caso contrário INCONCLUSIVO.
  Se o veredito mudar entre os três cenários estelares: INCONCLUSIVO.
T2 inclinação radial: apenas descritivo (NFW paramétrico).
T3 massa do hospedeiro dentro dos clusters: inclinação de Theil-Sen de (obs - prev) a 300 kpc contra S, IC 95% bootstrap.
  IC que exclui zero => a H-host nesta forma não captura a dependência.
T4 LCDM: se o RMS da mediana dos mocks (f=0,12) for <= RMS da H-host, não há evidência além de LCDM,
  mesmo que T1 passe.
Mesmo com T1 e T4 favoráveis, H-host não fica validada: exige ação não local explícita (RGD_Falsification §9).
