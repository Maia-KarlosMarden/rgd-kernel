# Eficiência de isolamento no SMDPL UniverseMachine (z=0,315): resultado PRELIMINAR

Catálogo: sfr_catalog_0.760300.bin (86.190.067 registros). Rodado pelo autor com extrair_umachine.py e medir_caixa.py.
Protocolo: PROTOCOLO.md + EMENDA_SMDPL.md. Não decide o critério registrado no preprint.

## epsilon, 0,3-1 Mpc (todas as lentes)
base 0,448±0,061 | pz018 0,441±0,063 | m0 0,324±0,053 | m020 0,487±0,056 | obs_sm 0,498±0,063
Referência: janela 0,43-0,68. Base na borda inferior; indecidido.

## Tendência com a massa (combinação das duas faixas baixas menos as duas altas)
base -0,08±0,13 | pz018 +0,16±0,10 | m0 +0,12±0,09 | m020 +0,13±0,09 | obs_sm +0,13±0,09
Nenhuma tendência significativa (no máximo ~1,7 sigma). Inconclusivo.

## Limitações
- Valores por faixa mudam até ~0,25 entre variantes próximas: instáveis; o erro de jackknife não cobre isso.
- Sem limite de magnitude: só ~1% das lentes da faixa mais baixa passam o corte (187 de 20000), contra ~37% na mais alta.
  A faixa baixa não representa a amostra real.
- Satélites: o halo hospedeiro da lente foi excluído do ambiente (falha do protocolo, notada após o resultado).
- Caixa num instante; massa de halos em cascas 3D, não Delta-Sigma; "todas" não segue a abundância real.
- Inclinação de xi não medida.

## Análise POSTERIOR (EMENDA2_SATELITES.md): halo hospedeiro das satélites incluído (variante B)
Rodada pelo autor com versão curta (medir_caixa.measure com host=-1 nas satélites).
base 0,557±0,022 | pz018 0,550±0,022 | m0 0,507±0,024 | m020 0,531±0,026 | obs_sm 0,526±0,024
O número registrado continua 0,448±0,061. Variante C, fração de satélites e erro pareado de B-A: não medidos.
Ressalva adicional: a massa do hospedeiro entra como um ponto no seu centro, o que pesa mais nesta variante.
