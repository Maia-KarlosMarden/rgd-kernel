# Emenda ao protocolo: SMDPL UniverseMachine (escrita antes de ver os dados)

MOTIVO: sem acesso ao MICE. NATUREZA: preliminar. Não decide o critério registrado no preprint
(não há cone de luz nem lenteamento). A janela 0,43-0,68 é mostrada só como referência.

Mudanças em relação ao PROTOCOLO.md:
- Caixa periódica, um instante (a=0,7603; z=0,315), h=0,6777, L=400 Mpc/h.
- Linha de visada = eixo z da caixa; erro de photo-z convertido em distância comóvel: c*sz*(1+z)/H(z),
  H(z) com Om=0,307. Velocidades peculiares ignoradas.
- Massa estelar: `sm` (verdadeira) + ruído gaussiano de 0,12 dex. Variantes: 0; 0,20; `obs_sm` sem ruído extra.
- Sem limite de magnitude. Vizinhos: sm > 10^7,9 Msun. Lentes: 9,0 <= log Mb < 11,5.
- Até 20000 lentes sorteadas por faixa; o número "todas" é a união delas (NÃO segue a abundância real).
- Massa ambiental: halos centrais com M > 10^11 Msun/h, excluído o halo da lente.
- Controle reponderado só em massa estelar observada. Jackknife em 50 fatias no eixo x.
- Resultado inconclusivo se o excesso do controle em 0,3-1 Mpc tiver S/N < 3.
