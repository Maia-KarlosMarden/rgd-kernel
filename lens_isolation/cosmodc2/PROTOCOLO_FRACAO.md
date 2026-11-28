# Análise POSTERIOR: por que o mock isola lentes demais? (escrito antes de rodar)

Conhecido: com vizinhos r<20 e erro de photo-z 0,02(1+z) o mock isola 44,7% das lentes com 9,0<=log Mb<11,5;
M24b informa ~16%; epsilon = 0,39+-0,06.

Hipótese H: o excesso de isoladas vem dos erros de photo-z independentes, que separam vizinhos reais.

Teste: vizinhos r<20; sz em {0,018; 0,010; 0,005; 0,002; 0} (0 = redshift do catálogo, com velocidade peculiar),
mesma semente, mesmo ruído de massa. Primeiro mede-se SÓ a fração isolada.
- H REFUTADA se nenhum sz levar a fração a 16% ou menos. Nesse caso nenhum epsilon novo é calculado.
- Se a fração cruzar 16%: sz* por interpolação linear (no cruzamento de maior sz) e só então epsilon em sz*,
  por faixa e no total, com a regra S/N<3.
- Cruzar 16% não prova H: só mostra que um erro efetivo menor basta para a fração total.
- Os números já declarados (0,543 e 0,391) não são substituídos.
