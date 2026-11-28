# Uso: py fracao_dc2.py   (precisa de medir_dc2.py e dos dc2_*.npz na pasta)
import sys, json, glob
import numpy as np
import medir_dc2 as M

def frac(S):
    k = (S['lmb']>=9.0)&(S['lmb']<11.5)
    return float(S['iso'][k].mean())

def main(files, alvo=0.16):
    d = M.load(files); area = M.Area(d['ra'], d['dec']); pos, Xinc, Xexc, _ = M.environment(d, area)
    grid = [0.018, 0.010, 0.005, 0.002, 0.0]; f = []
    for sz in grid:
        S = M.select(d, area, pos, sz=sz, rviz=20.0); f.append(frac(S))
        print('sz %.3f  fracao isolada %.1f%%  (lentes %d)' % (sz, 100*f[-1], len(S['lens'])), flush=True)
    res = {'alvo': alvo, 'sz': grid, 'fracao': f}
    szs = None
    if f[0] <= alvo: szs = grid[0]
    else:
        for i in range(1, len(grid)):
            if f[i] <= alvo:
                szs = grid[i-1] + (alvo-f[i-1])*(grid[i]-grid[i-1])/(f[i]-f[i-1]); break
    if szs is None:
        print('H REFUTADA: nenhum sz chega a %.0f%%' % (100*alvo))
    else:
        S = M.select(d, area, pos, sz=szs, rviz=20.0)
        a = M.analyse(S, Xinc[pos[S['lens']]], len(files))
        res.update(sz_estrela=szs, fracao_em_sz_estrela=frac(S), resultado=a)
        print('sz* = %.4f  fracao %.1f%%' % (szs, 100*frac(S)))
        for n,v in a.items():
            if isinstance(v, dict):
                print(n, 'eps', v.get('eps_0.3_1Mpc'), '+-', v.get('erro'), 'S/N', v.get('sn_excesso_controle'),
                      'isoladas', v.get('n_isoladas'), 'de', v.get('n_lentes', '-'))
            else: print(n, v)
    json.dump(res, open('fracao_resultado.json','w'), indent=1)

if __name__=='__main__':
    files = sorted(glob.glob('dc2_*.npz'))
    if len(files) != 30: sys.exit('esperava 30 recortes, achei %d' % len(files))
    main(files)
