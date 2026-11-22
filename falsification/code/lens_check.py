import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))  # paths relativos a v5/code
"""Checagens: (a) SPARC com n FIXO no valor exigido pelos clusters, CV completa; (b) outra semente;
(c) controle sem luz: a_eff ∝ (V_b,max²)^n  — mesma informação de 'profundidade de potencial' sem geodésica."""
import sys, io, contextlib, math, json
import numpy as np
from scipy.optimize import minimize
from sklearn.model_selection import KFold
sys.path.insert(0, '.')
with contextlib.redirect_stdout(io.StringIO()):
    import lens_test as LT
B = LT.B
NFIX = {}
orig_feat = LT.feat
def feat(model, i, th):
    if model.startswith('B2'): return orig_feat('B2 α galáxia', i, th)
    if model.startswith('B1'): return orig_feat('B1 α local', i, th)
    if model == 'CTRL Vb²max':
        gb, y, sig, r, vb2, w = B.observables(i, th); return np.full(len(w), np.log10(vb2.max()))
    return orig_feat(model, i, th)
LT.feat = feat
orig_fg = LT.fit_global
def fit_global(model, gals, nuis, center, g0):
    if model in NFIX:
        data = [(B.observables(i, nuis[i]), feat(model, i, nuis[i])) for i in gals]
        ob = lambda g: sum(np.sum(((LT.pred(np.array([g[0], NFIX[model]]), d[0], F, center) - d[1]) / d[2])**2) for d, F in data)
        return np.array([minimize(ob, g0[:1], method='L-BFGS-B', bounds=[(-12, -8)]).x[0], NFIX[model]])
    return orig_fg(model, gals, nuis, center, g0)
LT.fit_global = fit_global
def cv(models, seed):
    folds = list(KFold(5, shuffle=True, random_state=seed).split(np.arange(LT.ng)))
    chi = {m: 0.0 for m in models}; nn = {m: [] for m in models}
    for itr, ite in folds:
        for m in models:
            glob, center, _ = LT.train(m, list(itr)); nn[m].append(round(float(glob[1]), 3) if len(glob) > 1 else 0)
            for i in ite:
                th = LT.fit_gal(i, m, glob, center, B.th_default(i))
                gb, y, sig, r, vb2, w = B.observables(i, th)
                chi[m] += np.sum(((LT.pred(glob, gb, feat(m, i, th), center) - y) / sig)**2) + B.prior(i, th)
    return chi, nn
R = {}
NFIX.update({'B2 n=0.6 (exigido pelos clusters)': 0.6, 'B2 n=0.3': 0.3})
for seed, models in [(42, ['RAR', 'B2 n=0.6 (exigido pelos clusters)', 'B2 n=0.3', 'CTRL Vb²max']), (7, ['RAR', 'B1 α local', 'B2 α galáxia', 'CTRL Vb²max'])]:
    chi, nn = cv(models, seed)
    print(f'\n== seed {seed} ==')
    for m in models:
        print(f"  {m:36s} χ² {chi[m]:8.1f}  Δχ² vs RAR {chi[m]-chi['RAR']:+8.1f}   n {nn[m]}")
        R[f'{seed}:{m}'] = dict(chi2=chi[m], dchi2=chi[m] - chi['RAR'], n=nn[m])
json.dump(R, open('results/lens_check_resultados.json', 'w'), indent=1, ensure_ascii=False)
