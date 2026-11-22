import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))  # paths relativos a v5/code
"""Robustez do sinal a_eff = a0·(α_b/α_ref)^n  (n≈0.08 em SPARC) a sistemáticas.
Para cada configuração: CV 5-fold por galáxia, χ² + incômodos (Υd, Υb, D, i), σ_int=0.057;
compara 'RAR' vs 'B2 α galáxia' na MESMA configuração. Uso: python robust1.py <grupo>
"""
import sys, io, contextlib, math, json
import numpy as np
from scipy.optimize import minimize
from sklearn.model_selection import KFold
GROUP = sys.argv[1] if len(sys.argv) > 1 else 'A'
sys.path.insert(0, '.')
with contextlib.redirect_stdout(io.StringIO()):
    import lens_test as LT
B = LT.B; A = B.A; ng = LT.ng
logL = np.log10(np.array([A.meta[nm]['L'] for nm in A.names])); logL_med = np.median(logL)
fD = np.array([A.meta[nm]['fD'] for nm in A.names]); Q = np.array([A.meta[nm]['Q'] for nm in A.names])
GAL_FULL = [g.copy() for g in B.gal]
CFG = dict(sig_ud=0.1, beta=0.0)

def mu_ud(i): return math.log10(0.5) + CFG['beta'] * (logL[i] - logL_med)
def prior(i, th):
    lud, lub, D, inc = th
    return ((lud - mu_ud(i)) / CFG['sig_ud'])**2 + ((lub - math.log10(0.7)) / 0.1)**2 + ((D - B.D0[i]) / B.eD[i])**2 + ((inc - B.i0[i]) / B.ei[i])**2
def bounds_i(i):
    s = max(0.3, 3 * CFG['sig_ud'])
    return [(mu_ud(i) - s, mu_ud(i) + s), (math.log10(.7) - .3, math.log10(.7) + .3),
            (max(B.D0[i] - 3 * B.eD[i], 0.05 * B.D0[i]), B.D0[i] + 3 * B.eD[i]), (max(B.i0[i] - 3 * B.ei[i], 10), min(B.i0[i] + 3 * B.ei[i], 90))]
def th_default(i): return np.array([mu_ud(i), math.log10(.7), B.D0[i], B.i0[i]])
B.prior = prior; B.bounds_i = bounds_i; B.th_default = th_default

def run(label, gals=None, seed=42, drop_inner=0, **cfg):
    CFG.update(dict(sig_ud=0.1, beta=0.0)); CFG.update(cfg)
    for i in range(ng):
        w = GAL_FULL[i]; o = np.argsort(B.r0[w]); B.gal[i] = w[o][drop_inner:] if len(w) - drop_inner >= 3 else w[o]
    gals = np.arange(ng) if gals is None else np.asarray(gals)
    folds = list(KFold(5, shuffle=True, random_state=seed).split(gals))
    models = ['RAR', 'B2 α galáxia']
    chi = {m: np.zeros(ng) for m in models}; nn = []
    for itr, ite in folds:
        for m in models:
            glob, center, _ = LT.train(m, list(gals[itr]))
            if m != 'RAR': nn.append(float(glob[1]))
            for i in gals[ite]:
                th = LT.fit_gal(i, m, glob, center, B.th_default(i))
                gb, y, sig, r, vb2, w = B.observables(i, th)
                chi[m][i] = np.sum(((LT.pred(glob, gb, LT.feat(m, i, th), center) - y) / sig)**2) + B.prior(i, th)
    d = chi['B2 α galáxia'][gals] - chi['RAR'][gals]
    srt = np.sort(d)  # mais negativos = mais ajudam B2
    top10 = float(srt[:10].sum()); tot = float(d.sum())
    # jackknife: Δχ² sem as 10 galáxias que mais favorecem B2
    res = dict(n_gal=int(len(gals)), dchi2=tot, n_folds=[round(x, 3) for x in nn], n_mean=float(np.mean(nn)),
               gal_melhoram=int(np.sum(d < 0)), dchi2_sem_top10=tot - top10)
    print(f"  {label:44s} gal {len(gals):3d}  Δχ²(B2−RAR) {tot:+7.1f}  sem top-10 {tot-top10:+7.1f}  melhoram {res['gal_melhoram']:3d}  n {np.mean(nn):+.3f} {res['n_folds']}", flush=True)
    return res

R = {}
if GROUP == 'A':
    print('== Grupo A: M/L e sementes ==')
    R['base s42'] = run('base (Υd prior 0.1 dex)')
    R['base s7'] = run('base, seed 7', seed=7)
    R['base s123'] = run('base, seed 123', seed=123)
    R['prior largo'] = run('Υd prior largo (0.25 dex)', sig_ud=0.25)
    R['beta+0.15'] = run('Υd ∝ L^+0.15 (M/L cresce com massa)', beta=0.15)
    R['beta-0.15'] = run('Υd ∝ L^−0.15 (M/L cai com massa)', beta=-0.15)
else:
    print('== Grupo B: subamostras ==')
    R['fD>=2'] = run('distâncias precisas (fD ≥ 2)', gals=np.where(fD >= 2)[0])
    R['fD=1'] = run('distâncias por fluxo de Hubble (fD = 1)', gals=np.where(fD == 1)[0])
    R['Q=1'] = run('qualidade Q = 1', gals=np.where(Q == 1)[0])
    R['sem 2 internos'] = run('sem os 2 pontos mais internos', drop_inner=2)
    R['L baixo'] = run('metade de menor L (anãs)', gals=np.where(logL <= logL_med)[0])
    R['L alto'] = run('metade de maior L', gals=np.where(logL > logL_med)[0])
json.dump(R, open(f'results/robust1_{GROUP}.json', 'w'), indent=1, ensure_ascii=False)
