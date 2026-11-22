import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))  # paths relativos a v5/code
"""Fase 1c — perda χ² e marginalização (perfil MAP) de parâmetros de incômodo por galáxia.
Incômodos por galáxia, com priors (como Li, Lelli, McGaugh & Schombert 2018, A&A 615, A3):
  log10 Υ_disk ~ N(log10 0.5, 0.1);  log10 Υ_bul ~ N(log10 0.7, 0.1)
  D ~ N(D0, eD);  i ~ N(i0, e_i)  (truncados em ±3σ, i ∈ [10°, 90°])
Escalas: r ∝ D; V_obs ∝ sin i0 / sin i; V_b² ∝ D (g_b independe de D); g_obs ∝ 1/D.
σ_y = (2/ln10)·eV/V  [+ σ_int opcional].
CV 5-fold por galáxia: parâmetros globais (a0 e expoentes) ajustados no treino por otimização
alternada (incômodos ↔ globais); no teste, globais fixos e só os incômodos de cada galáxia
de teste são ajustados (mesma liberdade para todos os modelos). Métrica: χ² perfilado de
teste (χ² + prior) e RMSE por galáxia (dex) após os incômodos.
"""
import math, json, sys, io, contextlib
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from sklearn.model_selection import KFold
_ARGS = sys.argv[1:]; sys.argv = ['x']; sys.path.insert(0, '.')
with contextlib.redirect_stdout(io.StringIO()):
    import audit2 as A

KPC = 3.085677581491367e19; c = 299792458.0; LN10 = math.log(10)
n, ng, idx, names = A.n, A.ng, A.idx, A.names
r0, V0, eV0 = A.r, A.V, A.eV
VG = np.empty(n); VD = np.empty(n); VB = np.empty(n)
# reconstruir componentes na mesma ordem de A.rows
comp = {}
for line in open('data/SPARC_Table2.txt'):
    t = line.split()
    if len(t) == 10 and t[0] in A.meta:
        try: comp[(t[0], float(t[2]))] = list(map(float, t[1:]))
        except ValueError: pass
for j, row in enumerate(A.rows):
    v = comp[(row[0], row[1])]; VG[j], VD[j], VB[j] = v[4], v[5], v[6]
gal = [np.where(idx == i)[0] for i in range(ng)]
M = lambda key: np.array([A.meta[nm][key] for nm in names])
D0, eD, i0, ei = M('D'), M('eD'), M('inc'), np.maximum(M('einc'), 1.0)
Rd, Lg, SBd = A.Rd, A.Lg, A.SBd
SIG_INT = float(_ARGS[0]) if _ARGS else 0.0

def rar(a, g): return g / (-np.expm1(-np.sqrt(g / a)))

def observables(i, th, pts=None):
    """th = (log10 Υd, log10 Υb, D, inc) → (gb, y, sigma, r) para os pontos da galáxia i"""
    w = gal[i] if pts is None else pts
    lud, lub, D, inc = th
    fD = D / D0[i]; fi = math.sin(math.radians(i0[i])) / math.sin(math.radians(inc))
    r = r0[w] * fD; V = V0[w] * fi
    vb2 = fD * (VG[w] * np.abs(VG[w]) + 10**lud * VD[w]**2 + 10**lub * VB[w]**2)
    gb = np.maximum(vb2, 1e-6) * 1e6 / (r * KPC)
    y = np.log10(V**2 * 1e6 / (r * KPC))
    sig = np.sqrt((2 / LN10 * eV0[w] / V0[w])**2 + SIG_INT**2)
    return gb, y, sig, r, vb2, w

def prior(i, th):
    lud, lub, D, inc = th
    return ((lud - math.log10(0.5)) / 0.1)**2 + ((lub - math.log10(0.7)) / 0.1)**2 + ((D - D0[i]) / eD[i])**2 + ((inc - i0[i]) / ei[i])**2

def features(model, gb, r, vb2, w):
    if model == 'RAR': return None
    if model == 'RAR+K': return np.zeros((len(w), 1))
    if model == 'log r': return np.log10(r)[:, None]
    if model == 'log(eV/V)': return np.log10(np.maximum(eV0[w], .5) / V0[w])[:, None]
    if model == 'livre (C,G,L)':
        x = r / (Rd[w] * r / r0[w]); cf = np.maximum(1 - (1 + x) * np.exp(-x), 1e-8)
        return np.column_stack([np.log10(SBd[w] * cf / x**2), np.log10(4 * np.maximum(vb2, 1e-6) * 1e6 / c**2), np.log10(Lg[w] * 1e9 / (4 * math.pi * r0[w]**2))])
    raise ValueError(model)

def gal_obj(i, th, model, glob, center):
    gb, y, sig, r, vb2, w = observables(i, th)
    F = features(model, gb, r, vb2, w)
    if F is not None: F = F - center
    p = pred(model, glob, gb, F)
    return np.sum(((p - y) / sig)**2) + prior(i, th)

def pred(model, glob, gb, F):
    a = 10**glob[0]; base = rar(a, gb)
    if F is None: return np.log10(base)
    z = np.clip(glob[1] + F @ glob[2:], -6, 6)
    return np.log10(gb + (base - gb) * 10**z)

def bounds_i(i):
    return [(math.log10(.5) - .3, math.log10(.5) + .3), (math.log10(.7) - .3, math.log10(.7) + .3),
            (max(D0[i] - 3 * eD[i], 0.05 * D0[i]), D0[i] + 3 * eD[i]), (max(i0[i] - 3 * ei[i], 10), min(i0[i] + 3 * ei[i], 90))]

def fit_gal(i, model, glob, center, th0):
    res = minimize(lambda t: gal_obj(i, t, model, glob, center), th0, method='L-BFGS-B', bounds=bounds_i(i), options={'maxiter': 200})
    return res.x, res.fun

def th_default(i): return np.array([math.log10(.5), math.log10(.7), D0[i], i0[i]])

def fit_global(model, gals, nuis, center, glob0):
    data = [observables(i, nuis[i]) for i in gals]
    Fs = [features(model, d[0], d[3], d[4], d[5]) for d in data]
    def obj(g):
        s = 0.0
        for d, F in zip(data, Fs):
            p = pred(model, g, d[0], None if F is None else F - center)
            s += np.sum(((p - d[1]) / d[2])**2)
        return s
    k = 1 if model == 'RAR' else 2 + Fs[0].shape[1]
    bnds = [(-12, -9)] + ([] if k == 1 else [(-2, 2)] + [(-3, 3)] * (k - 2))
    return minimize(obj, glob0[:k], method='L-BFGS-B', bounds=bnds, options={'maxiter': 400, 'ftol': 1e-12}).x

def train(model, gals, nuisance=True, iters=4):
    nuis = {i: th_default(i) for i in gals}
    # centro dos proxies = mediana no treino (com incômodos nominais)
    F_all = [features(model, *[observables(i, nuis[i])[q] for q in (0, 3, 4, 5)]) for i in gals]
    center = 0 if F_all[0] is None else np.median(np.vstack(F_all), axis=0)
    glob = np.array([math.log10(1.2e-10), 0, 0, 0, 0, 0], float)
    for it in range(iters if nuisance else 1):
        glob = fit_global(model, gals, nuis, center, np.concatenate([glob, np.zeros(6)]))
        if not nuisance: break
        for i in gals: nuis[i], _ = fit_gal(i, model, glob, center, nuis[i])
    return glob, center

def evaluate(models, nuisance=True, seed=42):
    folds = list(KFold(5, shuffle=True, random_state=seed).split(np.arange(ng)))
    out = {m: dict(chi2=0.0, per=np.zeros(ng), glob=[]) for m in models}
    resid = {m: np.full(n, np.nan) for m in models}
    for itr, ite in folds:
        for m in models:
            glob, center = train(m, list(itr), nuisance)
            out[m]['glob'].append(glob.tolist())
            for i in ite:
                th = fit_gal(i, m, glob, center, th_default(i))[0] if nuisance else th_default(i)
                gb, y, sig, r, vb2, w = observables(i, th)
                F = features(m, gb, r, vb2, w); p = pred(m, glob, gb, None if F is None else F - center)
                out[m]['chi2'] += np.sum(((p - y) / sig)**2) + (prior(i, th) if nuisance else 0)
                out[m]['per'][i] = np.mean((p - y)**2); resid[m][w] = y - p
    return out, resid

def report(title, out):
    b = out['RAR']
    print(f'\n== {title} ==')
    for m, o in sorted(out.items(), key=lambda t: t[1]['chi2']):
        rm = math.sqrt(o['per'].mean()); g = np.mean(o['glob'], axis=0).round(4).tolist()
        print(f"  {m:16s} χ²_teste {o['chi2']:10.1f}  Δχ² vs RAR {o['chi2'] - b['chi2']:+9.1f}   RMSE/gal {rm:.4f} dex   globais médios {g}")
    return {m: dict(chi2=o['chi2'], dchi2=o['chi2'] - b['chi2'], rmse_gal=math.sqrt(o['per'].mean()), globais=o['glob'],
                    galaxias_melhoram=int(np.sum(o['per'] < b['per']))) for m, o in out.items()}

def slope_diag(resid, title):
    """tendência do resíduo de teste com log r: inclinação ponderada por galáxia, interno vs externo"""
    e = resid; lr = np.log10(A.r); inner = A.rank < 2
    wt = 1 / np.bincount(idx, minlength=ng)[idx]
    X = np.column_stack([np.ones(n), lr]); W = np.sqrt(wt)
    beta = np.linalg.lstsq(X * W[:, None], e * W, rcond=None)[0]
    print(f'  [{title}] resíduo de teste: inclinação d(res)/d(log r) = {beta[1]:+.4f} dex/dex;  média 2 pts internos {np.mean(e[inner]):+.4f}, demais {np.mean(e[~inner]):+.4f}')
    return dict(inclinacao=float(beta[1]), media_internos=float(np.mean(e[inner])), media_demais=float(np.mean(e[~inner])))

if __name__ == '__main__':
    models = ['RAR', 'RAR+K', 'log r', 'livre (C,G,L)', 'log(eV/V)']
    R = {}
    for nuis in [False, True]:
        tag = 'com_incomodos' if nuis else 'sem_incomodos'
        out, resid = evaluate(models, nuisance=nuis)
        R[tag] = report(f"χ² {'COM' if nuis else 'SEM'} incômodos por galáxia (Υd, Υb, D, i)  σ_int={SIG_INT}", out)
        R[tag + '_diag'] = {m: slope_diag(resid[m], m) for m in ['RAR', 'RAR+K', 'log r']}
    json.dump(R, open(f'results/auditoria_fase1c_resultados_sigint{SIG_INT}.json', 'w'), indent=1, ensure_ascii=False)
