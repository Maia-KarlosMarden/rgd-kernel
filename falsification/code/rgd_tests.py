import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))  # paths relativos a v5/code
"""Testes pré-manuscrito RGD v5.
T1 (item 5 da agenda v4, substituto do DESI): λ do kernel (n=1/2, g†=cH0/2π) ajustado em SPARC
    (a) com Υ fixo e (b) com Υd, Υb, D, i marginalizados (χ² + priors, σ_int=0.057); CV 5-fold por galáxia vs RAR.
    Previsão do artigo: λ = 0.999 se mantém e a dispersão cai.
T3 (agenda nova): n = 1/2 + A·ξ_b com A LIVRE nas galáxias → A_gal (IC por perfil de χ² em CV) vs A_b,clusters = 2.26.
T4 (item 6, ação): monotonicidade de y·K(y) (condição para QUMOND bem-posto) em função de n.
Uso: python rgd_tests.py <H0 km/s/Mpc>
"""
import sys, io, contextlib, math, json
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from sklearn.model_selection import KFold
H0_KMS = float(sys.argv[1]) if len(sys.argv) > 1 else 70.0
sys.argv = ['x', '0.057']; sys.path.insert(0, '.')
with contextlib.redirect_stdout(io.StringIO()):
    import audit3 as B
A_ = B.A; ng, gal = B.ng, B.gal
c, G, MSUN, MPC = 2.998e8, 6.674e-11, 1.989e30, 3.0857e22
H0 = H0_KMS * 1e3 / MPC; GD = c * H0 / (2 * math.pi); FCOS = 0.157
rhoc0 = 3 * H0**2 / (8 * math.pi * G)
Mb_g = np.array([(0.5 * A_.meta[nm]['L'] + 1.33 * A_.meta[nm]['MHI']) * 1e9 for nm in A_.names]) * MSUN
XIB = G * (4 * math.pi / 3 * 500 * FCOS * rhoc0)**(2 / 3) * Mb_g**(1 / 3) / GD   # ξ_b por galáxia (nominal)

def lgK(gb, lam, n):
    x = gb / GD; return np.log10(1 + lam * np.exp(-x) / (-np.expm1(-x))**n)
def pred(model, glob, gb, i):
    if model == 'RAR': return np.log10(B.rar(10**glob[0], gb))
    if model == 'RGD λ': return np.log10(gb) + lgK(gb, glob[0], 0.5)
    if model == 'RGD λ=0.999': return np.log10(gb) + lgK(gb, 0.999, 0.5)
    if model == 'v5 A livre': return np.log10(gb) + lgK(gb, 0.999, 0.5 + glob[0] * XIB[i])
    if model.startswith('v5 A='): return np.log10(gb) + lgK(gb, 0.999, 0.5 + float(model[5:]) * XIB[i])
G0 = {'RAR': ([math.log10(1.2e-10)], [(-11, -9)]), 'RGD λ': ([1.0], [(0.3, 3.0)]), 'v5 A livre': ([0.0], [(-20, 40)])}

def gal_obj(i, th, model, glob):
    gb, y, sig, r, vb2, w = B.observables(i, th)
    return np.sum(((pred(model, glob, gb, i) - y) / sig)**2) + B.prior(i, th)
def fit_gal(i, model, glob, th0):
    return minimize(lambda t: gal_obj(i, t, model, glob), th0, method='L-BFGS-B', bounds=B.bounds_i(i), options={'maxiter': 200}).x
def fit_glob(model, gals, nuis):
    if model not in G0: return np.array([])
    data = [(i, B.observables(i, nuis[i])) for i in gals]
    ob = lambda g: sum(np.sum(((pred(model, g, d[0], i) - d[1]) / d[2])**2) for i, d in data)
    x0, bd = G0[model]
    if len(x0) == 1:
        return np.array([minimize_scalar(lambda v: ob([v]), bounds=bd[0], method='bounded', options={'xatol': 1e-6}).x])
def train(model, gals, nuisance=True, iters=4):
    nuis = {i: B.th_default(i) for i in gals}; glob = np.array(G0.get(model, ([0],))[0], float)
    for _ in range(iters if nuisance else 1):
        glob = fit_glob(model, gals, nuis)
        if not nuisance: break
        for i in gals: nuis[i] = fit_gal(i, model, glob, nuis[i])
    return glob, nuis
def cv(models, nuisance=True, seed=42):
    folds = list(KFold(5, shuffle=True, random_state=seed).split(np.arange(ng)))
    out = {m: dict(chi2=0.0, per=np.zeros(ng), glob=[]) for m in models}
    for itr, ite in folds:
        for m in models:
            glob, _ = train(m, list(itr), nuisance); out[m]['glob'].append(glob.tolist())
            for i in ite:
                th = fit_gal(i, m, glob, B.th_default(i)) if nuisance else B.th_default(i)
                gb, y, sig, r, vb2, w = B.observables(i, th)
                p = pred(m, glob, gb, i)
                out[m]['chi2'] += np.sum(((p - y) / sig)**2) + (B.prior(i, th) if nuisance else 0); out[m]['per'][i] = np.mean((p - y)**2)
    return out
R = dict(H0=H0_KMS, g_dagger=GD)
L = np.log10
if __name__ == '__main__':
    print(f'H0 = {H0_KMS} km/s/Mpc → g† = cH0/2π = {GD:.3e} m/s²')
    # ---------- T1 ----------
    for nuis in (False, True):
        tag = 'marginalizado' if nuis else 'Υ fixo'
        o = cv(['RAR', 'RGD λ', 'RGD λ=0.999'], nuisance=nuis)
        lam = [g[0] for g in o['RGD λ']['glob']]; a0 = [10**g[0] for g in o['RAR']['glob']]
        Npts = sum(len(g) for g in gal)
        R[f'T1_{tag}'] = {m: dict(chi2=v['chi2'], chi2_N=v['chi2'] / Npts, dchi2_vs_RAR=v['chi2'] - o['RAR']['chi2'], rms_gal=math.sqrt(v['per'].mean()), glob=v['glob']) for m, v in o.items()}
        print(f'\n== T1 ({tag}) — CV 5-fold por galáxia ==')
        for m, v in o.items():
            print(f"  {m:14s} χ²/N {v['chi2']/Npts:.3f}  Δχ² vs RAR {v['chi2']-o['RAR']['chi2']:+8.1f}  RMS/gal {math.sqrt(v['per'].mean()):.4f} dex")
        print(f'  λ por fold: {np.round(lam,4).tolist()} (média {np.mean(lam):.4f});  a0 RAR por fold: {np.mean(a0):.3e};  λ² g† = {np.mean(lam)**2*GD:.3e}')
        R[f'T1_{tag}']['lambda_media'] = float(np.mean(lam)); R[f'T1_{tag}']['lambda_dp'] = float(np.std(lam))
    # ---------- T3 ----------
    o = cv(['v5 A livre', 'RGD λ=0.999'] + [f'v5 A={a}' for a in (2.26, 1.0)], nuisance=True)
    Agal = [g[0] for g in o['v5 A livre']['glob']]
    base = o['RGD λ=0.999']['chi2']
    R['T3'] = {m: dict(chi2=v['chi2'], dchi2_vs_n05=v['chi2'] - base, glob=v['glob']) for m, v in o.items()}
    R['T3']['A_gal_folds'] = Agal; R['T3']['xi_b_range'] = [float(XIB.min()), float(XIB.max())]
    print(f'\n== T3 — inclinação A de n = 1/2 + A·ξ_b medida nas galáxias (χ² + incômodos) ==')
    print(f'  ξ_b das galáxias {XIB.min():.4f}–{XIB.max():.4f};  A_gal por fold {np.round(Agal,2).tolist()}  (média {np.mean(Agal):+.2f} ± {np.std(Agal):.2f})')
    for m, v in o.items(): print(f"  {m:14s} Δχ² vs n=1/2: {v['chi2']-base:+8.1f}")
    # perfil de χ² em A no conjunto completo (incômodos reotimizados para cada A)
    allg = list(range(ng)); prof = {}
    for a in (-4, -2, 0, 1, 2.26, 4, 6):
        m = f'v5 A={a}'; _, nuis = train('RGD λ=0.999', allg, iters=1) if a == -4 else (None, nuis)
        prof[a] = sum(gal_obj(i, fit_gal(i, m, [], nuis[i]), m, []) for i in allg)
    pm = min(prof.values()); R['T3']['perfil_A'] = {str(k): v - pm for k, v in prof.items()}
    print('  perfil Δχ²(A), conjunto completo: ' + '  '.join(f'A={k}: {v-pm:.1f}' for k, v in prof.items()))
    # ---------- T4 ----------
    y = np.logspace(-4, 3, 4000)
    def mono(n, lam=0.999):
        f = y * (1 + lam * np.exp(-y) / (-np.expm1(-y))**n); d = np.diff(np.log(f)) / np.diff(np.log(y))
        return float(d.min()), float(y[np.argmin(d)])
    R['T4'] = {}
    print('\n== T4 — QUMOND exige y·K(y) estritamente crescente (g_obs monotônico em g_b) ==')
    for n in (0.5, 0.7, 0.8, 0.95, 1.0, 1.2, 2.5):
        dm, ym = mono(n); R['T4'][str(n)] = dict(dlog_min=dm, y_at=ym, ok=dm > 0)
        print(f'  n={n:4.2f}: min d ln(g_obs)/d ln(g_b) = {dm:+.3f} em y={ym:.2e}  → {"bem-posto" if dm > 0 else "NÃO invertível (lei de força mal definida)"}')
    json.dump(R, open(f'results/rgd_tests_H{H0_KMS:g}.json', 'w'), indent=1, ensure_ascii=False, default=float)
