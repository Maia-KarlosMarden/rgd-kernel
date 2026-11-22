import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))  # paths relativos a v5/code
"""Trilha do lenteamento — modelos que imitam 'gravidade induzida pela luz que percorre o espaço-tempo curvo'.

Ideia testada (do usuário): massa curva o espaço-tempo; a luz que o preenche, ao ser defletida,
induziria um campo gravitacional adicional acoplado ao original (análogo à indução E↔B).

Duas versões fenomenológicas, ambas sem usar a massa faltante como entrada:
  (A) QUANTIDADE: g_ind ∝ (energia radiativa) × (curvatura)   → já testada: K ∝ 1/L nos clusters (falha)
  (B) GEOMETRIA: a presença de luz basta; o efeito depende só da deflexão interna α_b que um fóton
      emitido dentro do sistema sofre ao sair, calculada SÓ com a massa bariônica (GR de campo fraco):
          α_b(r) = (2 r / c²) ∫_0^∞ g_b(r cosh t) dt        (fóton emitido tangencialmente em r, até o infinito)
      e a escala de aceleração passa a ser  a_eff = a0 · (α_b / α_ref)^n.
      n = 0 recupera a RAR (sem indução geométrica). A hipótese exige n ≠ 0 e o MESMO n em galáxias e clusters.
      B1: α_b local (no raio de cada ponto)     B2: α_b característico da galáxia (mediana)

Teste de falsificação: medir n nas galáxias (SPARC, χ² + incômodos, CV por galáxia) e usar esse n para
PREVER a escala dos clusters (g‡ = 2.02e-9 m/s², Tian+2020) a partir de α_b dos clusters.
"""
import sys, math, json, io, contextlib
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from sklearn.model_selection import KFold
sys.argv = ['x', '0.057']; sys.path.insert(0, '.')
with contextlib.redirect_stdout(io.StringIO()):
    import audit3 as B
A = B.A
KPC, c, G, MSUN = B.KPC, B.c, 6.67430e-11, 1.98847e30
n, ng, idx, gal = B.n, B.ng, B.idx, B.gal
TG = np.linspace(0, 14, 700)

def alpha_tangential(rq, rs, gs):
    """α para fótons emitidos tangencialmente em cada rq (m), perfil g(rs) (m/s²); cauda pontual além de rs[-1]"""
    out = np.empty(len(rq))
    for j, r0 in enumerate(rq):
        s = r0 * np.cosh(TG)
        g = np.interp(s, rs, gs)
        far = s > rs[-1]; g[far] = gs[-1] * (rs[-1] / s[far])**2
        out[j] = 2 * r0 / c**2 * np.trapezoid(g, TG)
    return out

# deflexões por componente (Υ=1, D=D0) — α é linear nas componentes e ∝ D (g_b independe de D, r ∝ D)
aG = np.zeros(n); aD = np.zeros(n); aB = np.zeros(n)
for i in range(ng):
    w = gal[i]; o = np.argsort(B.r0[w]); w = w[o]
    rs = B.r0[w] * KPC
    for comp, V in ((aG, B.VG), (aD, B.VD), (aB, B.VB)):
        gs = V[w] * np.abs(V[w]) * 1e6 / rs
        comp[w] = alpha_tangential(rs, rs, gs)

def alpha_of(i, th):
    lud, lub, D, inc = th
    a = (D / B.D0[i]) * (aG[gal[i]] + 10**lud * aD[gal[i]] + 10**lub * aB[gal[i]])
    return np.maximum(a, 1e-12)

def feat(model, i, th):
    if model == 'RAR': return None
    a = alpha_of(i, th)
    if model == 'B1 α local': return np.log10(a)
    if model == 'B2 α galáxia': return np.full(len(a), np.median(np.log10(a)))

def pred(glob, gb, F, center):
    la = glob[0] if F is None else glob[0] + glob[1] * (F - center)
    return np.log10(B.rar(10**la, gb))

def gal_obj(i, th, model, glob, center):
    gb, y, sig, r, vb2, w = B.observables(i, th)
    return np.sum(((pred(glob, gb, feat(model, i, th), center) - y) / sig)**2) + B.prior(i, th)

def fit_gal(i, model, glob, center, th0):
    r = minimize(lambda t: gal_obj(i, t, model, glob, center), th0, method='L-BFGS-B', bounds=B.bounds_i(i), options={'maxiter': 200})
    return r.x

def fit_global(model, gals, nuis, center, g0):
    data = [(B.observables(i, nuis[i]), feat(model, i, nuis[i])) for i in gals]
    def obj(g):
        return sum(np.sum(((pred(g, d[0], F, center) - d[1]) / d[2])**2) for d, F in data)
    k = 1 if model == 'RAR' else 2
    return minimize(obj, g0[:k], method='L-BFGS-B', bounds=[(-12, -8), (-3, 3)][:k], options={'maxiter': 300, 'ftol': 1e-12}).x

def train(model, gals, iters=4):
    nuis = {i: B.th_default(i) for i in gals}
    center = 0.0 if model == 'RAR' else float(np.median(np.concatenate([feat(model, i, nuis[i]) for i in gals])))
    glob = np.array([math.log10(1.2e-10), 0.0])
    for _ in range(iters):
        glob = fit_global(model, gals, nuis, center, np.concatenate([glob, [0.0]]))
        for i in gals: nuis[i] = fit_gal(i, model, glob, center, nuis[i])
    return glob, center, nuis

if __name__ == '__main__':
    R = {}
    la_nom = np.log10(np.maximum(B.r0 * 0 + (aG + .5 * aD + .7 * aB), 1e-12))
    gal_med = np.array([np.median(la_nom[gal[i]]) for i in range(ng)])
    print(f'α_b interno (bárions, nominal): galáxias de {10**gal_med.min():.1e} a {10**gal_med.max():.1e} rad  '
          f'(fator {10**(gal_med.max()-gal_med.min()):.0f});  pontos de {10**np.percentile(la_nom,1):.1e} a {10**np.percentile(la_nom,99):.1e}')
    R['alpha_galaxias'] = dict(min=float(10**gal_med.min()), max=float(10**gal_med.max()), mediana=float(10**np.median(gal_med)))

    # --- 1. CV por galáxia, χ² + incômodos
    models = ['RAR', 'B1 α local', 'B2 α galáxia']
    folds = list(KFold(5, shuffle=True, random_state=42).split(np.arange(ng)))
    chi = {m: 0.0 for m in models}; per = {m: np.zeros(ng) for m in models}; globs = {m: [] for m in models}
    for itr, ite in folds:
        for m in models:
            glob, center, _ = train(m, list(itr)); globs[m].append([*glob.tolist(), center])
            for i in ite:
                th = fit_gal(i, m, glob, center, B.th_default(i))
                gb, y, sig, r, vb2, w = B.observables(i, th)
                p = pred(glob, gb, feat(m, i, th), center)
                chi[m] += np.sum(((p - y) / sig)**2) + B.prior(i, th); per[m][i] = np.mean((p - y)**2)
    print('\n== 1. SPARC: χ² de teste (5-fold por galáxia, incômodos Υd, Υb, D, i; σ_int=0.057) ==')
    R['cv'] = {}
    for m in models:
        ns = [g[1] for g in globs[m]] if m != 'RAR' else [0] * 5
        R['cv'][m] = dict(chi2=chi[m], dchi2=chi[m] - chi['RAR'], rmse=math.sqrt(per[m].mean()), n_folds=ns, a0_folds=[10**g[0] for g in globs[m]])
        print(f"  {m:14s} χ² {chi[m]:8.1f}  Δχ² {chi[m]-chi['RAR']:+7.1f}   n por fold {np.round(ns,3).tolist()}   a0 {np.mean([10**g[0] for g in globs[m]]):.2e}")

    # --- 2. n no conjunto completo + perfil de verossimilhança em n (B2, que é o que importa para galáxia→cluster)
    allg = list(range(ng))
    globF, centF, nuisF = train('B2 α galáxia', allg)
    def prof_chi(nfix):
        g = np.array([globF[0], nfix])
        def ob(la):
            g[0] = la; return sum(gal_obj(i, nuisF[i], 'B2 α galáxia', g, centF) for i in allg)
        r = minimize_scalar(ob, bounds=(globF[0] - 1.5, globF[0] + 1.5), method='bounded')
        return r.fun
    grid = np.round(np.arange(-0.3, 0.71, 0.05), 3); prof = np.array([prof_chi(v) for v in grid])
    dch = prof - prof.min(); ok = grid[dch < 3.84]
    print(f'\n== 2. B2 no conjunto completo: n = {globF[1]:+.3f};  IC95 (perfil, Δχ²<3.84, incômodos fixos): [{ok.min():+.2f}, {ok.max():+.2f}]')
    for v, d in zip(grid, dch): print(f'     n={v:+.2f}  Δχ²={d:8.1f}')
    R['B2_full'] = dict(n=float(globF[1]), a0_ref=float(10**globF[0]), log_alpha_ref=centF, ic95=[float(ok.min()), float(ok.max())],
                        perfil=dict(zip(map(float, grid), map(float, dch))))

    # --- 3. clusters: α_b bariônico e n exigido
    # BCG (Hernquist, Tian+2020 Tab.1) + ICM (β-modelo β=2/3, normalizado por M_gas(<500 kpc)=M_tot,NFW/8, Donahue+2014)
    # NFW médio CLASH: M200c = 1.34e15 Msun, c200c = 4.0 (Umetsu+2014); estrelas fora da BCG: +12% do gás (f_c ~ 0.1, Chiu+2018)
    def nfw_M(rr, M200=1.34e15 * MSUN, c200=4.0, z=0.35):
        H = 70e3 / 3.0857e22 * math.sqrt(0.3 * (1 + z)**3 + 0.7); rho_c = 3 * H**2 / (8 * math.pi * G)
        r200 = (3 * M200 / (4 * math.pi * 200 * rho_c))**(1 / 3); rs = r200 / c200
        f = lambda x: np.log(1 + x) - x / (1 + x)
        return M200 * f(rr / rs) / f(c200), r200
    rs_grid = np.logspace(math.log10(1 * KPC), math.log10(5000 * KPC), 3000)
    Mt500, r200 = nfw_M(500 * KPC)
    cl = {}
    import runpy
    with contextlib.redirect_stdout(io.StringIO()):
        CT = runpy.run_path('clusters_test.py')
    Re_ = CT['Re']; Mstot = np.array([t[3] for t in CT['T']]) * 1e11 * MSUN; rB = CT['r']; aeff_obs = CT['a_eff']
    for rc_kpc in (50, 100, 200):
        rc = rc_kpc * KPC
        mshape = lambda rr: rr - rc * np.arctan(rr / rc)
        Mg = Mt500 / 8 * mshape(np.minimum(rs_grid, r200)) / mshape(500 * KPC)
        aB_bcg, aB_100, aB_400 = [], [], []
        for k in range(len(Mstot)):
            rh = 0.551 * Re_[k]
            Mb = 1.12 * Mg + Mstot[k] * rs_grid**2 / (rs_grid + rh)**2
            gs = G * Mb / rs_grid**2
            a3 = alpha_tangential(np.array([rB[k], 100 * KPC, 400 * KPC]), rs_grid, gs)
            aB_bcg.append(a3[0]); aB_100.append(a3[1]); aB_400.append(a3[2])
        cl[rc_kpc] = dict(alpha_bcg_med=float(np.median(aB_bcg)), alpha_100=float(np.median(aB_100)), alpha_400=float(np.median(aB_400)), alpha_bcg=aB_bcg)
    a0g = 10**globF[0]; aref = 10**centF
    a0_RAR = np.mean(R['cv']['RAR']['a0_folds'])
    print(f'\n== 3. Clusters: deflexão interna bariônica (ensemble CLASH) ==')
    print(f'  galáxias: α_b mediano {10**np.median(gal_med):.2e} rad;   a0 galáctico (RAR, χ²+incômodos) = {a0_RAR:.2e}')
    R['clusters'] = {}
    for rc_kpc, d in cl.items():
        ratio = d['alpha_100'] / 10**np.median(gal_med)
        n_req = math.log10(2.02e-9 / a0_RAR) / math.log10(ratio)
        pred_g = a0g * (d['alpha_100'] / aref)**globF[1]
        print(f"  r_c={rc_kpc:3d} kpc: α_b(BCG)={d['alpha_bcg_med']:.2e}  α_b(100 kpc)={d['alpha_100']:.2e}  α_b(400 kpc)={d['alpha_400']:.2e}"
              f"  | α_cl/α_gal={ratio:5.1f}  → n EXIGIDO p/ g‡=2.02e-9: {n_req:+.2f}   | g‡ PREVISTO com n_galáxias: {pred_g:.2e}")
        R['clusters'][rc_kpc] = dict(alpha_bcg=d['alpha_bcg_med'], alpha_100=d['alpha_100'], alpha_400=d['alpha_400'], razao=ratio, n_exigido=n_req, g_previsto_n_gal=pred_g)
    # dentro dos clusters: a_eff por cluster vs α_b(BCG) por cluster
    from scipy import stats
    rs_, ps_ = stats.spearmanr(np.log10(cl[100]['alpha_bcg']), np.log10(aeff_obs))
    print(f'  entre clusters (escala BCG): Spearman(log α_b, log a_eff) = {rs_:+.2f} (p={ps_:.2f})')
    R['clusters']['spearman_bcg'] = dict(rho=float(rs_), p=float(ps_))
    json.dump(R, open('results/lens_resultados.json', 'w'), indent=1, ensure_ascii=False, default=float)
