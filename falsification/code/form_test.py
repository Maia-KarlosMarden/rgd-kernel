import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))  # paths relativos a v5/code
"""Item 2 — existe ALGUMA lei local g_obs = F(g_b, α_b) que sirva a galáxias e clusters?
Teste livre de forma funcional: comparar SPARC e BCGs CLASH no MESMO (g_b, α_b).
α_b = deflexão interna tangencial calculada só com bárions (lens_test.alpha_tangential).
BCG: duas versões de α_b — (i) BCG + ICM do cluster (profundidade total, 'ambiental'); (ii) só a BCG (isolada).
"""
import sys, io, contextlib, math, json, runpy
import numpy as np
sys.path.insert(0, '.')
with contextlib.redirect_stdout(io.StringIO()):
    import lens_test as LT
    CT = runpy.run_path('clusters_test.py')
B = LT.B; A = B.A; KPC, G, MSUN = LT.KPC, LT.G, LT.MSUN
# SPARC (nominal Υd=0.5, Υb=0.7)
a_s = np.maximum(LT.aG + .5 * LT.aD + .7 * LT.aB, 1e-12)
gb_s, go_s = A.gb, A.go
Q = np.array([A.meta[nm]['Q'] for nm in A.names])[A.idx]
eps = (2 / math.log(10)) * A.eV / A.V
# BCG
gb_c, gt_c, r_c = CT['gbar'], CT['gtot'], CT['r']; Re = CT['Re']
Ms_tot = np.array([t[3] for t in CT['T']]) * 1e11 * MSUN
Mg_c = 0.1 * np.array([t[4] or 0.0 for t in CT['T']]) * 1e11 * MSUN
sig_c = CT['sig_tot']
rs = np.logspace(math.log10(1 * KPC), math.log10(5000 * KPC), 3000)
Mt500, r200 = None, None
def nfw_M(rr, M200=1.34e15 * MSUN, c200=4.0, z=0.35):
    H = 70e3 / 3.0857e22 * math.sqrt(0.3 * (1 + z)**3 + 0.7); rho = 3 * H**2 / (8 * math.pi * G)
    R200 = (3 * M200 / (4 * math.pi * 200 * rho))**(1 / 3); x = rr / (R200 / c200); f = lambda u: np.log(1 + u) - u / (1 + u)
    return M200 * f(x) / f(c200), R200
Mt500, r200 = nfw_M(500 * KPC); rc = 100 * KPC; shp = lambda rr: rr - rc * np.arctan(rr / rc)
Mg_icm = Mt500 / 8 * shp(np.minimum(rs, r200)) / shp(500 * KPC)
a_env, a_iso = [], []
for k in range(len(gb_c)):
    rh = 0.551 * Re[k]; Mst = Ms_tot[k] * rs**2 / (rs + rh)**2
    a_env.append(LT.alpha_tangential(np.array([r_c[k]]), rs, G * (1.12 * Mg_icm + Mst) / rs**2)[0])
    a_iso.append(LT.alpha_tangential(np.array([r_c[k]]), rs, G * Mst / rs**2)[0])
a_env, a_iso = np.array(a_env), np.array(a_iso)

L = np.log10
print(f'SPARC: α_b de {np.percentile(a_s,1):.1e} a {a_s.max():.1e}; g_b de {gb_s.min():.1e} a {gb_s.max():.1e}')
print(f'BCG  : α_b(BCG+ICM) de {a_env.min():.1e} a {a_env.max():.1e};  α_b(só BCG) de {a_iso.min():.1e} a {a_iso.max():.1e};  g_b de {gb_c.min():.1e} a {gb_c.max():.1e}')

def compare(a_cl, label, dlg=0.25, dla=0.25, qsel=None):
    """para cada BCG, pontos SPARC dentro de ±dlg dex em g_b e ±dla dex em α_b; compara log(g_obs/g_b)"""
    rows = []; diffs = []; wts = []
    for k in range(len(gb_c)):
        m = (np.abs(L(gb_s) - L(gb_c[k])) < dlg) & (np.abs(L(a_s) - L(a_cl[k])) < dla)
        if qsel is not None: m &= qsel
        ngal = len(np.unique(A.idx[m]))
        if ngal < 2: continue
        # mediana por galáxia, depois mediana entre galáxias (peso igual por galáxia)
        per = [np.median(L(go_s[m & (A.idx == g)] / gb_s[m & (A.idx == g)])) for g in np.unique(A.idx[m])]
        sp = float(np.median(per)); spread = float(np.std(per))
        bc = float(L(gt_c[k] / gb_c[k]))
        rows.append((CT['names'][k], L(gb_c[k]), L(a_cl[k]), bc, sp, spread, ngal))
        diffs.append(bc - sp); wts.append(1 / (sig_c[k]**2 + spread**2 / ngal + 0.057**2))
    print(f'\n== {label}: janela ±{dlg} dex em g_b, ±{dla} dex em α_b ==')
    if not rows: print('  nenhuma BCG tem vizinhos SPARC (sem sobreposição)'); return dict(n=0)
    print(f"  {'cluster':20s} {'log g_b':>7s} {'log α_b':>7s} {'BCG log(gt/gb)':>14s} {'SPARC mesmo ponto':>17s} {'Δ':>6s} {'n gal':>5s}")
    for r_ in rows: print(f'  {r_[0]:20s} {r_[1]:7.2f} {r_[2]:7.2f} {r_[3]:14.3f} {r_[4]:11.3f}±{r_[5]:.2f} {r_[3]-r_[4]:+6.2f} {r_[6]:5d}')
    d = np.array(diffs); w = np.array(wts); mu = np.sum(w * d) / w.sum(); se = 1 / math.sqrt(w.sum())
    print(f'  → diferença média ponderada BCG − SPARC = {mu:+.3f} ± {se:.3f} dex  ({mu/se:.1f}σ),  fator {10**mu:.2f}× em g_obs;  {len(d)} BCGs com vizinhos')
    return dict(n=len(d), delta_dex=float(mu), erro=float(se), sigma=float(mu / se), linhas=[list(map(lambda x: x if isinstance(x, str) else float(x), r_)) for r_ in rows])

R = {}
R['iso_025'] = compare(a_iso, 'α_b só da BCG (lei local sem ambiente)')
R['iso_025_Q1'] = compare(a_iso, 'idem, só SPARC Q=1', qsel=(Q == 1))
R['iso_040'] = compare(a_iso, 'idem, janela larga', dlg=0.4, dla=0.4)
R['env_025'] = compare(a_env, 'α_b BCG + ICM (profundidade total do ambiente)')
R['env_040'] = compare(a_env, 'idem, janela larga', dlg=0.4, dla=0.4)
R['env_060'] = compare(a_env, 'idem, janela muito larga', dlg=0.4, dla=0.6)
# só g_b (lei tipo RAR, sem α): mesma comparação ignorando α
R['so_gb'] = compare(np.full(len(gb_c), np.median(a_s)), 'só g_b (α ignorado: teste da RAR universal)', dla=9)
# janela de transição exigida: maior α de galáxia vs α ambiental das BCGs
amax_gal = np.percentile(a_s, 99.5)
R['janela'] = dict(alpha_gal_995=float(amax_gal), alpha_env_min=float(a_env.min()), alpha_env_med=float(np.median(a_env)),
                   dlog=float(L(np.median(a_env)) - L(amax_gal)), inclinacao_min=float(math.log10(17) / max(L(np.median(a_env)) - L(amax_gal), 1e-3)))
print(f"\nJanela para uma transição em α_b: maior α de galáxia (p99.5) {amax_gal:.2e} → α ambiental mediano das BCGs {np.median(a_env):.2e}"
      f" ({R['janela']['dlog']:.2f} dex);  d log a0/d log α ≥ {R['janela']['inclinacao_min']:.1f} nessa janela")
json.dump(R, open('results/form_resultados.json', 'w'), indent=1, ensure_ascii=False)
