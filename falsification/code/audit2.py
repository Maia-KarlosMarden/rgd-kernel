import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))  # paths relativos a v5/code
"""Fase 1b — origem do efeito radial no resíduo da RAR (SPARC).
Reconstrói a amostra do zero com parser corrigido (175 galáxias no catálogo),
mesmos filtros do pacote (Q<=2, >=5 pontos, g_b>0), mesma perda (log10, peso igual por galáxia),
mesmo modelo de correção multiplicativa sobre o resíduo positivo da RAR, 5-fold por galáxia.
Pergunta: o ganho de 'log r' é físico (raio) ou sistemático (distância, resolução angular,
raios internos, inclinação, erro de velocidade)?
"""
import math, json, sys
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from sklearn.model_selection import KFold

KPC = 3.085677581491367e19; c = 299792458.0
DATA = sys.argv[1] if len(sys.argv) > 1 else 'data'

# ---------- parser corrigido ----------
meta = {}
for line in open(f'{DATA}/SPARC_Lelli2016c.mrt'):
    t = line.split()
    if len(t) < 19: continue
    try: v = list(map(float, t[1:18]))
    except ValueError: continue
    T, D, eD, fD, inc, einc, L, eL, Re, SBe, Rd, SBd, MHI, RHI, Vf, eVf, Q = v
    meta[t[0]] = dict(T=T, D=D, eD=eD, fD=fD, inc=inc, einc=einc, L=L, Rd=Rd, SBd=SBd, MHI=MHI, Vf=Vf, Q=Q)
print('catálogo:', len(meta), 'galáxias')

rows = []
for line in open(f'{DATA}/SPARC_Table2.txt'):
    t = line.split()
    if len(t) != 10 or t[0] not in meta: continue
    try: D, r, v, ev, vg, vd, vbu, sbd, sbb = map(float, t[1:])
    except ValueError: continue
    m = meta[t[0]]
    if m['Q'] > 2 or r <= 0 or v <= 0 or m['L'] <= 0 or m['Rd'] <= 0: continue
    vb2 = vg * abs(vg) + .5 * vd**2 + .7 * vbu**2
    if vb2 <= 0: continue
    rows.append((t[0], r, v, ev, vb2))
from collections import Counter
cnt = Counter(x[0] for x in rows)
rows = [x for x in rows if cnt[x[0]] >= 5]
names = sorted(set(x[0] for x in rows)); ng = len(names)
idx = np.array([names.index(x[0]) for x in rows])
r = np.array([x[1] for x in rows]); V = np.array([x[2] for x in rows]); eV = np.array([x[3] for x in rows])
vb2 = np.array([x[4] for x in rows])
gb = vb2 * 1e6 / (r * KPC); go = V**2 * 1e6 / (r * KPC); y = np.log10(go); n = len(r)
G = lambda k: np.array([meta[nm][k] for nm in names])[idx]
Dist, eD, fD, inc, Rd, Lg, SBd = G('D'), G('eD'), G('fD'), G('inc'), G('Rd'), G('L'), G('SBd')
r_arcsec = r / (Dist * 1e3) * 206265.0
rank = np.zeros(n, int)
for i in range(ng):
    w = np.where(idx == i)[0]; rank[w[np.argsort(r[w])]] = np.arange(len(w))
rnorm = np.zeros(n)  # r / r_max da galáxia
for i in range(ng):
    w = idx == i; rnorm[w] = r[w] / r[w].max()
x = r / Rd; cf = np.maximum(1 - (1 + x) * np.exp(-x), 1e-8)
logC = np.log10(SBd * cf / x**2); logGp = np.log10(4 * vb2 * 1e6 / c**2); logLp = np.log10(Lg * 1e9 / (4 * math.pi * r**2))
print(f'amostra: {ng} galáxias, {n} pontos, g_obs<=g_b: {int(np.sum(go <= gb))}')

def rar(a, g): return g / (-np.expm1(-np.sqrt(g / a)))

def evaluate(feats, mask_pts=None, seed=42, return_coef=False):
    """feats: dict nome -> matriz (n,k) ou None (RAR). Avalia só nos pontos de mask_pts."""
    M = np.ones(n, bool) if mask_pts is None else mask_pts
    gal = np.unique(idx[M]); cts = np.bincount(idx[M], minlength=ng); wt = np.zeros(n); wt[M] = 1 / cts[idx[M]]
    folds = list(KFold(5, shuffle=True, random_state=seed).split(gal))
    pred = {k: np.full(n, np.nan) for k in feats}; coefs = {k: [] for k in feats}
    for itr, ite in folds:
        tr = M & np.isin(idx, gal[itr]); te = M & np.isin(idx, gal[ite])
        wm = lambda e, m: np.sum(wt[m] * e * e) / wt[m].sum()
        a = 10**minimize_scalar(lambda z: wm(np.log10(rar(10**z, gb[tr])) - y[tr], tr), bounds=(-12, -9), method='bounded', options={'xatol': 1e-8}).x
        base = rar(a, gb) - gb
        for k, F in feats.items():
            if F is None: pred[k][te] = np.log10(rar(a, gb[te])); continue
            F = F.reshape(n, -1); Fc = F - np.median(F[tr], axis=0)
            def obj(t):
                e = np.log10(gb[tr] + base[tr] * 10**np.clip(t[0] + Fc[tr] @ t[1:], -6, 6)) - y[tr]
                return wm(e, tr) + 0.002 * np.sum(t[1:]**2)
            th = minimize(obj, np.zeros(Fc.shape[1] + 1), method='L-BFGS-B', bounds=[(-2, 2)] + [(-3, 3)] * Fc.shape[1], options={'maxiter': 300, 'ftol': 1e-11}).x
            coefs[k].append(th[1:].tolist())
            pred[k][te] = np.log10(gb[te] + base[te] * 10**np.clip(th[0] + Fc[te] @ th[1:], -6, 6))
    def score(p):
        e = p - y; return math.sqrt(np.mean([np.mean(e[M & (idx == i)]**2) for i in gal]))
    s = {k: score(p) for k, p in pred.items()}
    b = s['RAR']; out = {k: dict(rmse=v, ganho=100 * (1 - v / b), coef=np.mean(coefs[k], axis=0).round(4).tolist() if coefs[k] else None) for k, v in s.items()}
    out['_n'] = dict(galaxias=int(len(gal)), pontos=int(M.sum()))
    return out

def show(title, res):
    print(f"\n== {title}  [{res['_n']['galaxias']} gal, {res['_n']['pontos']} pts] ==")
    for k, v in sorted(((k, v) for k, v in res.items() if k != '_n'), key=lambda t: t[1]['rmse']):
        print(f"  {k:34s} {v['rmse']:.6f}  {v['ganho']:+6.2f}%   coef {v['coef']}")

L = np.log10
R = {}
# 1. amostra corrigida: reproduz o pacote?
base_feats = {'RAR': None, 'livre (C,G,L)': np.column_stack([logC, logGp, logLp]), 'só G_proxy': logGp, 'log r [kpc]': L(r)}
R['1_amostra_corrigida'] = evaluate(base_feats); show('1. Amostra com parser corrigido', R['1_amostra_corrigida'])

# 2. o que é 'r'? decompor log r = log r_arcsec + log D
dec = {'RAR': None, 'log r [kpc]': L(r), 'log D (galáxia)': L(Dist), 'log r_arcsec (angular)': L(r_arcsec),
       'log r_arcsec + log D': np.column_stack([L(r_arcsec), L(Dist)]), 'log r/R_d': L(r / Rd),
       'log r/r_max': L(rnorm), 'log R_d (galáxia)': L(Rd), 'log r + log D': np.column_stack([L(r), L(Dist)])}
R['2_decomposicao'] = evaluate(dec); show('2. Decomposição do raio: físico × angular × distância', R['2_decomposicao'])

# 3. raios internos (beam smearing / movimentos não circulares)
for lab, msk in [('sem os 2 pontos mais internos', rank >= 2), ('r_arcsec > 30"', r_arcsec > 30), ('r > 2 R_d', r > 2 * Rd)]:
    k = '3_' + lab; R[k] = evaluate({'RAR': None, 'log r [kpc]': L(r), 'só G_proxy': logGp, 'log D': L(Dist)}, mask_pts=msk); show('3. ' + lab, R[k])

# 4. distância: método e erro
for lab, msk in [('distâncias precisas (fD = 2,3,4,5: TRGB/Cefeidas/UMa/SN)', fD >= 2), ('distância por fluxo de Hubble (fD = 1)', fD == 1),
                 ('erro relativo D < 10%', eD / Dist < 0.10)]:
    k = '4_' + lab; R[k] = evaluate({'RAR': None, 'log r [kpc]': L(r), 'log D': L(Dist)}, mask_pts=msk); show('4. ' + lab, R[k])

# 5. inclinação e erro de velocidade
R['5_inc'] = evaluate({'RAR': None, 'log r [kpc]': L(r), 'sin i': np.sin(np.radians(inc)), 'log r + sin i': np.column_stack([L(r), np.sin(np.radians(inc))]),
                       'log(eV/V)': L(np.maximum(eV, 0.5) / V)}); show('5. Inclinação e incerteza', R['5_inc'])
R['5b_inc>50'] = evaluate({'RAR': None, 'log r [kpc]': L(r)}, mask_pts=inc > 50); show('5b. Só inclinação > 50°', R['5b_inc>50'])

# 6. estabilidade em 20 partições dos principais
st = {}
for sd in range(20):
    o = evaluate({'RAR': None, 'log r [kpc]': L(r), 'log D': L(Dist), 'log r_arcsec': L(r_arcsec), 'só G_proxy': logGp}, seed=sd)
    for k, v in o.items():
        if k != '_n': st.setdefault(k, []).append(v['ganho'])
R['6_20_particoes'] = {k: dict(media=float(np.mean(v)), dp=float(np.std(v))) for k, v in st.items()}
print('\n== 6. 20 partições ==')
for k, v in R['6_20_particoes'].items(): print(f"  {k:20s} {v['media']:+.2f}% ± {v['dp']:.2f}")

json.dump(R, open('results/auditoria_fase1b_resultados.json', 'w'), indent=1, ensure_ascii=False)
