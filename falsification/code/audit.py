import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))  # paths relativos a v5/code
"""Fase 1 — auditoria estatística do teste de expoentes assinados (SPARC).
Mesmo modelo, mesma perda e mesmos folds que fit_signed_exponents.py; acrescenta:
  1. reprodução + modelos de UMA variável (qual proxy carrega o ganho)
  2. controles sem radiação: log g_b, log r, log V_b^2, galáxia (L, SBd, Rd)
  3. forma alternativa da RAR (1 parâmetro extra, sem proxies)
  4. validação aninhada (escolha do padrão de sinais dentro do treino)
  5. bootstrap por galáxia do ΔRMSE pareado + repetição com 20 sementes de folds
  6. controle de embaralhamento de L entre galáxias de massa/tamanho similares
  7. sensibilidade a M/L (Υ_disk, Υ_bul)
"""
import json, math, itertools, sys
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from sklearn.model_selection import KFold

KPC = 3.085677581491367e19; c = 299792458.0
D = json.load(open('data/SPARC_tripla_test_results.json'))
pts = D['points']; names = D['names']; meta = D['meta']
n = len(pts); ng = len(names)
idx = np.array([names.index(p['name']) for p in pts])
r = np.array([p['r'] for p in pts]); go = np.array([p['go'] for p in pts]); y = np.log10(go)
gb0 = np.array([p['gb'] for p in pts])
cts = np.bincount(idx, minlength=ng); wt = 1 / cts[idx]

# recuperar componentes de velocidade para sensibilidade M/L
comp = {}
for line in open('data/SPARC_Table2.txt'):
    nm = line[:11].strip()
    if nm not in meta: continue
    try: vals = list(map(float, line[11:].split()))
    except ValueError: continue
    if len(vals) != 9: continue
    comp[(nm, round(vals[1], 4))] = vals
VG = np.array([comp[(p['name'], round(p['r'], 4))][4] for p in pts])
VD = np.array([comp[(p['name'], round(p['r'], 4))][5] for p in pts])
VBu = np.array([comp[(p['name'], round(p['r'], 4))][6] for p in pts])

def gbar(ud=0.5, ub=0.7):
    vb2 = VG * np.abs(VG) + ud * VD**2 + ub * VBu**2
    return vb2 * 1e6 / (r * KPC), vb2

def proxies(gb, vb2, Lgal=None):
    Rd = np.array([meta[p['name']]['Rd'] for p in pts]); SBd = np.array([meta[p['name']]['SBd'] for p in pts])
    Lg = np.array([meta[p['name']]['L'] for p in pts]) if Lgal is None else Lgal
    x = r / Rd; cf = np.maximum(1 - (1 + x) * np.exp(-x), 1e-8)
    C = SBd * cf / x**2; G = 4 * vb2 * 1e6 / c**2; L = Lg * 1e9 / (4 * math.pi * r**2)
    return np.log10(np.column_stack([C, G, L]))

def rar(a, g): return g / (-np.expm1(-np.sqrt(g / a)))
def rar_n(a, g, m):  # família de forma: m=1 recupera a RAR
    return g / (-np.expm1(-(g / a)**(m / 2)))**(1 / m)

def wmse(e, mask): w = wt[mask]; return np.sum(w * e * e) / w.sum()

def fit_a(gb, mask):
    obj = lambda z: wmse(np.log10(rar(10**z, gb[mask])) - y[mask], mask)
    return 10**minimize_scalar(obj, bounds=(-12, -9), method='bounded', options={'xatol': 1e-8}).x

def fit_corr(F, gb, base, mask, bounds, ridge=0.002):
    k = F.shape[1]
    def obj(t):
        z = np.clip(t[0] + F[mask] @ t[1:], -6, 6)
        e = np.log10(gb[mask] + base[mask] * 10**z) - y[mask]
        return wmse(e, mask) + ridge * np.sum(t[1:]**2)
    return minimize(obj, np.zeros(k + 1), method='L-BFGS-B', bounds=[(-2, 2)] + bounds,
                    options={'maxiter': 300, 'ftol': 1e-11}).x

def pat_bounds(pat):
    return [(-3, 3)] * 3 if pat == 'free' else [(0, 3) if s == '+' else (-3, 0) for s in pat]

def rmse_gal(pred_log, sel=None):
    e = pred_log - y
    per = np.array([np.mean(e[idx == i]**2) for i in range(ng)])
    return math.sqrt(per.mean()), per

PATS = ['free'] + [''.join(v) for v in itertools.product('+-', repeat=3)]

def run_cv(feature_sets, gb, vb2, seed=42, nested=False, rar_shape=False):
    """feature_sets: dict nome -> (função(train_mask)->F matriz completa centrada, bounds)"""
    folds = list(KFold(5, shuffle=True, random_state=seed).split(np.arange(ng)))
    out = {'RAR': np.full(n, np.nan)}
    for k in feature_sets: out[k] = np.full(n, np.nan)
    if nested: out['nested'] = np.full(n, np.nan); chosen = []
    if rar_shape: out['RAR_forma'] = np.full(n, np.nan)
    for itr, ite in folds:
        tr = np.isin(idx, itr); te = np.isin(idx, ite)
        a = fit_a(gb, tr); base = rar(a, gb) - gb
        out['RAR'][te] = np.log10(rar(a, gb[te]))
        for k, (mk, bnds) in feature_sets.items():
            F = mk(tr); th = fit_corr(F, gb, base, tr, bnds)
            out[k][te] = np.log10(gb[te] + base[te] * 10**np.clip(th[0] + F[te] @ th[1:], -6, 6))
        if rar_shape:
            obj = lambda t: wmse(np.log10(rar_n(10**t[0], gb[tr], t[1])) - y[tr], tr)
            t = minimize(obj, [math.log10(a), 1.0], method='L-BFGS-B', bounds=[(-12, -9), (0.2, 5)]).x
            out['RAR_forma'][te] = np.log10(rar_n(10**t[0], gb[te], t[1]))
        if nested:
            X = proxies(gb, vb2)
            inner = list(KFold(4, shuffle=True, random_state=seed + 1).split(itr))
            score = {p: 0.0 for p in ['RAR'] + PATS}
            for jtr, jte in inner:
                itr2 = np.isin(idx, itr[jtr]); ite2 = np.isin(idx, itr[jte])
                a2 = fit_a(gb, itr2); b2 = rar(a2, gb) - gb
                F2 = X - np.median(X[itr2], axis=0)
                score['RAR'] += rmse_sub(np.log10(rar(a2, gb)), ite2)
                for p in PATS:
                    th = fit_corr(F2, gb, b2, itr2, pat_bounds(p))
                    score[p] += rmse_sub(np.log10(gb + b2 * 10**np.clip(th[0] + F2 @ th[1:], -6, 6)), ite2)
            best = min(score, key=score.get); chosen.append(best)
            if best == 'RAR': out['nested'][te] = out['RAR'][te]
            else:
                F = X - np.median(X[tr], axis=0); th = fit_corr(F, gb, base, tr, pat_bounds(best))
                out['nested'][te] = np.log10(gb[te] + base[te] * 10**np.clip(th[0] + F[te] @ th[1:], -6, 6))
    if nested: out['_chosen'] = chosen
    return out

def rmse_sub(pred_log, mask):
    e = pred_log - y; gs = np.unique(idx[mask])
    return math.sqrt(np.mean([np.mean(e[(idx == i)][:]**2) for i in gs]))

def centered(Fraw):
    return lambda tr: Fraw - np.median(Fraw[tr], axis=0)

def feature_library(gb, vb2):
    X = proxies(gb, vb2)
    Rd = np.array([meta[p['name']]['Rd'] for p in pts]); SBd = np.array([meta[p['name']]['SBd'] for p in pts])
    Lg = np.array([meta[p['name']]['L'] for p in pts])
    lib = {
        'livre (C,G,L)': (centered(X), [(-3, 3)] * 3),
        '+++ (C,G,L)': (centered(X), [(0, 3)] * 3),
        'só C': (centered(X[:, [0]]), [(-3, 3)]),
        'só G_proxy': (centered(X[:, [1]]), [(-3, 3)]),
        'só L_proxy': (centered(X[:, [2]]), [(-3, 3)]),
        'CTRL só log g_b': (centered(np.log10(gb)[:, None]), [(-3, 3)]),
        'CTRL só log r': (centered(np.log10(r)[:, None]), [(-3, 3)]),
        'CTRL log g_b + log r': (centered(np.column_stack([np.log10(gb), np.log10(r)])), [(-3, 3)] * 2),
        'CTRL galáxia (L,SBd,Rd)': (centered(np.log10(np.column_stack([Lg, SBd, Rd]))), [(-3, 3)] * 3),
    }
    return lib, X

def summarize(out, title):
    base, perB = rmse_gal(out['RAR'])
    rows = []
    for k, v in out.items():
        if k.startswith('_'): continue
        s, per = rmse_gal(v)
        rows.append((k, s, 100 * (1 - s / base)))
    print(f'\n== {title} ==')
    for k, s, g in sorted(rows, key=lambda t: t[1]): print(f'  {k:28s} {s:.6f}  ganho {g:+.2f}%')
    return {k: dict(rmse=s, ganho_pct=g) for k, s, g in rows}

if __name__ == '__main__':
    R = {}
    gb, vb2 = gbar()
    assert np.allclose(gb, gb0, rtol=1e-9)
    lib, X = feature_library(gb, vb2)
    # identidade algébrica: log G = log g_b + log r + const
    print('corr(logG, log g_b + log r) =', np.corrcoef(X[:, 1], np.log10(gb) + np.log10(r))[0, 1])

    # 1+2+3+4: seed 42
    out = run_cv(lib, gb, vb2, seed=42, nested=True, rar_shape=True)
    R['seed42'] = summarize(out, 'CV 5-fold por galáxia, seed 42 (mesmos folds do pacote)')
    R['nested_chosen'] = out['_chosen']; print('  padrões escolhidos no aninhado:', out['_chosen'])

    # 5a: bootstrap pareado por galáxia (seed 42)
    rng = np.random.default_rng(0)
    _, perB = rmse_gal(out['RAR']); boots = {}
    for k in ['livre (C,G,L)', '+++ (C,G,L)', 'nested', 'só G_proxy', 'CTRL log g_b + log r', 'CTRL só log g_b', 'RAR_forma']:
        _, perK = rmse_gal(out[k]); d = []
        for _ in range(4000):
            s = rng.integers(0, ng, ng)
            d.append(100 * (1 - math.sqrt(perK[s].mean()) / math.sqrt(perB[s].mean())))
        d = np.array(d); boots[k] = dict(mediana=float(np.median(d)), ic95=[float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))], frac_pos=float(np.mean(d > 0)))
        # galáxias que melhoram
        boots[k]['galaxias_melhoram'] = int(np.sum(perK < perB))
    R['bootstrap'] = boots
    print('\n== Bootstrap por galáxia (ganho % sobre RAR, 4000 reamostragens) ==')
    for k, v in boots.items(): print(f"  {k:24s} mediana {v['mediana']:+.2f}%  IC95 [{v['ic95'][0]:+.2f}, {v['ic95'][1]:+.2f}]  P(>0)={v['frac_pos']:.3f}  galáxias que melhoram {v['galaxias_melhoram']}/{ng}")

    # 5b: 20 sementes de partição
    sub = {k: lib[k] for k in ['livre (C,G,L)', '+++ (C,G,L)', 'só G_proxy', 'CTRL log g_b + log r', 'CTRL só log g_b']}
    seeds = {}
    for sd in range(20):
        o = run_cv(sub, gb, vb2, seed=sd, rar_shape=True)
        b, _ = rmse_gal(o['RAR'])
        for k in o:
            seeds.setdefault(k, []).append(100 * (1 - rmse_gal(o[k])[0] / b))
    R['seeds'] = {k: dict(media=float(np.mean(v)), dp=float(np.std(v)), min=float(np.min(v)), max=float(np.max(v))) for k, v in seeds.items()}
    print('\n== 20 partições diferentes (ganho % sobre RAR) ==')
    for k, v in R['seeds'].items(): print(f"  {k:24s} média {v['media']:+.2f}%  dp {v['dp']:.2f}  [{v['min']:+.2f}, {v['max']:+.2f}]")

    # 6: embaralhar L entre galáxias com M_b e R_d similares
    Lg = np.array([meta[nm]['L'] for nm in names]); Rd = np.array([meta[nm]['Rd'] for nm in names])
    Mb = 0.5 * Lg + 1.33 * np.array([meta[nm]['MHI'] for nm in names])
    order = np.lexsort((np.log10(Rd), np.digitize(np.log10(Mb), np.quantile(np.log10(Mb), np.linspace(0, 1, 9)[1:-1]))))
    real = R['seed42']['livre (C,G,L)']['ganho_pct']; perm = []
    for t in range(60):
        newL = Lg.copy()
        for blk in np.array_split(order, ng // 4):
            newL[blk] = Lg[rng.permutation(blk)]
        Xp = proxies(gb, vb2, Lgal=newL[idx])
        o = run_cv({'perm': (centered(Xp), [(-3, 3)] * 3)}, gb, vb2, seed=42)
        perm.append(100 * (1 - rmse_gal(o['perm'])[0] / rmse_gal(o['RAR'])[0]))
    perm = np.array(perm)
    R['shuffle_L'] = dict(real=real, perm_media=float(perm.mean()), perm_dp=float(perm.std()), p_valor=float((1 + np.sum(perm >= real)) / (1 + len(perm))))
    print(f"\n== Embaralhamento de L (blocos de 4 por M_b e R_d, 60 perm.) ==\n  real {real:+.2f}%  permutado {perm.mean():+.2f}% ± {perm.std():.2f}  p={R['shuffle_L']['p_valor']:.3f}")

    # 7: sensibilidade M/L
    R['ML'] = {}
    print('\n== Sensibilidade a M/L (seed 42) ==')
    for ud, ub in [(0.4, 0.6), (0.5, 0.7), (0.6, 0.8)]:
        g2, v2 = gbar(ud, ub); lib2, _ = feature_library(g2, v2)
        o = run_cv({k: lib2[k] for k in ['livre (C,G,L)', 'só G_proxy', 'CTRL log g_b + log r']}, g2, v2, seed=42, rar_shape=True)
        b = rmse_gal(o['RAR'])[0]
        R['ML'][f'{ud}/{ub}'] = {k: dict(rmse=rmse_gal(v)[0], ganho=100 * (1 - rmse_gal(v)[0] / b)) for k, v in o.items()}
        print(f'  Υd={ud} Υb={ub}: RAR {b:.6f} | ' + ' | '.join(f"{k} {R['ML'][f'{ud}/{ub}'][k]['ganho']:+.2f}%" for k in o if k != 'RAR'))
    json.dump(R, open('results/auditoria_fase1_resultados.json', 'w'), indent=1, ensure_ascii=False)
