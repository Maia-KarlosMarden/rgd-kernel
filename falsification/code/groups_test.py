import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))  # paths relativos a v5/code
"""Teste decisivo do v5: grupos (Sun+2009) + clusters X-COP (Eckert+2019), com variável bariônica ξ_b.
Modelos para g_obs(R_b), validação deixa-um-fora em 25 sistemas:
  v5    : n = 1/2 + A_b ξ_b                  (1 parâmetro, só bárions)
  nulo1 : M_tot = M_b / f_u                  (1 parâmetro, 'amostra justa')
  nulo2 : M_tot = M_b / [f0 (M_b/1e13)^β]    (2 parâmetros, tendência empírica de f_bar com a massa)
  v4    : n = 1/2 + A ξ(massa total)          (circular, referência)
Dados de Sun+2009 transcritos de extração parcial (tabela do PDF ilegível em partes): conferir contra o original.
Uso: python groups_test.py <H0> <f_star extra (fração de M500 em estrelas nos grupos; 0 = só ×1.09)>
"""
import sys, math, json
import numpy as np
from scipy.optimize import minimize_scalar, minimize
H0_KMS = float(sys.argv[1]) if len(sys.argv) > 1 else 70.0
FSTAR_G = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
c, G, MSUN, MPC = 2.998e8, 6.674e-11, 1.989e30, 3.0857e22
H0 = H0_KMS * 1e3 / MPC; E = lambda z: math.sqrt(.3 * (1 + z)**3 + .7)
gd = lambda z: c * H0 * E(z) / (2 * math.pi); rhoc = lambda z: 3 * (H0 * E(z))**2 / (8 * math.pi * G)
LAM, FCOS = 0.999, 0.157
def K(g, z, n): x = g / gd(z); return 1 + LAM * math.exp(-x) / (-math.expm1(-x))**n
# grupos: nome, z, M500 [1e13], fgas500, fgas2500
GR = [('NGC1550', .0124, 3.18, .097, .056), ('MKW4', .0200, 4.85, .086, .047), ('3C442A', .0263, 3.90, .068, .028),
      ('NGC5098', .0368, 2.00, .108, .048), ('A160', .0447, 7.90, .085, .043), ('UGC842', .0452, 5.60, .056, .032),
      ('RXJ1159', .0808, 8.30, .065, .042), ('A1692', .0848, 9.70, .090, .043), ('A2550', .122, 7.90, .093, .067),
      ('A2717', .0498, 12.9, .076, .053), ('RXJ1022', .0543, 8.00, .075, .038), ('AS1101', .0564, 14.1, .114, .085),
      ('A2462', .0733, 8.80, .099, .052)]
# X-COP: nome, z, M500 [1e14], fgas500, fgas200
XC = [('A1644', .0473, 3.48, .128, .126), ('A1795', .0622, 4.63, .139, .144), ('A2029', .0766, 8.65, .141, .152),
      ('A2142', .0909, 8.95, .158, .168), ('A2255', .0809, 5.26, .153, .146), ('A2319', .0557, 7.31, .189, .237),
      ('A3158', .0597, 4.26, .145, .155), ('A3266', .0589, 8.80, .132, .108), ('A644', .0704, 5.66, .132, .139),
      ('A85', .0555, 5.65, .150, .159), ('RXC1825', .0650, 4.08, .133, .155), ('ZwCl1215', .0766, 7.66, .106, .092)]
rows = []
def add(nm, z, M5, fb5, fb_other, r_ratio, m_ratio, kind):
    """fb_other e r_ratio/m_ratio: fração de bárions, R/R500 e M/M500 num segundo raio (para a inclinação local)"""
    R5 = (3 * M5 / (4 * math.pi * 500 * rhoc(z)))**(1 / 3); Mb5 = fb5 * M5
    s_b = math.log(fb_other * m_ratio / fb5) / math.log(r_ratio); s_b = min(max(s_b, 0.5), 1.8)
    s_t = math.log(m_ratio) / math.log(r_ratio)
    Rb = R5 * (fb5 / FCOS)**(1 / (3 - s_b)); Mb = Mb5 * (Rb / R5)**s_b; Mt = M5 * (Rb / R5)**s_t
    gb, gobs = G * Mb / Rb**2, G * Mt / Rb**2
    rows.append(dict(nome=nm, tipo=kind, z=z, Mb=Mb / MSUN, gb=gb, gobs=gobs, xi_b=gb / gd(z), xi_t=(G * M5 / R5**2) / gd(z), fb5=fb5, Rb_R500=Rb / R5))
for nm, z, M, f5, f25 in GR:
    fs = 1.09 * f5 + FSTAR_G
    add(nm, z, M * 1e13 * MSUN, fs, 1.09 * f25 + FSTAR_G, 0.471, 2500 / 500 * 0.471**3, 'grupo')
for nm, z, M, f5, f2 in XC:
    add(nm, z, M * 1e14 * MSUN, 1.09 * f5, 1.09 * f2, 1.53, 1.45, 'X-COP')
N = len(rows); gb = np.array([r['gb'] for r in rows]); gobs = np.array([r['gobs'] for r in rows]); zz = np.array([r['z'] for r in rows])
xb = np.array([r['xi_b'] for r in rows]); xt = np.array([r['xi_t'] for r in rows]); Mb = np.array([r['Mb'] for r in rows])
y = np.log10(gobs); L = np.log10
def pv(i, a, x): return L(gb[i] * K(gb[i], zz[i], .5 + a * x[i]))
def fit(model, m):
    idx = np.where(m)[0]
    if model in ('v5', 'v4'):
        x = xb if model == 'v5' else xt
        return [minimize_scalar(lambda a: sum((pv(i, a, x) - y[i])**2 for i in idx), bounds=(-30, 60), method='bounded').x]
    if model == 'nulo1': return [np.median(L(gb[idx] / gobs[idx]))]
    if model == 'nulo2':
        X = np.column_stack([np.ones(len(idx)), L(Mb[idx] / 1e13)]); return np.linalg.lstsq(X, L(gb[idx] / gobs[idx]), rcond=None)[0].tolist()
def predict(model, p, i):
    if model == 'v5': return pv(i, p[0], xb)
    if model == 'v4': return pv(i, p[0], xt)
    if model == 'nulo1': return L(gb[i]) - p[0]
    if model == 'nulo2': return L(gb[i]) - (p[0] + p[1] * L(Mb[i] / 1e13))
R = dict(H0=H0_KMS, fstar_grupos=FSTAR_G, linhas=rows, loo={}, ajuste_total={})
tipo = np.array([r['tipo'] for r in rows])
print(f'H0={H0_KMS}, f_star extra nos grupos={FSTAR_G}:  ξ_b grupos {xb[tipo=="grupo"].min():.3f}–{xb[tipo=="grupo"].max():.3f}, X-COP {xb[tipo=="X-COP"].min():.3f}–{xb[tipo=="X-COP"].max():.3f};'
      f'  f_bar(R500) grupos {np.median([r["fb5"] for r in rows if r["tipo"]=="grupo"]):.3f}, X-COP {np.median([r["fb5"] for r in rows if r["tipo"]=="X-COP"]):.3f}')
for model in ('v5', 'nulo1', 'nulo2', 'v4'):
    e = np.array([predict(model, fit(model, np.arange(N) != i), i) - y[i] for i in range(N)])
    p_all = fit(model, np.ones(N, bool))
    R['loo'][model] = dict(rms=float(np.sqrt(np.mean(e**2))), rms_grupos=float(np.sqrt(np.mean(e[tipo == 'grupo']**2))),
                           rms_xcop=float(np.sqrt(np.mean(e[tipo == 'X-COP']**2))), vies_grupos=float(np.mean(e[tipo == 'grupo'])), vies_xcop=float(np.mean(e[tipo == 'X-COP'])))
    R['ajuste_total'][model] = [float(v) for v in p_all]
    v = R['loo'][model]
    print(f"  LOO {model:6s} RMS {v['rms']:.3f} dex | grupos {v['rms_grupos']:.3f} (viés {v['vies_grupos']:+.3f}) | X-COP {v['rms_xcop']:.3f} (viés {v['vies_xcop']:+.3f}) | parâmetros {np.round(p_all,3).tolist()}")
# previsão cruzada: calibrar só em X-COP, prever grupos
for model in ('v5', 'nulo1', 'nulo2'):
    p = fit(model, tipo == 'X-COP'); e = np.array([predict(model, p, i) - y[i] for i in np.where(tipo == 'grupo')[0]])
    R.setdefault('xcop_para_grupos', {})[model] = dict(rms=float(np.sqrt(np.mean(e**2))), vies=float(e.mean()))
    print(f"  calibrado em X-COP → grupos: {model:6s} RMS {np.sqrt(np.mean(e**2)):.3f}  viés {e.mean():+.3f}")
json.dump(R, open(f'results/groups_resultados_H{H0_KMS:g}_fs{FSTAR_G:g}.json', 'w'), indent=1, ensure_ascii=False, default=float)
