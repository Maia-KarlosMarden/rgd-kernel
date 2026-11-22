import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))  # paths relativos a v5/code
"""v5 — teste corrigido: 23 grupos Sun+2009 (Tier 1+2, Tabela 3 completa) + 12 X-COP (Eckert+2019),
contabilidade bariônica idêntica para todos: M_b = M_gas + f_*(M500)·M_tot(<r)  (estrelas seguem a massa total),
r2500 medido por grupo (M2500 = 5 M500 (r2500/r500)^3), e MOCK ΛCDM do próprio v5 (sem física nova).
Cenários de f_*:  G13 = Gonzalez+2013 3D: 0.032 (M500/1e14)^-0.48 ;  B14 = Budzynski+2014: 0.013 (≈ constante, Chabrier)
Uso: python groups_test2.py <H0>
"""
import sys, math, json
import numpy as np
from scipy.optimize import minimize_scalar
H0_KMS = float(sys.argv[1]) if len(sys.argv) > 1 else 67.4
c, G, MSUN, MPC = 2.998e8, 6.674e-11, 1.989e30, 3.0857e22
H0 = H0_KMS * 1e3 / MPC; E = lambda z: math.sqrt(.3 * (1 + z)**3 + .7)
gd = lambda z: c * H0 * E(z) / (2 * math.pi); rhoc = lambda z: 3 * (H0 * E(z))**2 / (8 * math.pi * G)
LAM, FCOS = 0.999, 0.157; L = np.log10
def K(g, z, n): x = g / gd(z); return 1 + LAM * math.exp(-x) / (-math.expm1(-x))**n
FSTAR = {'G13': lambda M: 0.032 * (M / 1e14)**-0.48, 'B14': lambda M: 0.013}
# Sun+2009 Tab.3, Tier 1+2 (M500 e f_gas,500 medidos): nome, z, r500[kpc], M500[1e13], r2500[kpc], fgas2500, fgas500
GR = [('NGC1550', .0124, 465, 3.18, 222, .056, .097), ('MKW4', .0200, 538, 4.85, 259, .047, .086),
      ('3C442A', .0263, 495, 3.90, 277, .028, .068), ('UGC5088', .0274, 364, 1.48, 163, .029, .049),
      ('NGC4104', .0282, 535, 4.85, 274, .036, .069), ('A1177', .0316, 550, 5.28, 264, .037, .060),
      ('NGC6269', .0348, 645, 8.49, 232, .044, .076), ('ESO306-017', .0358, 690, 10.3, 337, .052, .081),
      ('NGC5098', .0368, 398, 2.00, 195, .048, .108), ('A160', .0447, 626, 7.90, 286, .043, .085),
      ('UGC842', .0452, 570, 5.60, 276, .032, .056), ('A2717', .0498, 732, 12.9, 342, .053, .076),
      ('RXJ1022', .0543, 631, 8.00, 326, .038, .075), ('AS1101', .0564, 768, 14.1, 362, .085, .114),
      ('ESO351-021', .0571, 437, 3.22, 194, .036, .074), ('A3880', .0581, 799, 14.9, 309, .074, .088),
      ('A1991', .0587, 749, 13.4, 348, .066, .094), ('A1275', .0637, 592, 6.90, 202, .069, .094),
      ('A2092', .0669, 659, 8.95, 283, .050, .078), ('A2462', .0733, 646, 8.80, 327, .052, .099),
      ('RXJ1159', .0808, 630, 8.30, 273, .042, .065), ('A1692', .0848, 658, 9.70, 356, .043, .090),
      ('A2550', .122, 617, 7.90, 286, .067, .093)]
XC = [('A1644', .0473, 3.48, .128, .126), ('A1795', .0622, 4.63, .139, .144), ('A2029', .0766, 8.65, .141, .152),
      ('A2142', .0909, 8.95, .158, .168), ('A2255', .0809, 5.26, .153, .146), ('A2319', .0557, 7.31, .189, .237),
      ('A3158', .0597, 4.26, .145, .155), ('A3266', .0589, 8.80, .132, .108), ('A644', .0704, 5.66, .132, .139),
      ('A85', .0555, 5.65, .150, .159), ('RXC1825', .0650, 4.08, .133, .155), ('ZwCl1215', .0766, 7.66, .106, .092)]

def system(z, M5, fg5, fg_o, r_ratio, fs):
    """M5 [kg]; fg5, fg_o: f_gas em R500 e no 2º raio; r_ratio = R_o/R500. Retorna g_b, g_obs, ξ_b, f_bar500"""
    R5 = (3 * M5 / (4 * math.pi * 500 * rhoc(z)))**(1 / 3)
    Delta_o = 2500 if r_ratio < 1 else 200
    m_ratio = (Delta_o / 500) * r_ratio**3
    s_t = math.log(m_ratio) / math.log(r_ratio)
    fb5 = fg5 + fs; fbo = fg_o + fs
    s_b = min(max(math.log(fbo * m_ratio / fb5) / math.log(r_ratio), 0.5), 1.8)
    Rb = R5 * (fb5 / FCOS)**(1 / (3 - s_b)); Mb = fb5 * M5 * (Rb / R5)**s_b; Mt = M5 * (Rb / R5)**s_t
    gb, gobs = G * Mb / Rb**2, G * Mt / Rb**2
    return gb, gobs, gb / gd(z), fb5, Mb

def build(fsk):
    rows = []
    for nm, z, r5, M, r25, f25, f5 in GR:
        M5 = M * 1e13 * MSUN; fs = FSTAR[fsk](M * 1e13)
        rows.append((nm, 'grupo', z) + system(z, M5, f5, f25, r25 / r5, fs))
    for nm, z, M, f5, f2 in XC:
        M5 = M * 1e14 * MSUN; fs = FSTAR[fsk](M * 1e14)
        # f_gas Eckert já é gás; X-COP: estrelas pela mesma relação
        rows.append((nm, 'X-COP', z) + system(z, M5, f5, f2, 1.53, fs))
    return rows

def evaluate(rows):
    tipo = np.array([r[1] for r in rows]); zz = np.array([r[2] for r in rows]); gb = np.array([r[3] for r in rows])
    gobs = np.array([r[4] for r in rows]); xb = np.array([r[5] for r in rows]); Mb = np.array([r[7] for r in rows]); y = L(gobs); N = len(rows)
    def pv(i, a): return L(gb[i] * K(gb[i], zz[i], .5 + a * xb[i]))
    def fit(model, idx):
        if model == 'v5': return [minimize_scalar(lambda a: sum((pv(i, a) - y[i])**2 for i in idx), bounds=(-30, 60), method='bounded').x]
        if model == 'nulo1': return [np.median(L(gb[idx] / gobs[idx]))]
        X = np.column_stack([np.ones(len(idx)), L(Mb[idx] / MSUN / 1e13)]); return np.linalg.lstsq(X, L(gb[idx] / gobs[idx]), rcond=None)[0].tolist()
    def pred(model, p, i):
        if model == 'v5': return pv(i, p[0])
        if model == 'nulo1': return L(gb[i]) - p[0]
        return L(gb[i]) - (p[0] + p[1] * L(Mb[i] / MSUN / 1e13))
    out = {}
    for model in ('v5', 'nulo1', 'nulo2'):
        e = np.array([pred(model, fit(model, np.where(np.arange(N) != i)[0]), i) - y[i] for i in range(N)])
        p = fit(model, np.where(tipo == 'X-COP')[0]); ec = np.array([pred(model, p, i) - y[i] for i in np.where(tipo == 'grupo')[0]])
        out[model] = dict(loo=float(np.sqrt(np.mean(e**2))), cruzada=float(np.sqrt(np.mean(ec**2))), vies_cruzada=float(ec.mean()),
                          param_total=[float(v) for v in fit(model, np.arange(N))])
    out['_xi_b'] = dict(grupos=[float(xb[tipo == 'grupo'].min()), float(xb[tipo == 'grupo'].max())], xcop=[float(xb[tipo == 'X-COP'].min()), float(xb[tipo == 'X-COP'].max())])
    out['_fbar'] = dict(grupos=float(np.median([r[6] for r in rows if r[1] == 'grupo'])), xcop=float(np.median([r[6] for r in rows if r[1] == 'X-COP'])))
    return out

def mock(fsk, seed, Ng=23, Nc=12):
    """ΛCDM sem física nova: f_gas,500 = 0.088 (M/1e14)^0.26 (Gonzalez+2013), 0.05 dex de dispersão;
    perfis de f_gas como nos dados (razão 2500/500 e 200/500 medianas); massas e z como nas amostras reais."""
    rng = np.random.default_rng(seed); rows = []
    r25 = np.median([g[5] / g[6] for g in GR]); r200 = np.median([x[4] / x[3] for x in XC])
    rr = np.median([g[4] / g[2] for g in GR])
    for k in range(Ng):
        g = GR[rng.integers(len(GR))]; M = g[3] * 1e13 * 10**rng.normal(0, 0.1)
        f5 = min(0.088 * (M / 1e14)**0.26 * 10**rng.normal(0, .05), .16)
        rows.append(('m', 'grupo', g[1]) + system(g[1], M * MSUN, f5, f5 * r25, rr, FSTAR[fsk](M)))
    for k in range(Nc):
        x = XC[rng.integers(len(XC))]; M = x[2] * 1e14 * 10**rng.normal(0, 0.1)
        f5 = min(0.088 * (M / 1e14)**0.26 * 10**rng.normal(0, .05), .16)
        rows.append(('m', 'X-COP', x[1]) + system(x[1], M * MSUN, f5, f5 * r200, 1.53, FSTAR[fsk](M)))
    return evaluate(rows)

R = dict(H0=H0_KMS)
print(f'H0 = {H0_KMS}')
for fsk in ('G13', 'B14'):
    o = evaluate(build(fsk)); R[f'dados_{fsk}'] = o
    print(f"\n== Dados reais, f_* = {fsk}: ξ_b grupos {o['_xi_b']['grupos'][0]:.3f}–{o['_xi_b']['grupos'][1]:.3f}, X-COP {o['_xi_b']['xcop'][0]:.3f}–{o['_xi_b']['xcop'][1]:.3f};"
          f" f_bar(R500) grupos {o['_fbar']['grupos']:.3f}, X-COP {o['_fbar']['xcop']:.3f}")
    for m in ('v5', 'nulo1', 'nulo2'):
        print(f"  {m:6s} LOO (35) {o[m]['loo']:.3f} dex | X-COP→grupos {o[m]['cruzada']:.3f} dex (viés {o[m]['vies_cruzada']:+.3f}) | parâmetros {np.round(o[m]['param_total'],3).tolist()}")
    ms = [mock(fsk, s) for s in range(200)]
    R[f'mock_{fsk}'] = {m: dict(loo_med=float(np.median([x[m]['loo'] for x in ms])), cruz_med=float(np.median([x[m]['cruzada'] for x in ms])))
                        for m in ('v5', 'nulo1', 'nulo2')}
    win = np.mean([x['v5']['cruzada'] < min(x['nulo1']['cruzada'], x['nulo2']['cruzada']) for x in ms])
    gain = np.array([min(x['nulo1']['cruzada'], x['nulo2']['cruzada']) - x['v5']['cruzada'] for x in ms])
    real_gain = min(o['nulo1']['cruzada'], o['nulo2']['cruzada']) - o['v5']['cruzada']
    R[f'mock_{fsk}']['frac_v5_vence'] = float(win); R[f'mock_{fsk}']['ganho_real'] = float(real_gain)
    R[f'mock_{fsk}']['p_ganho'] = float(np.mean(gain >= real_gain))
    print(f"  MOCK ΛCDM (200 réplicas, sem física nova): mediana X-COP→grupos v5 {R[f'mock_{fsk}']['v5']['cruz_med']:.3f}, nulo1 {R[f'mock_{fsk}']['nulo1']['cruz_med']:.3f}, nulo2 {R[f'mock_{fsk}']['nulo2']['cruz_med']:.3f};"
          f" v5 vence em {100*win:.0f}% das réplicas;  ganho real do v5 sobre o melhor nulo {real_gain:+.3f} dex (P_mock ≥ real = {R[f'mock_{fsk}']['p_ganho']:.2f})")
json.dump(R, open(f'results/groups2_resultados_H{H0_KMS:g}.json', 'w'), indent=1, ensure_ascii=False)
