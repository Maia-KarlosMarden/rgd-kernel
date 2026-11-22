import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))  # paths relativos a v5/code
"""RGD v5 — aprimoramentos propostos ao kernel RGD v4 (Maia 2026) e seus testes.

Kernel (inalterado):  K_eff = 1 + λ q/(1-q)^n,  q = exp(-g_b/g†(z)),  g† = cH(z)/2π,  λ = 0.999

(1) CONTROLE DE CIRCULARIDADE: mocks sem física nova (só f_bar(M) conhecido) → qual r(n_emp, ξ) e qual A saem?
(2) VARIÁVEL DO HOSPEDEIRO SÓ COM BÁRIONS:
      R_b  : raio onde a densidade bariônica média = Δ_b ρ_c(z),  Δ_b = 500·f_cos  (f_cos = Ω_b/Ω_m = 0.157)
      ξ_b  = G M_b(<R_b) / (R_b² g†(z))  ∝ M_b^{1/3}   (não usa massa total)
    lei:  n(ξ_b) = 1/2 + A_b ξ_b, calibrada nos 12 X-COP em R_b; LOO; comparada com v4 (ξ com M_HSE)
    e com o modelo nulo "amostra justa" (M_tot = M_b / f_b,universal).
(3) PREVISÃO SEM REAJUSTE: (a) SPARC com n_i = 1/2 + A_b ξ_b,i;  (b) BCGs CLASH com n do hospedeiro.
Dados X-COP: Eckert+2019 (A&A 621, A40) Tab. 2: M_HSE,500, f_gas,500, f_gas,200; Ettori+2019 Tab.1: M200/M500.
z dos X-COP: valores de catálogo (≈, efeito em E(z) < 4%).
"""
import sys, io, contextlib, math, json, runpy
import numpy as np
from scipy.optimize import minimize_scalar, brentq
H0_KMS = float(sys.argv[1]) if len(sys.argv) > 1 else 70.0
sys.argv = ['x']; sys.path.insert(0, '.')
with contextlib.redirect_stdout(io.StringIO()):
    import audit2 as A
    CT = runpy.run_path('clusters_test.py')
c = 2.998e8; G = 6.674e-11; MSUN = 1.989e30; MPC = 3.0857e22; KPC = MPC / 1e3
H0 = H0_KMS * 1e3 / MPC; E = lambda z: math.sqrt(.3 * (1 + z)**3 + .7)
gd = lambda z: c * H0 * E(z) / (2 * math.pi); rhoc = lambda z: 3 * (H0 * E(z))**2 / (8 * math.pi * G)
LAM, FCOS = 0.999, 0.157
def K(g, z, n): x = g / gd(z); return 1 + LAM * np.exp(-x) / (-np.expm1(-x))**n
def n_invert(Kobs, x):  # n tal que K(x;n)=Kobs
    return math.log(LAM * math.exp(-x) / (Kobs - 1)) / math.log(-math.expm1(-x))
L = np.log10; R = {}

# ---------------- (1) mocks ----------------
rng = np.random.default_rng(3)
def mock(trend=True, N=5000, scatter=0.05):
    M = 10**rng.uniform(13.3, 15.3, N); z = rng.uniform(0.03, 0.7, N)
    fb = (1.09 * 0.0708 * (M / 1e13)**0.135 if trend else np.full(N, 0.13)) * 10**rng.normal(0, scatter, N)
    fb = np.minimum(fb, 0.17)
    g500 = G * (4 * math.pi / 3 * 500 * np.array([rhoc(zz) for zz in z]))**(2 / 3) * (M * MSUN)**(1 / 3)
    xi = g500 / np.array([gd(zz) for zz in z])
    n = np.array([n_invert(1 / f, f * x) for f, x in zip(fb, xi)])
    A_ = np.sum((n - .5) * xi) / np.sum(xi**2)
    return float(np.corrcoef(n, xi)[0, 1]), float(A_), float(np.sqrt(np.mean((n - .5 - A_ * xi)**2)))
print('== (1) Controle de circularidade: mocks com bárions comuns, NENHUMA física nova ==')
for lab, tr in [('f_bar(M) de Sun+2009 (escala usada no artigo)', True), ('f_bar constante = 0.13', False)]:
    r_, A_, s_ = mock(tr); R[f'mock_{tr}'] = dict(r=r_, A=A_, rms_n=s_)
    print(f'  {lab:48s} r(n_emp, ξ) = {r_:.3f}   A ajustado = {A_:.4f}   dispersão em n = {s_:.3f}')
print('  artigo v4:                                       r = 0.97          A = 0.4068')

# ---------------- (2) X-COP com variável bariônica ----------------
X = [  # nome, z, M_HSE500 [1e14], f_gas500, f_gas200
 ('A1644', .0473, 3.48, .128, .126), ('A1795', .0622, 4.63, .139, .144), ('A2029', .0766, 8.65, .141, .152),
 ('A2142', .0909, 8.95, .158, .168), ('A2255', .0809, 5.26, .153, .146), ('A2319', .0557, 7.31, .189, .237),
 ('A3158', .0597, 4.26, .145, .155), ('A3266', .0589, 8.80, .132, .108), ('A644', .0704, 5.66, .132, .139),
 ('A85', .0555, 5.65, .150, .159), ('RXC1825', .0650, 4.08, .133, .155), ('ZwCl1215', .0766, 7.66, .106, .092)]
M200_M500, R200_R500 = 1.45, 1.53          # NFW c≈3.3 (Ettori+2019): razões típicas
s_tot = math.log(M200_M500) / math.log(R200_R500)
rows = []
for nm, z, M5, f5, f2 in X:
    M5 *= 1e14 * MSUN; R5 = (3 * M5 / (4 * math.pi * 500 * rhoc(z)))**(1 / 3)
    Mb5 = 1.09 * f5 * M5                                        # gás + estrelas (Ettori+2019: f_bar = 1.09 f_gas)
    s_b = math.log(f2 * M200_M500 / f5) / math.log(R200_R500)   # inclinação local de M_b(<r) entre R500 e R200
    s_b = min(max(s_b, 0.6), 1.6)
    # raio bariônico: Mb(<R)/(4π/3 R³) = Δ_b ρc,  com Mb ∝ r^s_b
    Rb = R5 * ((1.09 * f5) / FCOS)**(1 / (3 - s_b)); Mb = Mb5 * (Rb / R5)**s_b
    Mt = M5 * (Rb / R5)**s_tot                                  # massa hidrostática no mesmo raio (alvo)
    gb, gobs = G * Mb / Rb**2, G * Mt / Rb**2
    xi_b = gb / gd(z)                                           # = g_b(R_b)/g†  — só bárions
    xi_t = (G * M5 / R5**2) / gd(z)                             # ξ do v4 (massa total)
    rows.append(dict(nome=nm, z=z, Rb_R500=Rb / R5, xi_b=xi_b, xi_v4=xi_t, gb=gb, gobs=gobs, Mb=Mb / MSUN,
                     n_emp=n_invert(gobs / gb, gb / gd(z))))
xi_b = np.array([r['xi_b'] for r in rows]); xi_t = np.array([r['xi_v4'] for r in rows])
gb = np.array([r['gb'] for r in rows]); gobs = np.array([r['gobs'] for r in rows]); zz = np.array([r['z'] for r in rows])
Mb = np.array([r['Mb'] for r in rows]); nemp = np.array([r['n_emp'] for r in rows])
print('\n== (2) X-COP: variável do hospedeiro só com bárions ==')
print(f"  R_b/R500 = {np.min([r['Rb_R500'] for r in rows]):.2f}–{np.max([r['Rb_R500'] for r in rows]):.2f};  ξ_b = {xi_b.min():.3f}–{xi_b.max():.3f};  ξ_v4 = {xi_t.min():.2f}–{xi_t.max():.2f};  n_emp(R_b) = {nemp.min():.3f}–{nemp.max():.3f}")
print(f'  corr(n_emp, ξ_b) = {np.corrcoef(nemp, xi_b)[0,1]:+.3f}   corr(n_emp, ξ_v4) = {np.corrcoef(nemp, xi_t)[0,1]:+.3f}   corr(n_emp, f_bar) = {np.corrcoef(nemp, [1.09*x[3] for x in X])[0,1]:+.3f}')
def fitA(xv, msk):
    f = lambda a: np.sum((L(gb[msk] * np.array([K(g, z_, .5 + a * x) for g, z_, x in zip(gb[msk], zz[msk], xv[msk])])) - L(gobs[msk]))**2)
    return minimize_scalar(f, bounds=(-20, 60), method='bounded').x
def loo(model):
    e = []
    for i in range(len(rows)):
        m = np.arange(len(rows)) != i
        if model == 'nulo':
            fb_u = np.median(gb[m] / gobs[m]); p = L(gb[i] / fb_u)
        else:
            xv = xi_b if model == 'v5' else xi_t; a = fitA(xv, m); p = L(gb[i] * K(gb[i], zz[i], .5 + a * xv[i]))
        e.append(p - L(gobs[i]))
    return float(np.sqrt(np.mean(np.square(e))))
allm = np.ones(len(rows), bool)
A_b = fitA(xi_b, allm); A_v4 = fitA(xi_t, allm)
R['xcop'] = dict(A_b=A_b, A_v4_refit=A_v4, linhas=rows,
                 loo=dict(v5_xi_bar=loo('v5'), v4_xi_total=loo('v4'), nulo_amostra_justa=loo('nulo'), galactico_n05=float(np.sqrt(np.mean((L(gb * K(gb, 0, .5)) - L(gobs))**2)))))
print(f'  A_b (v5, ξ bariônico) = {A_b:.3f};   A (v4 reajustado neste raio) = {A_v4:.3f}')
print('  LOO RMS em g_obs(R_b):  ' + '  '.join(f'{k} {v:.3f} dex' for k, v in R['xcop']['loo'].items()))

# ---------------- (3a) SPARC: previsão sem reajuste ----------------
meta = A.meta; names = A.names
Mb_g = np.array([(0.5 * meta[n]['L'] + 1.33 * meta[n]['MHI']) * 1e9 for n in names]) * MSUN
xi_g = G * (4 * math.pi / 3 * 500 * FCOS * rhoc(0))**(2 / 3) * Mb_g**(1 / 3) / gd(0)
n_g = 0.5 + A_b * xi_g
def rmsg(p): e = p - A.y; return math.sqrt(np.mean([np.mean(e[A.idx == i]**2) for i in range(A.ng)]))
base = rmsg(L(A.gb * K(A.gb, 0, 0.5))); v5 = rmsg(L(A.gb * K(A.gb, 0, n_g[A.idx])))
# BTFR: inclinação implicada para as galáxias mais massivas
R['sparc'] = dict(xi_b_min=float(xi_g.min()), xi_b_max=float(xi_g.max()), n_min=float(n_g.min()), n_max=float(n_g.max()), rms_n05=base, rms_v5=v5)
print(f'\n== (3a) SPARC, previsão sem reajuste ==\n  ξ_b das galáxias = {xi_g.min():.4f}–{xi_g.max():.4f}  →  n = {n_g.min():.3f}–{n_g.max():.3f}')
print(f'  RMS/gal: n=1/2 {base:.4f} dex | v5 n(ξ_b) {v5:.4f} dex  (Δ = {v5-base:+.4f})')
print(f'  galáxia mais massiva: v² ∝ R^(2n−1) = R^{2*n_g.max()-1:+.3f} no regime profundo (BTFR exige 0)')

# ---------------- (3b) BCGs CLASH ----------------
# M_b do hospedeiro: f_bar(R500) da mediana X-COP × M500 de kT (Arnaud+2005) — estimativa; entra como M_b^{1/3}
fb_med = float(np.median([1.09 * x[3] for x in X]))
out = []
for k, t in enumerate(CT['T']):
    z = t[1]; M5 = 3.84e14 * (t[7] / 5)**1.71 / E(z) * MSUN
    xi_h = G * (4 * math.pi / 3 * 500 * FCOS * rhoc(z))**(2 / 3) * (fb_med * M5)**(1 / 3) / gd(z)
    n_h = 0.5 + A_b * xi_h; g_ = CT['gbar'][k]
    out.append((L(CT['gtot'][k] / g_), L(K(g_, z, n_h)), n_h, xi_h))
o = np.array(out)
R['bcg'] = dict(xi_host=[float(o[:, 3].min()), float(o[:, 3].max())], n_host=[float(o[:, 2].min()), float(o[:, 2].max())],
                obs_med=float(np.median(o[:, 0])), pred_med=float(np.median(o[:, 1])), rms=float(np.sqrt(np.mean((o[:, 1] - o[:, 0])**2))))
# n necessário em cada BCG
n_need = [n_invert(10**ob, g_ / gd(t[1])) for (ob, *_), g_, t in zip(out, CT['gbar'], CT['T'])]
R['bcg']['n_necessario_mediano'] = float(np.median(n_need))
print(f"\n== (3b) BCGs CLASH com n do hospedeiro (sem lei radial) ==\n  ξ_b(hospedeiro) {o[:,3].min():.3f}–{o[:,3].max():.3f} → n {o[:,2].min():.3f}–{o[:,2].max():.3f}")
print(f"  log(g_obs/g_b): observado {R['bcg']['obs_med']:.3f}  previsto {R['bcg']['pred_med']:.3f}  RMS {R['bcg']['rms']:.3f} dex;  n necessário mediano {R['bcg']['n_necessario_mediano']:.2f}")
json.dump(R, open(f'results/rgd_v5_resultados_H{H0_KMS:g}.json', 'w'), indent=1, ensure_ascii=False, default=float)
