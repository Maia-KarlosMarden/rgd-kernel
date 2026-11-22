import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))  # paths relativos a v5/code
"""Fase 5 — teste radiativo em 20 clusters CLASH (escala da BCG, mesma abertura).
Dados transcritos de fontes publicadas:
  Tian, Umetsu, Ko, Donahue & Chiu (2020), ApJ 896, 70, Tabela 1: r (kpc), M*_BCG, M_gas, M_tot em 1e11 Msun
     (M_tot de lenteamento forte+fraco, Umetsu+2016; M_gas Chandra, Donahue+2014; M* Cooke+2016)
  Postman et al. (2012), ApJS 199, 25, Tabela 4: kT (keV), L_bol 0.1–100 keV (1e44 erg/s)
Pergunta: a discrepância de massa dos clusters depende da luminosidade, como exige g_ind ∝ f(radiação)?
"""
import math, json
import numpy as np
from scipy import stats
from scipy.optimize import minimize_scalar

G = 6.67430e-11; MSUN = 1.98847e30; KPC = 3.085677581491367e19; c = 299792458.0
# nome, z, r_kpc, Mstar, Mgas(None=sem medida), Mtot, eMtot   [1e11 Msun]   | kT, eT, Lbol, eL  (Postman)
T = [
 ('Abell 383',        .187, 14.3,  4.45, 1.78,  7.55, 2.23,  6.5, .24,  6.7, .2),
 ('Abell 209',        .206, 14.3,  4.85, None,  3.87, 0.79,  7.3, .54, 12.7, .3),
 ('Abell 2261',       .224, 23.6, 12.30, 0.48,  6.44, 1.48,  7.6, .30, 18.0, .2),
 ('RX J2129.7+0005',  .234, 14.3,  5.81, 2.18,  6.65, 1.91,  5.8, .40, 11.4, 2.0),
 ('Abell 611',        .288, 22.2,  6.58, 0.48,  6.24, 1.81,  7.9, .35, 11.7, .2),
 ('MS2137-2353',      .313, 14.3,  3.65, 2.94,  3.98, 1.57,  5.9, .30,  9.9, .3),
 ('RX J2248.7-4431',  .348, 30.3,  8.09, 1.01,  6.19, 2.16, 12.4, .60, 69.5, .1),
 ('MACS J1115.9+0129',.355, 16.2,  3.00, 5.80,  6.25, 1.53,  8.0, .40, 21.1, .4),
 ('MACS J1931.8-2635',.352, 14.3,  6.92, 1.47,  7.21, 2.90,  6.7, .40, 20.9, .6),
 ('RX J1532.9+3021',  .362, 14.3,  3.34, 1.13,  6.80, 4.18,  5.5, .40, 20.5, .9),
 ('MACS J1720.3+3536',.387, 23.6,  6.59, 1.15,  6.83, 2.07,  6.6, .40, 13.3, .5),
 ('MACS J0416.1-2403',.397, 14.3,  3.14, None,  4.22, 0.94,  7.5, .80, 16.0, .9),
 ('MACS J0429.6-0253',.399, 17.2, 11.90, 6.71,  9.98, 3.40,  6.0, .44, 11.2, .5),
 ('MACS J1206.2-0847',.439, 14.3,  3.13, None,  6.90, 2.07, 10.8, .60, 43.0, 1.0),
 ('MACS J0329.7-0211',.450, 22.2,  8.47, 13.1, 25.30, 6.50,  8.0, .50, 17.0, .6),
 ('RX J1347.5-1145',  .451, 14.3,  4.52, 5.11,  7.28, 1.88, 15.5, .60, 90.8, 1.0),
 ('MACS J1149.5+2223',.544, 14.3,  4.72, None,  4.54, 1.12,  8.7, .90, 30.2, 1.2),
 ('MACS J0717.5+3745',.548, 14.3,  2.19, None,  4.07, 0.73, 12.5, .70, 55.8, 1.1),
 ('MACS J0647.7+7015',.584, 14.3, 14.70, None,  7.71, 2.77, 13.3, 1.8, 32.5, 2.1),
 ('MACS J0744.9+3927',.686, 14.3,  7.74, None,  7.65, 2.45,  8.9, .80, 29.1, 1.2),
]
names = [t[0] for t in T]
r = np.array([t[2] for t in T]) * KPC
Ms = np.array([t[3] for t in T]); Mg = 0.1 * np.array([t[4] or 0.0 for t in T]);  # coluna M_gas em 1e10 Msun: validado pelo máximo publicado g_bar=2.1e-10 e por 'gás subdominante em r<R_e'
has_gas = np.array([t[4] is not None for t in T])
Mt = np.array([t[5] for t in T]); eMt = np.array([t[6] for t in T])
kT = np.array([t[7] for t in T]); L = np.array([t[9] for t in T]) * 1e37  # W
z = np.array([t[1] for t in T])
Re = np.array([17.0,22.1,18.7,41.4,30.4,15.0,34.5,52.9,33.2,21.8,17.2,56.2,29.3,44.8,22.9,21.7,34.3,13.2,56.9,14.7]) * KPC  # R_e (Tian+2020 Tab.1)
rh = 0.551 * Re
Ms = Ms * r**2 / (r + rh)**2          # Hernquist: massa estelar da BCG contida em r
gbar = G * (Ms + Mg) * 1e11 * MSUN / r**2
gtot = G * Mt * 1e11 * MSUN / r**2
sig = eMt / Mt / math.log(10)                          # incerteza log10 g_tot
sig_tot = np.sqrt(sig**2 + (0.1 / math.log(10))**2)    # + 10% em M* (Tian+2020)

# --- RAR de cluster, inclinação 1/2: g_tot = sqrt(g‡ g_bar)
ya, xa = np.log10(gtot), np.log10(gbar)
w = 1 / sig_tot**2
lg = np.sum(w * (2 * ya - xa)) / np.sum(w)            # log g‡ (média ponderada)
res = ya - 0.5 * (xa + lg)
print(f'g‡ (inclinação 1/2, escala BCG, esta reconstrução) = {10**lg:.3e} m/s²   [Tian+2020 com todos os raios: 2.02e-9]')
print(f'espalhamento dos resíduos: {np.std(res):.3f} dex;  g_tot/g_bar mediano = {np.median(gtot/gbar):.2f}')
a_eff = gtot**2 / gbar                                 # g‡ efetivo por cluster (limite profundo)

# --- grandezas radiativas na mesma abertura
g_gamma = G * L / (c**3 * r)                           # monopolo de E/c² em trânsito (fonte central, limite superior de ordem)
dg = gtot - gbar
K_req = dg / g_gamma
out = dict(g_dagger_cluster=float(10**lg), scatter_dex=float(np.std(res)), clusters={})
print(f'\n{"cluster":20s} {"L_bol(1e44)":>11s} {"g_bar":>9s} {"g_tot":>9s} {"a_eff":>9s} {"g_gamma":>9s} {"K exigido":>10s}')
for i, nm in enumerate(names):
    print(f'{nm:20s} {L[i]/1e37:11.1f} {gbar[i]:9.2e} {gtot[i]:9.2e} {a_eff[i]:9.2e} {g_gamma[i]:9.2e} {K_req[i]:10.2e}')
    out['clusters'][nm] = dict(L_bol_1e44=float(L[i] / 1e37), kT=float(kT[i]), r_kpc=float(r[i] / KPC), g_bar=float(gbar[i]), g_tot=float(gtot[i]),
                               a_eff=float(a_eff[i]), g_gamma=float(g_gamma[i]), K_req=float(K_req[i]), gas_medido=bool(has_gas[i]))

# --- testes de dependência
def corr(xv, yv, lab):
    rs, ps = stats.spearmanr(xv, yv); rp, pp = stats.pearsonr(xv, yv)
    sl = np.polyfit(xv, yv, 1)[0]
    print(f'  {lab:46s} Spearman {rs:+.2f} (p={ps:.3f})  Pearson {rp:+.2f} (p={pp:.3f})  inclinação {sl:+.3f}')
    return dict(spearman=float(rs), p_spearman=float(ps), pearson=float(rp), p_pearson=float(pp), slope=float(sl))
lL = np.log10(L); pos = K_req > 0; lK = np.where(pos, np.log10(np.abs(K_req)), np.nan)
print(f'\nvalidação: max g_bar = {gbar.max():.2e} (Tian+2020: ≈2.1e-10);  clusters com Δg>0: {pos.sum()}/20')
print('\n== Dependências (20 clusters) ==')
out['corr'] = {
  'res_RAR_vs_logL': corr(lL, res, 'resíduo da RAR de cluster vs log L_bol'),
  'log_aeff_vs_logL': corr(lL, np.log10(a_eff), 'log a_eff (a0′ por cluster) vs log L_bol'),
  'log_aeff_vs_logT': corr(np.log10(kT), np.log10(a_eff), 'log a_eff vs log kT'),
  'logK_vs_logL': corr(lL[pos], lK[pos], 'log K exigido vs log L_bol  (hipótese: 0)'),
  'res_vs_L_per_Mtot': corr(np.log10(L / Mt), res, 'resíduo vs log(L_bol/M_tot)'),
}
m = has_gas
print('  -- só clusters com gás medido na região da BCG (n=%d)' % m.sum())
out['corr']['res_vs_logL_gas'] = corr(lL[m], res[m], 'resíduo vs log L_bol (gás medido)')
mp = m & pos; out['corr']['logK_vs_logL_gas'] = corr(lL[mp], lK[mp], 'log K vs log L_bol (gás medido)')

# --- comparação de modelos (validação deixa-um-fora, χ² em log g_tot)
def loo(pred_fn, nparam):
    e = []
    for i in range(len(T)):
        tr = np.arange(len(T)) != i
        e.append(pred_fn(tr, i))
    e = np.array(e); return float(np.sum((e / sig_tot)**2)), float(np.sqrt(np.mean(e**2)))
def rar_pred(tr, i):
    lg_ = np.sum(w[tr] * (2 * ya[tr] - xa[tr])) / np.sum(w[tr]); return ya[i] - 0.5 * (xa[i] + lg_)
def rad_pred(tr, i):   # g_tot = g_bar + K g_gamma, K universal
    f = lambda lk: np.sum(w[tr] * (ya[tr] - np.log10(gbar[tr] + 10**lk * g_gamma[tr]))**2)
    lk = minimize_scalar(f, bounds=(5, 15), method='bounded').x
    return ya[i] - np.log10(gbar[i] + 10**lk * g_gamma[i])
def rad_s_pred(tr, i):  # g_tot = g_bar + K g_gamma^s g_bar^(1-s) ... forma livre: K * (L/Lref)^s * g_bar^(1/2)
    Lr = np.median(L[tr])
    def f(t): return np.sum(w[tr] * (ya[tr] - (0.5 * (xa[tr] + t[0]) + t[1] * np.log10(L[tr] / Lr)))**2)
    from scipy.optimize import minimize
    t = minimize(f, [lg, 0.0]).x
    return ya[i] - (0.5 * (xa[i] + t[0]) + t[1] * np.log10(L[i] / Lr))
def newton_pred(tr, i): return ya[i] - xa[i]
print('\n== Modelos, validação deixa-um-fora (20 clusters) ==')
out['loo'] = {}
for lab, fn in [('Newton (só bárions)', newton_pred), ('RAR de cluster, g‡ universal', rar_pred),
                ('Radiativo: g_bar + K·g_gamma, K universal', rad_pred), ('RAR de cluster × (L/L_ref)^s', rad_s_pred)]:
    chi, rm = loo(fn, 0); out['loo'][lab] = dict(chi2=chi, rmse_dex=rm)
    print(f'  {lab:44s} χ² {chi:8.1f}   RMSE {rm:.3f} dex')
json.dump(out, open('results/clusters_resultados.json', 'w'), indent=1, ensure_ascii=False)
