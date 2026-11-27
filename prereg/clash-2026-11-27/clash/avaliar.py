"""Aplica T1-T4 do PROTOCOLO.md a previsao congelada. Nao altera modelos nem previsoes."""
import json, hashlib, numpy as np, pandas as pd
from pathlib import Path
from astropy.cosmology import FlatLambdaCDM
from scipy.stats import theilslopes
import astropy.units as u
ROOT = Path(__file__).resolve().parents[1]; PRE = ROOT / 'runs' / 'clash_pre'; OUT = ROOT / 'runs' / 'clash_eval'
for line in (PRE / 'CONGELADO.sha256').read_text().split('\n'):
    if line.strip():
        h, f = line.split(); assert hashlib.sha256((PRE / f).read_bytes()).hexdigest() == h, f
G, MSUN, KPC = 6.67430e-11, 1.98847e30, 3.0856775814913673e19
c27, c30 = FlatLambdaCDM(70, 0.27), FlatLambdaCDM(70, 0.30)
T = pd.read_csv(OUT / 'umetsu2016_tabela2.csv'); P = pd.read_csv(PRE / 'previsoes_CLASH.csv')
mu = lambda y: np.log1p(y) - y / (1 + y)
def r200(M200, z):
    rho = c27.critical_density(z).to_value(u.kg / u.m**3)
    return (3 * M200 * 1e14 * MSUN / (4 * np.pi * 200 * rho))**(1 / 3) / KPC
def gobs(r_kpc, M200, c, z):
    f = float(c30.angular_diameter_distance(z) / c27.angular_diameter_distance(z))
    R = r200(M200, z); M = M200 * 1e14 * mu(r_kpc / f * c / R) / mu(c) * f
    return G * M * MSUN / (r_kpc * KPC)**2
P = P.merge(
    T[['cluster', 'M200', 'c200', 'eM', 'ec']], left_on='system_id', right_on='cluster', how='inner')
T['rm2_calc'] = [r200(m, P.z[P.system_id == k].iloc[0]) / c / 1000 if k in set(P.system_id) else np.nan for k, m, c in zip(T.cluster, T.M200, T.c200)]
q = (T.rm2_calc / T.rm2)
print('consistencia r-2 calc/tabela: mediana %.3f min %.3f max %.3f n=%d' % (q.median(), q.min(), q.max(), q.notna().sum()))
print(T.loc[(q - 1).abs() > 0.1, ['nome', 'rm2', 'rm2_calc']].round(3).to_string(index=False))
P['logK_obs'] = [np.log10(gobs(r.r_kpc, r.M200, r.c200, r.z) / r.gbar_ms2) for r in P.itertuples()]
MODS = ['H-host', 'CTRL-local', 'RAR_SPARC+XCOP', 'RAR_SPARC', 'RGD_SPARC', 'mock_f10', 'mock_f12', 'mock_f14']
col = lambda m: m + '_med' if m.startswith('mock') else 'logK_' + m
rng = np.random.default_rng(20261127); rows = []
for am in ('primaria', 'estresse'):
    for st in ('baixo', 'mediano', 'alto'):
        d = P[(P.amostra == am) & (P.estrelas == st)]
        for m in MODS:
            e = d[col(m)] - d.logK_obs
            per = pd.DataFrame(dict(s=d.system_id, e=e, e2=e**2)).groupby('s')[['e', 'e2']].mean()
            b = rng.integers(0, len(per), (5000, len(per)))
            br = np.sqrt(per.e2.values[b].mean(1)); bb = per.e.values[b].mean(1)
            rows.append(dict(amostra=am, estrelas=st, modelo=m, n=len(per), rms=np.sqrt(per.e2.mean()), rms_lo=np.percentile(br, 2.5), rms_hi=np.percentile(br, 97.5),
                             vies=per.e.mean(), vies_lo=np.percentile(bb, 2.5), vies_hi=np.percentile(bb, 97.5)))
M = pd.DataFrame(rows); M.to_csv(OUT / 'metricas.csv', index=False)
def t1(st):
    k = M[(M.amostra == 'primaria') & (M.estrelas == st)].set_index('modelo'); H = k.loc['H-host']
    if H.rms > 0.20 or abs(H.vies) > 0.15: return 'FALHA'
    ok = H.rms <= 0.15 and abs(H.vies) <= 0.10 and all(k.loc[c].rms - H.rms >= 0.05 for c in ('CTRL-local', 'RAR_SPARC+XCOP'))
    return 'PASSA' if ok else 'INCONCLUSIVO'
v = {s: t1(s) for s in ('baixo', 'mediano', 'alto')}
T1 = v['mediano'] if len(set(v.values())) == 1 else 'INCONCLUSIVO'
d = P[(P.amostra == 'primaria') & (P.estrelas == 'mediano')]
k = d[d.r_kpc == 300]; res = (k.logK_obs - k['logK_H-host']).values; S = k.S.values
sl = theilslopes(res, S)[0]; boot = []
for _ in range(5000):
    i = rng.integers(0, len(S), len(S))
    if len(np.unique(S[i])) > 2: boot.append(theilslopes(res[i], S[i])[0])
ci = np.percentile(boot, [2.5, 97.5])
k12 = M[(M.amostra == 'primaria') & (M.estrelas == 'mediano')].set_index('modelo')
T4 = bool(k12.loc['mock_f12'].rms <= k12.loc['H-host'].rms)
def I1(c):
    out = []
    for _, g in d.groupby('system_id'):
        K = 10**g[c].values; out.append(np.polyfit(np.log(g.r_kpc), np.log(K - 1), 1)[0])
    return float(np.median(out))
byr = d.assign(**{'e_' + m: d[col(m)] - d.logK_obs for m in ('H-host', 'CTRL-local', 'mock_f12')}).groupby('r_kpc')[['logK_obs', 'logK_H-host', 'e_H-host', 'e_CTRL-local', 'e_mock_f12']].median()
out = dict(T1=T1, T1_por_cenario=v, T3=dict(inclinacao=float(sl), ic95=ci.tolist(), n=int(len(S)), exclui_zero=bool(ci[0] > 0 or ci[1] < 0)),
           T4_LCDM_explica_igual_ou_melhor=T4, dlnKm1_dlnr=dict(obs=I1('logK_obs'), H=I1('logK_H-host')))
json.dump(out, open(OUT / 'veredito.json', 'w'), indent=1)
P.to_csv(OUT / 'pontos_avaliados.csv', index=False)
print(M[(M.estrelas == 'mediano')].drop(columns='estrelas').round(3).to_string(index=False))
print(M[(M.modelo == 'H-host') & (M.amostra == 'primaria')][['estrelas', 'rms', 'vies']].round(3).to_string(index=False))
print(byr.round(3).to_string()); print(json.dumps(out, indent=1))
