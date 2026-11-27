"""Calibra em SPARC+X-COP, congela previsoes do CLASH (so barions) e mocks LCDM. Nao le massas do CLASH."""
import sys, json, warnings, hashlib, numpy as np, pandas as pd
from pathlib import Path
from scipy.optimize import brentq
ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / 'xcop_multi', ROOT / 'correcoes'): sys.path.insert(0, str(p))
import calibrar as K, common as C
from src import fit as F, mocks as M, validation as V
from src.constants import G, MSUN, KPC, rhoc, gdag
warnings.filterwarnings('ignore')
OUT = ROOT / 'runs' / 'clash_pre'
RAD = np.array([100, 150, 200, 300, 400, 500, 600.])
STRESS = ['MACS_J0416.1-2403', 'MACS_J0647.8+7015', 'MACS_J0717+3745', 'MACS_J1149.6+2223']
ll = lambda r, rr, y: np.exp(np.interp(np.log(r), np.log(rr), np.log(y)))
# ---- calibracao
sp = M.Base().sp.copy(); xc = K.xcop_table()
S_xc = {}
for c in C.ALL:
    b = C.baryons(c); S_xc[c] = np.log10(ll(500, b['r'], b['mgas']))
for t in (sp, xc):
    t['logPhi'] = np.log10(t.gbar_ms2 * t.radius_m / 299792458.0**2)
    t['physical_system_id'] = t.catalog + ':' + t.system_id.astype(str)
sp['S'] = sp.logMb_out; xc['S'] = xc.system_id.map(S_xc)
allc = pd.concat([sp, xc], ignore_index=True)
def fitH(tr): return F.fit_h5(tr, 'S')
H = fitH(allc)
ctl = K.fit(allc, 'CTRL-local'); rar = F.fit(allc, 'RAR'); rar_s = F.fit(sp, 'RAR'); rgd_s = F.fit(sp, 'RGD')
MOD = {'H-host': lambda d: F.predict(H, d) if False else F.h5_predict(H, d),
       'CTRL-local': lambda d: K.predict(ctl, d), 'RAR_SPARC+XCOP': lambda d: F.predict(rar, d),
       'RAR_SPARC': lambda d: F.predict(rar_s, d), 'RGD_SPARC': lambda d: F.predict(rgd_s, d)}
# diagnosticos internos (dados ja vistos)
eg = np.empty(len(sp))
for tr, te in V.group_folds(sp):
    eg[te] = F.h5_predict(fitH(pd.concat([sp.iloc[tr], xc], ignore_index=True)), sp.iloc[te]) - sp.y.values[te]
ex = np.empty(len(xc))
for c in C.ALL:
    te = (xc.system_id == c).values
    ex[te] = F.h5_predict(fitH(pd.concat([sp, xc[~te]], ignore_index=True)), xc[te]) - xc.y.values[te]
tS = lambda S: 1 / (1 + np.exp(-(np.asarray(S) - H['S0']) / H['s']))
diag = dict(rms_SPARC_cv=K.rms_sys(sp, eg), rms_XCOP_loco=K.rms_sys(xc, ex), vies_XCOP=float(ex.mean()),
            min_slope=K.min_slope(dict(H, model='H5:S'), allc), S0=float(H['S0']), theta=H['theta'].tolist(),
            S_SPARC_max=float(sp.S.max()), S_XCOP=[float(min(S_xc.values())), float(max(S_xc.values()))],
            t_XCOP_min=float(tS(min(S_xc.values()))), frac_SPARC_t_gt_01=float((tS(sp.groupby('system_id').S.first()) > 0.1).mean()))
# ---- estrelas
ratios = []
for c in C.ALL:
    b = C.baryons(c)
    if b['measured']: ratios.append(ll(RAD, b['r'], b['mstar']) / ll(RAD, b['r'], b['mgas']))
ratios = np.array(ratios); assert len(ratios) == 7
STAR = {'baixo': ratios.min(0), 'mediano': np.median(ratios, 0), 'alto': ratios.max(0)}
# ---- CLASH, so barions
g = pd.read_csv(ROOT / 'runs' / 'clash' / 'mgas_accept2_H70.csv')
cov = pd.read_csv(ROOT / 'runs' / 'clash' / 'cobertura.csv'); cov = cov[cov.grupo == 'CLASH']
assert len(cov) == 20
rows = []
for _, cr in cov.iterrows():
    q = g[g.system_id == cr.cluster]; r, mg = q.r_kpc.values, q.Mgas_msun.values
    S = np.log10(ll(500, r, mg)); assert r.min() <= 100 and r.max() >= 500
    for k, R in enumerate(RAD):
        if R > r.max(): continue
        m = ll(R, r, mg)
        for sc, fr in STAR.items():
            gb = G * m * (1 + fr[k]) * MSUN / (R * KPC)**2
            rows.append(dict(catalog='CLASH', system_id=cr.cluster, z=cr.z, amostra='estresse' if cr.cluster in STRESS else 'primaria',
                             estrelas=sc, r_kpc=R, radius_m=R * KPC, gbar_ms2=gb, x=gb / gdag(), Mgas=m, S=S,
                             logPhi=np.log10(gb * R * KPC / 299792458.0**2)))
P = pd.DataFrame(rows)
for n, f in MOD.items(): P['logK_' + n] = f(P)
diag['t_CLASH_min'] = float(tS(P.S.min())); diag['S_CLASH'] = [float(P.S.min()), float(P.S.max())]
# ---- mocks LCDM
def halo(r, mg, fg, z):
    n = max(3, len(r) // 3); sl = float(np.clip(np.polyfit(np.log(r[-n:]), np.log(mg[-n:]), 1)[0], 0, 1.5))
    gas = lambda R: ll(R, r, mg) if R <= r[-1] else mg[-1] * (R / r[-1])**sl
    f = lambda lr: np.log(gas(np.exp(lr)) / (fg * 4 / 3 * np.pi * (np.exp(lr) * KPC)**3 * 500 * rhoc(z) / MSUN))
    R = np.exp(brentq(f, np.log(150), np.log(6000)))
    return 4 / 3 * np.pi * (R * KPC)**3 * 500 * rhoc(z) / MSUN
for fg in (0.10, 0.12, 0.14):
    rng = np.random.default_rng(int(fg * 1000))
    sims = np.full((2000, len(P)), np.nan)
    for c, idx in P.groupby('system_id', sort=True).indices.items():
        q = g[g.system_id == c]; z = P.z.values[idx[0]]
        M5 = halo(q.r_kpc.values, q.Mgas_msun.values, fg, z)
        M2 = 1.45 * M5 * 10**rng.normal(0, 0.15, 2000); cc = M.conc(M2, rng); noise = rng.normal(0, 0.10 / np.log(10), (2000, len(idx)))
        for j, i in enumerate(idx):
            md = (1 - 0.157) * M.nfw_enclosed(P.radius_m.values[i], M2, cc, z)
            gb = P.gbar_ms2.values[i]
            sims[:, i] = np.log10((gb + G * md * MSUN / P.radius_m.values[i]**2) / gb) + noise[:, j]
    tag = f'mock_f{int(fg*100)}'
    P[tag + '_med'] = np.median(sims, 0); P[tag + '_lo'] = np.percentile(sims, 2.5, 0); P[tag + '_hi'] = np.percentile(sims, 97.5, 0)
P.to_csv(OUT / 'previsoes_CLASH.csv', index=False)
json.dump(dict(diag=diag, estrelas={k: v.round(4).tolist() for k, v in STAR.items()}, raios=RAD.tolist(),
               CTRL=np.asarray(ctl['theta']).tolist(), RAR=rar['theta'].tolist(), RAR_SPARC=rar_s['theta'].tolist(),
               RGD_SPARC=rgd_s['theta'].tolist(), H=dict(theta=H['theta'].tolist(), S0=float(H['S0']), s=H['s'])),
          open(OUT / 'modelos_congelados.json', 'w'), indent=1)
print(json.dumps(diag, indent=1))
m = P[(P.estrelas == 'mediano') & (P.amostra == 'primaria')]
print(m.groupby('r_kpc')[['x'] + ['logK_' + n for n in MOD] + ['mock_f10_med', 'mock_f12_med', 'mock_f14_med']].median().round(3).T.to_string())
print('pontos', len(m), 'clusters', m.system_id.nunique())
for sc in STAR: print(sc, STAR[sc].round(3), 'H mediana', round(P[P.estrelas == sc]['logK_H-host'].median(), 3))
