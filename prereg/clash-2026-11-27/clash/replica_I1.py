"""Replicacao de I1: CLASH (lenteamento) x X-COP (hidrostatico) nos mesmos raios. Segue runs/replica_I1/PROTOCOLO.md."""
import sys, json, hashlib, warnings, numpy as np, pandas as pd
from pathlib import Path
from scipy.optimize import brentq
ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / 'xcop_multi'): sys.path.insert(0, str(p))
import common as C
from src import mocks as M
from src.constants import G, MSUN, KPC, rhoc, gdag
warnings.filterwarnings('ignore')
OUT = ROOT / 'runs' / 'replica_I1'
assert hashlib.sha256((OUT / 'PROTOCOLO.md').read_bytes()).hexdigest() == (OUT / 'PROTOCOLO.sha256').read_text().split()[0]
RAD = np.array([100, 150, 200, 300, 400, 500, 600.])
lnlam = json.load(open(ROOT / 'runs' / 'xcop_multi' / 'E1_config.json'))['lnlam']
ll = lambda r, rr, y: np.exp(np.interp(np.log(r), np.log(rr), np.log(y)))
# CLASH
cl = pd.read_csv(ROOT / 'runs' / 'clash_eval' / 'pontos_avaliados.csv')
cl = cl[cl.amostra == 'primaria'].copy(); cl['logK'] = cl.logK_obs
# X-COP nos mesmos raios
rows = []
for c in C.PRIM:
    b = C.baryons(c); assert b['measured']
    with C.fits.open(C.one(c, 'hydro_mass')) as h:
        t = h[1].data; rm, mm = np.array(t['RADIUS'], float), np.array(t['M_FORW'], float)
    mg = ll(RAD, b['r'], b['mgas']); mb = mg + ll(RAD, b['r'], np.maximum(b['mstar'], 1e-30))
    gb = G * mb * MSUN / (RAD * KPC)**2
    for k, R in enumerate(RAD):
        rows.append(dict(catalog='XCOP', system_id=c, r_kpc=R, radius_m=R * KPC, gbar_ms2=gb[k], x=gb[k] / gdag(), z=b['z'],
                         logK=np.log10(ll(R, rm, mm)[()] / mb[k]) if np.ndim(ll(R, rm, mm)) == 0 else np.log10(ll(R, rm, mm) / mb[k])))
xc = pd.DataFrame(rows)
rng = np.random.default_rng(1)
def med(df, y):
    s = C.slopes(df, y, lnlam)
    out = {}
    for q in ('I1', 'R1'):
        v = s[q].dropna().values; b = np.median(v[rng.integers(0, len(v), (10000, len(v)))], 1)
        out[q] = dict(med=float(np.median(v)), ic=[float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))], n=int(len(v)),
                      negativos=int((v < 0).sum()))
    return out, s
res = {}
res['XCOP'], sx = med(xc, xc.logK.values)
for st in ('baixo', 'mediano', 'alto'):
    d = cl[cl.estrelas == st]; res['CLASH_' + st], s = med(d, d.logK.values)
    if st == 'mediano': per = pd.concat([s.assign(amostra='CLASH'), sx.assign(amostra='XCOP')])
def verdict(a, b):
    lo, hi = a['I1']['ic']; xlo, xhi = b['I1']['ic']
    if xlo <= 0 <= xhi: return 'INCONCLUSIVO'
    over = lo <= xhi and xlo <= hi
    if hi < 0 and over: return 'REPLICA'
    if hi >= 0 and not over: return 'NAO REPLICA'
    return 'INCONCLUSIVO'
v = {st: verdict(res['CLASH_' + st], res['XCOP']) for st in ('baixo', 'mediano', 'alto')}
res['vereditos'] = v; res['veredito'] = v['mediano'] if len(set(v.values())) == 1 else 'INCONCLUSIVO'
# mocks
g = pd.read_csv(ROOT / 'runs' / 'clash' / 'mgas_accept2_H70.csv')
gas = {c: (g[g.system_id == c].r_kpc.values, g[g.system_id == c].Mgas_msun.values) for c in cl.system_id.unique()}
for c in C.PRIM:
    b = C.baryons(c); gas[c] = (b['r'], b['mgas'])
def halo(r, mg, fg, z):
    n = max(3, len(r) // 3); sl = float(np.clip(np.polyfit(np.log(r[-n:]), np.log(mg[-n:]), 1)[0], 0, 1.5))
    gm = lambda R: ll(R, r, mg) if R <= r[-1] else mg[-1] * (R / r[-1])**sl
    f = lambda lr: np.log(gm(np.exp(lr)) / (fg * 4 / 3 * np.pi * (np.exp(lr) * KPC)**3 * 500 * rhoc(z) / MSUN))
    R = np.exp(brentq(f, np.log(150), np.log(6000)))
    return 4 / 3 * np.pi * (R * KPC)**3 * 500 * rhoc(z) / MSUN
mocks = []
for lab, d in (('CLASH', cl[cl.estrelas == 'mediano'].reset_index(drop=True)), ('XCOP', xc)):
    for fg in (0.10, 0.12, 0.14):
        r2 = np.random.default_rng(int(fg * 1000) + (0 if lab == 'CLASH' else 7))
        Y = np.empty((2000, len(d)))
        for c, idx in d.groupby('system_id', sort=True).indices.items():
            z = d.z.values[idx[0]]; M2 = 1.45 * halo(*gas[c], fg, z) * 10**r2.normal(0, 0.15, 2000); cc = M.conc(M2, r2)
            noise = r2.normal(0, 0.10 / np.log(10), (2000, len(idx)))
            for j, i in enumerate(idx):
                md = (1 - 0.157) * M.nfw_enclosed(d.radius_m.values[i], M2, cc, z); gb = d.gbar_ms2.values[i]
                Y[:, i] = np.log10((gb + G * md * MSUN / d.radius_m.values[i]**2) / gb) + noise[:, j]
        I = np.array([C.slopes(d, Y[k], lnlam).I1.median() for k in range(2000)])
        obs = res['CLASH_mediano' if lab == 'CLASH' else 'XCOP']['I1']['med']
        p = min(1, 2 * min((1 + (I <= obs).sum()) / 2001, (1 + (I >= obs).sum()) / 2001))
        mocks.append(dict(amostra=lab, f=fg, mock_med=float(np.median(I)), lo=float(np.percentile(I, 2.5)), hi=float(np.percentile(I, 97.5)),
                          obs=obs, p=p, dentro_95=bool(np.percentile(I, 2.5) <= obs <= np.percentile(I, 97.5))))
res['mocks'] = mocks
per.to_csv(OUT / 'por_cluster.csv'); xc.to_csv(OUT / 'xcop_mesmos_raios.csv', index=False)
json.dump(res, open(OUT / 'resultado.json', 'w'), indent=1)
for k in ('XCOP', 'CLASH_baixo', 'CLASH_mediano', 'CLASH_alto'): print(k, {q: (round(a['med'], 3), np.round(a['ic'], 3).tolist(), f"{a['negativos']}/{a['n']} neg") for q, a in res[k].items()})
print(res['veredito'], v); print(pd.DataFrame(mocks).round(3).to_string(index=False))
