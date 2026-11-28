# Uso: py medir_dc2.py            (le os dc2_*.npz da pasta; grava dc2_resultado.json)
import sys, json, glob, time
import numpy as np
from scipy.spatial import cKDTree
from scipy.ndimage import distance_transform_edt

SHELLS = np.array([0.1,0.3,0.5,0.7,1.0,1.5,2.0,3.0]); CORE = (SHELLS[:-1]>=0.3)&(SHELLS[1:]<=1.0)
EDGES = [9.0,10.5,10.8,11.1,11.5]; RCUT = 4.0*70/73; CELL = 0.05
MB = np.linspace(8,13,26); ZB = np.linspace(0.1,0.5,9)
G = np.linspace(0,3.0,6001); INV = 1/np.sqrt(0.2793*(1+G)**3+1-0.2793)
CHI = np.r_[0, np.cumsum((INV[1:]+INV[:-1])/2*np.diff(G))]*299792.458/73.0
def chi(z): return np.interp(z,G,CHI)

def xyz(ra,dec,d):
    a,b = np.radians(ra), np.radians(dec)
    return np.c_[d*np.cos(b)*np.cos(a), d*np.cos(b)*np.sin(a), d*np.sin(b)]

def load(files):
    keys = ['ra','dec','z','zr','mag','ms','cen','mh','tag','s_ra','s_dec','s_z','s_g1','s_g2']
    d = {k: [] for k in keys}; reg = []
    for j,f in enumerate(files):
        a = np.load(f)
        for k in keys: d[k].append(a[k])
        reg.append(np.full(len(a['ra']), j))
    d = {k: np.concatenate(v) for k,v in d.items()}; d['reg'] = np.concatenate(reg)
    return d

class Area:
    def __init__(self, ra, dec):
        self.ra0 = float(np.median(ra))
        x,y = self.xy(ra,dec)
        self.x0, self.y0 = x.min()-2*CELL, y.min()-2*CELL
        i = ((x-self.x0)/CELL).astype(int); j = ((y-self.y0)/CELL).astype(int)
        occ = np.zeros((i.max()+3, j.max()+3), bool); occ[i,j] = True
        self.occ = occ; self.dist = distance_transform_edt(occ)*CELL
        self.deg2 = occ.sum()*CELL**2
    def xy(self, ra, dec): return (ra-self.ra0)*np.cos(np.radians(dec)), dec
    def edge(self, ra, dec):
        x,y = self.xy(ra,dec)
        i = np.clip(((x-self.x0)/CELL).astype(int), 0, self.occ.shape[0]-1)
        j = np.clip(((y-self.y0)/CELL).astype(int), 0, self.occ.shape[1]-1)
        return self.dist[i,j] - CELL
    def randoms(self, marg, rng):
        ci,cj = np.nonzero(self.occ); dd = self.dist[ci,cj]-CELL
        o = np.argsort(-dd); ci,cj,dd = ci[o],cj[o],dd[o]
        n = np.searchsorted(-dd, -marg, side='left')
        k = (rng.random(len(marg))*np.maximum(n,1)).astype(int)
        x = self.x0+(ci[k]+rng.random(len(k)))*CELL; y = self.y0+(cj[k]+rng.random(len(k)))*CELL
        return self.ra0 + x/np.cos(np.radians(y)), y

def shell_mass(tree, mass, tag, pts, own):
    out = np.zeros((len(pts), len(SHELLS)-1)); found = np.zeros(len(pts), bool)
    for a in range(0, len(pts), 5000):
        for i,idx in enumerate(tree.query_ball_point(pts[a:a+5000], SHELLS[-1]), a):
            if not idx: continue
            idx = np.asarray(idx); mine = tag[idx]==own[i]; found[i] = mine.any(); idx = idx[~mine]
            r = np.linalg.norm(tree.data[idx]-pts[i], axis=1)
            out[i] = np.histogram(r, SHELLS, weights=mass[idx])[0]
    return out, found

def environment(d, area, nrand=5, seed=7):
    rng = np.random.default_rng(seed)
    c = np.where(d['cen'] & (d['mh']>1e11))[0]
    tree = cKDTree(xyz(d['ra'][c], d['dec'][c], chi(d['z'][c])))
    mh, tag = d['mh'][c], d['tag'][c]
    cand = np.where((d['mag']<20)&(d['z']>0.05)&(d['z']<0.6))[0]
    marg = np.degrees(4.0/chi(np.maximum(d['z'][cand],0.08)))
    cand = cand[area.edge(d['ra'][cand], d['dec'][cand]) > 0.6*marg]
    marg = np.degrees(4.0/chi(np.maximum(d['z'][cand],0.08)))
    pts = xyz(d['ra'][cand], d['dec'][cand], chi(d['z'][cand]))
    cent = d['cen'][cand]
    exc, found = shell_mass(tree, mh, tag, pts, d['tag'][cand])      # sem o proprio halo
    inc = exc.copy(); s = np.where(~cent)[0]
    if len(s): inc[s] = shell_mass(tree, mh, tag, pts[s], np.full(len(s), -999999999))[0]
    R = np.zeros_like(exc)
    for k in range(nrand):
        rra,rdec = area.randoms(marg, rng)
        R += shell_mass(tree, mh, tag, xyz(rra, rdec, chi(d['z'][cand])), np.full(len(cand), -999999999))[0]
    R /= nrand
    pos = np.full(len(d['z']), -1); pos[cand] = np.arange(len(cand))
    info = {'frac_sat_com_hospedeiro_a_3Mpc': float(found[s].mean()) if len(s) else None}
    return pos, inc-R, exc-R, info

def select(d, area, pos, sz=0.02, sn=0.12, rviz=22.5, seed=1, nocut=False):
    rng = np.random.default_rng(seed); n = len(d['z'])
    zp = d['zr'] + rng.normal(0, 1, n)*sz*(1+d['zr'])
    lms = np.log10(np.maximum(d['ms'],1.0)) + rng.normal(0, 1, n)*sn
    marg = np.degrees(4.0/chi(np.maximum(np.minimum(zp, d['z']), 0.08)))
    lens = np.where((d['mag']<20)&(zp>0.1)&(zp<0.5)&(pos>=0)&(area.edge(d['ra'],d['dec'])>marg))[0]
    iso = np.ones(len(lens), bool)
    if not nocut:
        nb_all = np.where(d['mag']<rviz)[0]
        P = xyz(d['ra'][nb_all], d['dec'][nb_all], chi(np.clip(zp[nb_all],0,None)))
        tree = cKDTree(P); m = lms[nb_all]; back = np.full(n,-1); back[nb_all] = np.arange(len(nb_all))
        Q = xyz(d['ra'][lens], d['dec'][lens], chi(zp[lens]))
        for k,i in enumerate(lens):
            idx = np.asarray(tree.query_ball_point(Q[k], RCUT), int)
            idx = idx[idx!=back[i]]
            iso[k] = not np.any(m[idx] > lms[i]-1.0)
    mst = lms[lens]; lmb = mst + np.log10(1+10**(-0.69*mst+6.63))
    return {'lens':lens, 'iso':iso, 'mst':mst, 'zp':zp[lens], 'lmb':lmb, 'reg':d['reg'][lens]}

def weights(S, sel, keep):
    a = sel&S['iso']&keep; b = sel&keep
    if a.sum()<20: return None
    Ha = np.histogram2d(S['mst'][a], S['zp'][a], [MB,ZB])[0]; Hb = np.histogram2d(S['mst'][b], S['zp'][b], [MB,ZB])[0]
    W = np.divide(Ha, Hb, out=np.zeros_like(Ha), where=Hb>0)
    i = np.clip(np.digitize(S['mst'][b], MB)-1, 0, len(MB)-2); j = np.clip(np.digitize(S['zp'][b], ZB)-1, 0, len(ZB)-2)
    return a, b, W[i,j]

def jk(vals):
    v = np.array(vals, float); v = np.where(np.isfinite(v), v, np.nan)
    n = np.isfinite(v).sum(0); mean = np.nansum(v, 0)/np.maximum(n,1)
    return np.where(n>1, np.sqrt((n-1)/np.maximum(n,1)*np.nansum((v-mean)**2, 0)), np.nan)

def bins(S):
    out = {'todas': np.ones(len(S['lens']), bool)}
    for lo,hi in zip(EDGES[:-1],EDGES[1:]): out['logMb_%.1f-%.1f'%(lo,hi)] = (S['lmb']>=lo)&(S['lmb']<hi)
    return out

def analyse(S, X, nreg):
    np.seterr(invalid='ignore', divide='ignore')
    res = {}; every = np.ones(len(S['lens']), bool); J = {}
    for name,sel in bins(S).items():
        def calc(keep):
            w = weights(S, sel, keep)
            if w is None: return np.full(len(SHELLS)+1, np.nan)
            a,b,wb = w
            Ei = X[a].mean(0); Ec = (X[b]*wb[:,None]).sum(0)/wb.sum()
            return np.r_[1-Ei/Ec, 1-Ei[CORE].sum()/Ec[CORE].sum(), Ec[CORE].sum()]
        full = calc(every)
        if not np.isfinite(full).any(): res[name] = {'n_isoladas': int((sel&S['iso']).sum())}; continue
        Js = np.array([calc(S['reg']!=q) for q in range(nreg)]); J[name] = Js[:,-2]
        err = jk(Js); snr = float(full[-1]/err[-1]) if err[-1]>0 else 0.0
        res[name] = {'eps_0.3_1Mpc': round(float(full[-2]),3), 'erro': round(float(err[-2]),3),
            'sn_excesso_controle': round(snr,1), 'inconclusivo': bool(snr<3),
            'eps_por_casca': np.round(full[:-2],3).tolist(), 'erro_por_casca': np.round(err[:-2],3).tolist(),
            'n_lentes': int(sel.sum()), 'n_isoladas': int((sel&S['iso']).sum())}
    names = list(bins(S))[1:]
    if all(n in J for n in names):
        f = [res[n]['eps_0.3_1Mpc'] for n in names]
        d = (J[names[0]]+J[names[1]])/2 - (J[names[2]]+J[names[3]])/2
        res['baixa_menos_alta'] = [round((f[0]+f[1]-f[2]-f[3])/2,3), round(float(jk(d)),3)]
    return res

def pairs(d, area, tree, ra, dec, zl, rmin, rmax):
    cl = chi(zl); tmax = np.degrees(rmax/cl)
    x,y = area.xy(ra, dec)
    idx = np.asarray(tree.query_ball_point([x,y], 1.15*tmax), int)
    if not len(idx): return None
    idx = idx[d['s_z'][idx] > zl+0.2]
    dx = (d['s_ra'][idx]-ra)*np.cos(np.radians(dec)); dy = d['s_dec'][idx]-dec
    R = np.radians(np.hypot(dx,dy))*cl
    ok = (R>=rmin)&(R<rmax); idx,dx,dy = idx[ok],dx[ok],dy[ok]
    cs = chi(d['s_z'][idx])
    sc = 1.6624e18*cs/(cl*(cs-cl)*(1+zl))          # Msun/Mpc2, comovel
    return idx, np.arctan2(dy,dx), sc

def shear_sign(d, area, tree, nmax=3000):
    h = np.where(d['cen']&(d['mh']>1e13)&(d['z']>0.1)&(d['z']<0.5))[0]
    marg = np.degrees(3.0/chi(d['z'][h])); h = h[area.edge(d['ra'][h],d['dec'][h])>marg]
    h = h[np.random.default_rng(3).permutation(len(h))[:nmax]]
    acc = np.zeros(4); n = 0
    for i in h:
        p = pairs(d, area, tree, d['ra'][i], d['dec'][i], d['z'][i], 0.5, 2.0)
        if p is None: continue
        idx,phi,sc = p; c2,s2 = np.cos(2*phi), np.sin(2*phi)
        for k,s in enumerate((1,-1)):
            acc[2*k] += (-(d['s_g1'][idx]*c2 + s*d['s_g2'][idx]*s2)).sum()
            acc[2*k+1] += (d['s_g1'][idx]*s2 - s*d['s_g2'][idx]*c2).sum()
        n += len(idx)
    acc /= max(n,1)
    sign = 1 if acc[0]>acc[2] else -1
    return sign, {'n_halos': int(len(h)), 'g_t(+)': float(acc[0]), 'g_x(+)': float(acc[1]),
                  'g_t(-)': float(acc[2]), 'g_x(-)': float(acc[3]), 'sinal': sign}

def delta_sigma(d, area, tree, S, sign, nreg):
    n = len(S['lens']); num = np.zeros(n); den = np.zeros(n)
    for k,i in enumerate(S['lens']):
        p = pairs(d, area, tree, d['ra'][i], d['dec'][i], d['z'][i], 1.0, 3.0)
        if p is not None:
            idx,phi,sc = p
            gt = -(d['s_g1'][idx]*np.cos(2*phi) + sign*d['s_g2'][idx]*np.sin(2*phi))
            num[k] = (gt/sc).sum(); den[k] = (1/sc**2).sum()
        if k%20000==0: print('  Delta-Sigma %d/%d' % (k,n), flush=True)
    res = {}; every = np.ones(n, bool)
    for name,sel in bins(S).items():
        def calc(keep):
            w = weights(S, sel, keep)
            if w is None: return np.full(3, np.nan)
            a,b,wb = w
            Di = num[a].sum()/den[a].sum(); Dc = (num[b]*wb).sum()/(den[b]*wb).sum()
            return np.array([1-Di/Dc, Di, Dc])
        full = calc(every)
        if not np.isfinite(full).all(): res[name] = None; continue
        err = jk(np.array([calc(S['reg']!=q) for q in range(nreg)]))
        res[name] = {'eps_1_3Mpc': round(float(full[0]),3), 'erro': round(float(err[0]),3),
                     'DS_isoladas_Msun_pc2': round(float(full[1])/1e12,3), 'DS_controle_Msun_pc2': round(float(full[2])/1e12,3),
                     'sn_controle': round(float(full[2]/err[2]),1) if err[2]>0 else 0.0}
    return res

VAR = {'base':{}, 'pz018':dict(sz=0.018), 'm0':dict(sn=0.0), 'm020':dict(sn=0.2), 'viz_r20':dict(rviz=20.0)}

def run(files, out):
    t0 = time.time(); d = load(files); nreg = len(files); area = Area(d['ra'], d['dec'])
    diag = {'recortes': nreg, 'area_deg2': round(float(area.deg2),1),
            'r20_por_deg2': round(float((d['mag']<20).sum()/area.deg2)), 'r225_por_deg2': round(float((d['mag']<22.5).sum()/area.deg2))}
    print(diag, flush=True)
    tree = cKDTree(np.c_[area.xy(d['s_ra'], d['s_dec'])])
    sign, diag['sinal_shear'] = shear_sign(d, area, tree); print('sinal do shear:', diag['sinal_shear'], flush=True)
    pos, Xinc, Xexc, info = environment(d, area); diag.update(info); print(info, '%.0f s' % (time.time()-t0), flush=True)
    res = {'diagnostico': diag, 'variantes': {}}
    for tag,kw in VAR.items():
        S = select(d, area, pos, **kw); k = pos[S['lens']]
        frac = float((~d['cen'][S['lens']]).mean())
        V = {'n_lentes': len(S['lens']), 'n_isoladas': int(S['iso'].sum()), 'frac_satelites_controle': round(frac,3),
             'frac_satelites_isoladas': round(float((~d['cen'][S['lens']][S['iso']]).mean()),3) if S['iso'].any() else None,
             'principal': analyse(S, Xinc[k], nreg), 'secundaria_sem_hospedeiro': analyse(S, Xexc[k], nreg)}
        res['variantes'][tag] = V; json.dump(res, open(out,'w'), indent=1)
        a = V['principal']['todas']; b = V['secundaria_sem_hospedeiro']['todas']
        print(tag, 'isoladas %d/%d | PRINCIPAL %s +- %s (S/N %s) | sem hosp. %s +- %s' % (V['n_isoladas'], V['n_lentes'],
              a.get('eps_0.3_1Mpc'), a.get('erro'), a.get('sn_excesso_controle'), b.get('eps_0.3_1Mpc'), b.get('erro')), flush=True)
        if tag=='base': Sbase = S
    print('faixas (base):', [(n, v.get('eps_0.3_1Mpc'), v.get('erro'), v.get('n_isoladas')) for n,v in res['variantes']['base']['principal'].items() if isinstance(v,dict)])
    print('baixa - alta:', res['variantes']['base']['principal'].get('baixa_menos_alta'))
    res['delta_sigma_secundario'] = delta_sigma(d, area, tree, Sbase, sign, nreg)
    print('Delta-Sigma:', res['delta_sigma_secundario']['todas'])
    json.dump(res, open(out,'w'), indent=1); print('fim: %.0f min' % ((time.time()-t0)/60))
    return res

if __name__=='__main__':
    files = sorted(glob.glob('dc2_*.npz'))
    if len(files)<30 and '--parcial' not in sys.argv:
        sys.exit('so %d de 30 recortes; espere a extracao terminar' % len(files))
    run(files, 'dc2_resultado.json')
