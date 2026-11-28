"""Uso: medir_epsilon.py catalogo.(parquet|csv|fits) saida.json [--teste]"""
import sys, json, numpy as np
from scipy.spatial import cKDTree
C=299792.458; H_MOCK=0.7
SHELLS=np.array([0.1,0.3,0.5,0.7,1.0,1.5,2.0,3.0]); MB_EDGES=[9.0,10.5,10.8,11.1,11.5]

def chi(z,H0=73.0,Om=0.2793):
    g=np.linspace(0,2.0,4001); inv=1/np.sqrt(Om*(1+g)**3+1-Om)
    t=np.concatenate([[0],np.cumsum((inv[1:]+inv[:-1])/2*np.diff(g))])*C/H0
    return np.interp(z,g,t)

def xyz(ra,dec,d):
    a,b=np.radians(ra),np.radians(dec)
    return np.c_[d*np.cos(b)*np.cos(a),d*np.cos(b)*np.sin(a),d*np.sin(b)]

def shell_mass(tree,mass,own,pts,pts_own):
    out=np.zeros((len(pts),len(SHELLS)-1))
    for i,idx in enumerate(tree.query_ball_point(pts,SHELLS[-1])):
        idx=np.asarray(idx,int)
        if idx.size==0: continue
        idx=idx[own[idx]!=pts_own[i]]
        r=np.linalg.norm(tree.data[idx]-pts[i],axis=1)
        out[i]=np.histogram(r,SHELLS,weights=mass[idx])[0]
    return out

def measure(c,sz=0.02,sm=0.12,seed=1,nj=50,sem_corte=False):
    rng=np.random.default_rng(seed); n=len(c['ra'])
    lms=c['lmstellar']-2*np.log10(H_MOCK); Mh=10**(c['lmhalo']-np.log10(H_MOCK))
    z=c['z']; zp=z+rng.normal(0,sz*(1+z),n); lms_o=lms+rng.normal(0,sm,n)
    true=xyz(c['ra'],c['dec'],chi(z)); phot=xyz(c['ra'],c['dec'],chi(np.clip(zp,0,None)))
    ra0,ra1,d0,d1=c['ra'].min(),c['ra'].max(),c['dec'].min(),c['dec'].max()
    marg=np.degrees(4.0/np.maximum(chi(np.clip(zp,0.05,None)),1))
    mra=marg/np.cos(np.radians(c['dec']))
    lens=np.where((c['rmag']<20)&(zp>0.1)&(zp<0.5)&(c['ra']>ra0+mra)&(c['ra']<ra1-mra)&
                  (c['dec']>d0+marg)&(c['dec']<d1-marg))[0]
    R=4.0*70/73; tp=cKDTree(phot); iso=np.ones(len(lens),bool)
    for k,i in enumerate(lens):
        nb=np.asarray(tp.query_ball_point(phot[i],R),int); nb=nb[nb!=i]
        iso[k]=not np.any(lms_o[nb]>lms_o[i]-1.0)
    if sem_corte: iso[:]=True
    cen=np.where(c['flag_central']==0)[0]
    tc=cKDTree(true[cen]); hid=c['halo_id'][cen]
    E=shell_mass(tc,Mh[cen],hid,true[lens],c['halo_id'][lens])
    # aleatórios: mesmo z, posição angular uniforme na área interna
    rra=rng.uniform(c['ra'][lens].min(),c['ra'][lens].max(),len(lens))
    s=np.sin(np.radians([c['dec'][lens].min(),c['dec'][lens].max()]))
    rdec=np.degrees(np.arcsin(rng.uniform(s[0],s[1],len(lens))))
    Er=shell_mass(tc,Mh[cen],hid,xyz(rra,rdec,chi(z[lens])),np.full(len(lens),-1))
    X=E-Er
    mst=lms_o[lens]; f=10**(-0.69*mst+6.63); lmb=mst+np.log10(1+f)
    mb=np.linspace(8,13,26); zb=np.linspace(0.1,0.5,9)
    jk=np.minimum((nj*(c['ra'][lens]-ra0)/(ra1-ra0+1e-9)).astype(int),nj-1)
    core=(SHELLS[:-1]>=0.3)&(SHELLS[1:]<=1.0)
    def eps(sel,keep):
        a=sel&iso&keep; b=sel&keep
        if a.sum()<20 or b.sum()<20: return None
        Ha=np.histogram2d(mst[a],zp[lens][a],[mb,zb])[0]; Hb=np.histogram2d(mst[b],zp[lens][b],[mb,zb])[0]
        W=np.divide(Ha,Hb,out=np.zeros_like(Ha),where=Hb>0)
        i=np.clip(np.digitize(mst[b],mb)-1,0,len(mb)-2); j=np.clip(np.digitize(zp[lens][b],zb)-1,0,len(zb)-2)
        w=W[i,j]; Ei=X[a].mean(0); Ec=(X[b]*w[:,None]).sum(0)/w.sum()
        return np.r_[1-Ei/Ec, 1-Ei[core].sum()/Ec[core].sum()]
    res={'n_lentes':int(len(lens)),'n_isoladas':int(iso.sum()),'cascas_Mpc':SHELLS.tolist(),'amostras':{}}
    sels={'todas':np.ones(len(lens),bool)}
    for lo,hi in zip(MB_EDGES[:-1],MB_EDGES[1:]): sels[f'logMb_{lo}-{hi}']=(lmb>=lo)&(lmb<hi)
    for name,sel in sels.items():
        full=eps(sel,np.ones(len(lens),bool))
        if full is None: res['amostras'][name]=None; continue
        J=np.array([e for e in (eps(sel,jk!=q) for q in range(nj)) if e is not None])
        err=np.sqrt((len(J)-1)/len(J)*((J-J.mean(0))**2).sum(0))
        b=sel
        tot=X[b][:,core].sum(1); snr=float(tot.mean()/(tot.std()/np.sqrt(len(tot))))
        res['amostras'][name]={'sn_excesso_controle':round(snr,1),'inconclusivo':bool(snr<3),'eps_por_casca':full[:-1].round(3).tolist(),'erro_por_casca':err[:-1].round(3).tolist(),
                               'eps_0.3_1Mpc':round(float(full[-1]),3),'erro':round(float(err[-1]),3),
                               'n_isoladas':int((sel&iso).sum())}
    return res

def synthetic(seed=0):
    """Catálogo artificial com halos reunidos em superestruturas. Só testa o código; não diz nada sobre CDM."""
    rng=np.random.default_rng(seed); rows=[]; h=0
    for p in range(1200):
        pra,pdec,pz=rng.uniform(0,12),rng.uniform(0,10),rng.uniform(0.08,0.42)
        d=chi(pz); deg=np.degrees(1.0/d); dz=1.0*73/299792.458*np.sqrt(0.28*(1+pz)**3+0.72)
        for m in range(1+rng.poisson(12)):
            lm=rng.normal(12.2,0.5)
            rows.append((pra+rng.normal(0,deg),pdec+rng.normal(0,deg),pz+rng.normal(0,dz),h,0,lm,
                         10.4+rng.normal(0,0.35)+0.4*(lm-12.2),rng.uniform(17,21)))
            h+=1
    a=np.array(rows).T
    return dict(zip(['ra','dec','z','halo_id','flag_central','lmhalo','lmstellar','rmag'],a))

def load(p):
    import pandas as pd
    if p.endswith('.parquet'): d=pd.read_parquet(p)
    elif p.endswith('.fits'):
        from astropy.table import Table; d=Table.read(p).to_pandas()
    else: d=pd.read_csv(p,comment='#')
    names={'ra':'ra_gal','dec':'dec_gal','z':'z_cgal','halo_id':'unique_halo_id','flag_central':'flag_central',
           'lmhalo':'lmhalo','lmstellar':'lmstellar','rmag':'sdss_r_true'}
    return {k:d[v].to_numpy() for k,v in names.items()}

if __name__=='__main__':
    if '--teste' in sys.argv:
        c=synthetic()
        r0=measure(c,nj=10,sem_corte=True); print('sem corte (deve dar 0):',r0['amostras']['todas']['eps_por_casca'],r0['amostras']['todas']['eps_0.3_1Mpc'])
        r=measure(c,sz=0.0,sm=0.0,nj=10); r2=measure(c,nj=10)
        print('sem erros :',r['n_isoladas'],'/',r['n_lentes'],r['amostras']['todas'])
        print('com erros :',r2['n_isoladas'],'/',r2['n_lentes'],r2['amostras']['todas'])
    else:
        c=load(sys.argv[1]); out={}
        for tag,(sz,sm) in {'base':(0.02,0.12),'pz018':(0.018,0.12),'m0':(0.02,0.0),'m020':(0.02,0.2)}.items():
            out[tag]=measure(c,sz,sm)
        json.dump(out,open(sys.argv[2],'w'),indent=1); print(out['base']['amostras']['todas'])
