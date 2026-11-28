"""Uso: py medir_satelites.py smdpl_z032.npz satelites.json   (precisa de medir_caixa.py na mesma pasta)"""
import sys, json
import numpy as np
from scipy.spatial import cKDTree
from medir_caixa import H, L, Z, OM, EDGES, CORE, shell_mass, synthetic

def measure(c,sz=0.02,sm_noise=0.12,use_obs=False,seed=1,nmax=20000,nj=50):
    rng=np.random.default_rng(seed)
    pos=np.mod(c['pos'].astype(float)/H,L); pos[pos>=L]=0
    gal=np.where((c['sm']>10**7.9)|(c['obs_sm']>10**7.9))[0]
    base=c['obs_sm'] if use_obs else c['sm']
    lms=np.log10(np.maximum(base[gal],1.0))+(0 if use_obs else rng.normal(0,sm_noise,len(gal)))
    sig=299792.458*sz*(1+Z)/(100*H*np.sqrt(OM*(1+Z)**3+1-OM))
    obs=pos[gal].copy(); obs[:,2]=np.mod(obs[:,2]+rng.normal(0,sig,len(gal)),L); obs[obs>=L]=0
    lmb=lms+np.log10(1+10**(-0.69*lms+6.63))
    sel=[]
    for lo,hi in zip(EDGES[:-1],EDGES[1:]):
        k=np.where((lmb>=lo)&(lmb<hi))[0]
        sel.append(rng.choice(k,min(nmax,len(k)),replace=False))
    lens=np.concatenate(sel); bin_of=np.concatenate([np.full(len(s),i) for i,s in enumerate(sel)])
    t=cKDTree(obs,boxsize=L); R=4.0*70/73; iso=np.ones(len(lens),bool)
    for k,i in enumerate(lens):
        nb=np.asarray(t.query_ball_point(obs[i],R),int); nb=nb[nb!=i]
        iso[k]=not np.any(lms[nb]>lms[i]-1.0)
    cen=np.where(c['central']&(c['m']>1e11))[0]
    tc=cKDTree(pos[cen],boxsize=L); mh=c['m'][cen]/H; hh=c['host'][cen]
    tp=pos[gal][lens]; own=c['host'][gal][lens]; is_cen=c['central'][gal][lens]
    XA=shell_mass(tc,mh,hh,tp,own)
    XR=shell_mass(tc,mh,hh,rng.uniform(0,L,tp.shape),np.full(len(lens),-1))
    XB=XA.copy(); s=np.where(~is_cen)[0]
    if len(s): XB[s]=shell_mass(tc,mh,hh,tp[s],np.full(len(s),-1))
    XA=XA-XR; XB=XB-XR
    m=lms[lens]; mbins=np.linspace(7.5,13,56); jk=np.minimum((tp[:,0]/L*nj).astype(int),nj-1)
    def eps(X,s,keep):
        a=s&iso&keep; b=s&keep
        if a.sum()<20: return np.nan
        ha=np.histogram(m[a],mbins)[0]; hb=np.histogram(m[b],mbins)[0]
        w=np.divide(ha,hb,out=np.zeros(len(ha)),where=hb>0)[np.clip(np.digitize(m[b],mbins)-1,0,54)]
        Eb=(X[b]*w[:,None]).sum(0)/w.sum()
        return 1-X[a].mean(0)[CORE].sum()/Eb[CORE].sum()
    def err(v):
        v=np.array(v); v=v[np.isfinite(v)]
        return float(np.sqrt((len(v)-1)/len(v)*((v-v.mean())**2).sum()))
    names=['todas']+[f'logMb_{lo}-{hi}' for lo,hi in zip(EDGES[:-1],EDGES[1:])]
    res={}
    for j,name in enumerate(names):
        s=np.ones(len(lens),bool) if j==0 else bin_of==j-1
        all_=np.ones(len(lens),bool)
        out={}
        for tag,X,sub in (('A_original',XA,s),('B_com_hospedeiro',XB,s),('C_so_centrais',XA,s&is_cen)):
            out[tag]=[round(float(eps(X,sub,all_)),3),round(err([eps(X,sub,jk!=q) for q in range(nj)]),3)]
        d=[eps(XB,s,jk!=q)-eps(XA,s,jk!=q) for q in range(nj)]
        out['B_menos_A']=[round(float(eps(XB,s,all_)-eps(XA,s,all_)),3),round(err(d),3)]
        out['frac_satelites_isoladas']=round(float((~is_cen)[s&iso].mean()),3) if (s&iso).any() else None
        out['frac_satelites_controle']=round(float((~is_cen)[s].mean()),3)
        res[name]=out
    return res

if __name__=='__main__':
    variants={'base':{},'pz018':dict(sz=0.018),'m0':dict(sm_noise=0),'m020':dict(sm_noise=0.2),'obs_sm':dict(use_obs=True)}
    if sys.argv[1]=='--teste':
        from medir_caixa import measure as orig
        c=synthetic(); n=len(c['m']); rng=np.random.default_rng(5)
        sat=rng.random(n)<0.3; c['central']=~sat
        cid=np.where(~sat)[0]; c['host']=np.arange(n)
        for i in np.where(sat)[0]:
            h=cid[min(np.searchsorted(cid,i),len(cid)-1)]; c['host'][i]=h; c['pos'][i]=c['pos'][h]+rng.normal(0,0.2,3)
        c['pos']=np.mod(c['pos'],400)
        a=orig(c,nmax=3000)['amostras']['todas']['eps_0.3_1Mpc']; r=measure(c,nmax=3000,nj=10)['todas']
        print('medir_caixa:',a,'| A:',r['A_original'][0],'(devem ser iguais)'); print(r)
    else:
        c=dict(np.load(sys.argv[1])); out={}
        for tag,kw in variants.items():
            out[tag]=measure(c,**kw); print(tag,out[tag]['todas'],flush=True)
        json.dump(out,open(sys.argv[2],'w'),indent=1)
