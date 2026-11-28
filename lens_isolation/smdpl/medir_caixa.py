"""Uso: medir_caixa.py smdpl_z032.npz saida.json  |  medir_caixa.py --teste x"""
import sys, json
import numpy as np
from scipy.spatial import cKDTree
H=0.6777; L=400/H; Z=0.315; OM=0.307
SHELLS=np.array([0.1,0.3,0.5,0.7,1.0,1.5,2.0,3.0]); EDGES=[9.0,10.5,10.8,11.1,11.5]
CORE=(SHELLS[:-1]>=0.3)&(SHELLS[1:]<=1.0)

def shell_mass(tree,mass,host,pts,own):
    out=np.zeros((len(pts),len(SHELLS)-1))
    for i,idx in enumerate(tree.query_ball_point(pts,SHELLS[-1])):
        idx=np.asarray(idx,int)
        idx=idx[host[idx]!=own[i]]
        d=tree.data[idx]-pts[i]; d-=L*np.round(d/L)
        out[i]=np.histogram(np.linalg.norm(d,axis=1),SHELLS,weights=mass[idx])[0]
    return out

def measure(c,sz=0.02,sm_noise=0.12,use_obs=False,seed=1,nmax=20000,sem_corte=False,nj=50):
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
    if not sem_corte:
        for k,i in enumerate(lens):
            nb=np.asarray(t.query_ball_point(obs[i],R),int); nb=nb[nb!=i]
            iso[k]=not np.any(lms[nb]>lms[i]-1.0)
    cen=np.where(c['central']&(c['m']>1e11))[0]
    tc=cKDTree(pos[cen],boxsize=L); mh=c['m'][cen]/H; hh=c['host'][cen]
    tp=pos[gal][lens]
    X=shell_mass(tc,mh,hh,tp,c['host'][gal][lens])
    X=X-shell_mass(tc,mh,hh,rng.uniform(0,L,tp.shape),np.full(len(lens),-1))
    m=lms[lens]; mbins=np.linspace(7.5,13,56); jk=np.minimum((tp[:,0]/L*nj).astype(int),nj-1)
    def eps(s,keep):
        a=s&iso&keep; b=s&keep
        if a.sum()<20: return None
        ha=np.histogram(m[a],mbins)[0]; hb=np.histogram(m[b],mbins)[0]
        w=np.divide(ha,hb,out=np.zeros(len(ha)),where=hb>0)[np.clip(np.digitize(m[b],mbins)-1,0,54)]
        Ea=X[a].mean(0); Eb=(X[b]*w[:,None]).sum(0)/w.sum()
        return np.r_[1-Ea/Eb,1-Ea[CORE].sum()/Eb[CORE].sum()]
    res={'n_lentes':len(lens),'n_isoladas':int(iso.sum()),'cascas_Mpc':SHELLS.tolist(),'amostras':{}}
    names=['todas']+[f'logMb_{lo}-{hi}' for lo,hi in zip(EDGES[:-1],EDGES[1:])]
    for j,name in enumerate(names):
        s=np.ones(len(lens),bool) if j==0 else bin_of==j-1
        full=eps(s,np.ones(len(lens),bool))
        if full is None: res['amostras'][name]=None; continue
        J=np.array([e for e in (eps(s,jk!=q) for q in range(nj)) if e is not None])
        err=np.sqrt((len(J)-1)/len(J)*((J-J.mean(0))**2).sum(0))
        tot=X[s][:,CORE].sum(1); snr=float(tot.mean()/(tot.std()/np.sqrt(len(tot))))
        res['amostras'][name]={'eps_0.3_1Mpc':round(float(full[-1]),3),'erro':round(float(err[-1]),3),
            'eps_por_casca':np.round(full[:-1],3).tolist(),'erro_por_casca':np.round(err[:-1],3).tolist(),
            'sn_excesso_controle':round(snr,1),'inconclusivo':bool(snr<3),'n_isoladas':int((s&iso).sum())}
    return res

def synthetic(seed=0):
    rng=np.random.default_rng(seed); P=[];M=[];S=[]
    for p in range(6000):
        c0=rng.uniform(0,400,3); n=1+rng.poisson(10)
        P.append(c0+rng.normal(0,1.0,(n,3))); lm=rng.normal(11.8,0.5,n); M.append(10**lm)
        S.append(10**(10.0+0.8*(lm-11.8)+rng.normal(0,0.3,n)))
    P=np.mod(np.vstack(P),400); M=np.concatenate(M); S=np.concatenate(S); n=len(M)
    return {'pos':P.astype('f4'),'m':M,'sm':S,'obs_sm':S,'central':np.ones(n,bool),'host':np.arange(n)}

if __name__=='__main__':
    if sys.argv[1]=='--teste':
        c=synthetic()
        for tag,kw in {'sem corte (deve dar 0)':dict(sem_corte=True),'posições exatas':dict(sz=0,sm_noise=0),'com erros':{}}.items():
            r=measure(c,nmax=3000,**kw); a=r['amostras']['todas']
            print(tag,r['n_isoladas'],'/',r['n_lentes'],a and (a['eps_0.3_1Mpc'],a['erro']))
    else:
        c=dict(np.load(sys.argv[1])); out={}
        for tag,kw in {'base':{},'pz018':dict(sz=0.018),'m0':dict(sm_noise=0),'m020':dict(sm_noise=0.2),'obs_sm':dict(use_obs=True)}.items():
            out[tag]=measure(c,**kw); print(tag,out[tag]['amostras']['todas'])
        json.dump(out,open(sys.argv[2],'w'),indent=1)
