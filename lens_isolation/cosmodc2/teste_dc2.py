# catalogo ARTIFICIAL: testa o codigo, nao diz nada sobre CDM
import numpy as np, medir_dc2 as M
def make(flip, seed=0):
    rng=np.random.default_rng(seed); G=[]; H=[]; tag=1
    for p in range(500):
        c=(rng.uniform(0,6),rng.uniform(-38,-32),rng.uniform(0.1,0.45))
        for q in range(1+rng.poisson(6)):
            deg=np.degrees(1.5/M.chi(c[2]))
            ra=c[0]+rng.normal(0,deg)/np.cos(np.radians(c[1])); dec=c[1]+rng.normal(0,deg); z=c[2]+rng.normal(0,0.0004)
            lm=rng.normal(12.6,0.5); G.append((ra,dec,z,10**(10.7+0.4*(lm-12.6)+rng.normal(0,0.2)),True,10**lm,tag,rng.uniform(17,20)))
            if lm>13: H.append((ra,dec,z,10**lm))
            for s in range(rng.poisson(2)):
                dg=np.degrees(0.3/M.chi(z))
                G.append((ra+rng.normal(0,dg),dec+rng.normal(0,dg),z,10**rng.normal(10.0,0.3),False,10**lm,tag,rng.uniform(18,22)))
            tag+=1
    for k in range(1500):
        z=rng.uniform(0.1,0.45); G.append((rng.uniform(0,6),rng.uniform(-38,-32),z,10**rng.normal(10.7,0.3),True,10**rng.normal(12.2,0.3),tag,rng.uniform(17,20))); tag+=1
    for k in range(60000):
        G.append((rng.uniform(0,6),rng.uniform(-38,-32),rng.uniform(0.05,0.6),10**rng.normal(9.3,0.3),False,0.0,-1,rng.uniform(21,22.5)))
    a=np.array(G,float).T
    d=dict(ra=a[0],dec=a[1],z=a[2],zr=a[2].copy(),ms=a[3],cen=a[4].astype(bool),mh=a[5],tag=a[6].astype(np.int64),mag=a[7])
    ns=700000; sra=rng.uniform(0,6,ns); sdec=rng.uniform(-38,-32,ns); sz=rng.uniform(0.4,1.0,ns); g1=np.zeros(ns); g2=np.zeros(ns)
    for ra,dec,z,m in H:
        dx=(sra-ra)*np.cos(np.radians(dec)); dy=sdec-dec
        k=np.where((np.abs(dx)<1)&(np.abs(dy)<1)&(sz>z+0.1))[0]
        cl=M.chi(z); R=np.maximum(np.radians(np.hypot(dx[k],dy[k]))*cl,0.1); cs=M.chi(sz[k])
        sc=1.6624e18*cs/(cl*(cs-cl)*(1+z)); gt=(m/(2*np.pi*R*1.0))/sc; phi=np.arctan2(dy[k],dx[k])
        g1[k]+=-gt*np.cos(2*phi); g2[k]+=-flip*gt*np.sin(2*phi)
    d.update(s_ra=sra,s_dec=sdec,s_z=sz,s_g1=g1,s_g2=g2)
    return d
for flip in (1,-1):
    d=make(flip); area=M.Area(d['ra'],d['dec']); tree=M.cKDTree(np.c_[area.xy(d['s_ra'],d['s_dec'])])
    s,info=M.shear_sign(d,area,tree); print('flip',flip,'->',info)
d=make(1); edges=np.linspace(0,6,5)
for j in range(4):
    k=(d['ra']>=edges[j]-(j==0)*9)&(d['ra']<edges[j+1]+(j==3)*9)
    gk={x:d[x][k] for x in ['ra','dec','z','zr','mag','ms','cen','mh','tag']}
    ks=(d['s_ra']>=edges[j])&(d['s_ra']<edges[j+1])
    gk.update({x:d[x][ks] for x in ['s_ra','s_dec','s_z','s_g1','s_g2']})
    np.savez('teste_tile_%d.npz'%j,**gk)
files=['teste_tile_%d.npz'%j for j in range(4)]
print('--- completo'); r=M.run(files,'teste.json')
print('--- sem erros de photo-z e massa (eps deve ser claramente positivo)')
d=M.load(files); area=M.Area(d['ra'],d['dec']); pos,Xi,Xe,_=M.environment(d,area)
S=M.select(d,area,pos,sz=0.0,sn=0.0); a=M.analyse(S,Xi[pos[S['lens']]],4)['todas']; b=M.analyse(S,Xe[pos[S['lens']]],4)['todas']
print(S['iso'].sum(),len(S['iso']),a['eps_0.3_1Mpc'],a['erro'],'| sem hosp',b['eps_0.3_1Mpc'],b['erro'])
tree=M.cKDTree(np.c_[area.xy(d['s_ra'],d['s_dec'])]); print(M.delta_sigma(d,area,tree,S,1,4)['todas'])
print('--- sem corte (eps deve ser 0)')
d=M.load(files); area=M.Area(d['ra'],d['dec']); pos,Xi,Xe,_=M.environment(d,area)
S=M.select(d,area,pos,nocut=True); print(M.analyse(S,Xi[pos[S['lens']]],4)['todas'])
tree=M.cKDTree(np.c_[area.xy(d['s_ra'],d['s_dec'])]); print(M.delta_sigma(d,area,tree,S,1,4)['todas'])
