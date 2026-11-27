import json, numpy as np, pandas as pd, warnings
from scipy.optimize import brentq
from common import *
warnings.filterwarnings('ignore')
B={c:baryons(c) for c in ALL}
P=pd.concat([profile(B[c]) for c in ALL],ignore_index=True); P.to_csv(OUT/'E1_barions.csv',index=False)
sp=sparc(); L=fit_laws(sp); lnlam=float(L['RGD']['theta'][0])
pred=P[['system_id','r_kpc','x']].copy()
for m,mm in L.items(): pred[m]=F.predict(mm,P)
pred.to_csv(OUT/'E1_previsoes_logK.csv',index=False)
prim=P[P.system_id.isin(PRIM)]; lo,hi=np.percentile(prim.x,[5,95])
w=sp[(sp.x>=lo)&(sp.x<=hi)]; w=w[w.groupby('system_id').x.transform('size')>=4]
gal=slopes(w,w.y.values,lnlam)[['I1','R1']].dropna(); gal.to_csv(OUT/'E1_galaxias.csv')
rows=[]
for m in L:
    s=slopes(P,pred[m].values,lnlam)[['I1','R1']]
    for q in ['I1','R1']: rows.append(dict(law=m,q=q,median_primaria=s.loc[PRIM,q].median()))
ps=pd.DataFrame(rows); ps.to_csv(OUT/'E1_previsao_inclinacoes.csv',index=False)
# ---- mocks LCDM: halo definido so pelo gas
def halo(b,fg,z):
    f=lambda lr: np.log(loglog(np.exp(lr),b['r'],b['mgas'])/(fg*4/3*np.pi*(np.exp(lr)*KPC)**3*500*rhoc(z)/MSUN))
    # extrapolacao constante do gas fora da grade
    rr=brentq(f,np.log(200),np.log(5000)); R=np.exp(rr)
    return 4/3*np.pi*(R*KPC)**3*500*rhoc(z)/MSUN
mock=[]
for fg in (0.10,0.12,0.14):
    M500={c:halo(B[c],fg,B[c]['z']) for c in ALL}; rng=np.random.default_rng(int(fg*1000))
    for it in range(2000):
        y=np.empty(len(P))
        for k,c in enumerate(ALL):
            idx=np.where(P.system_id.values==c)[0]; M2=1.45*M500[c]*10**rng.normal(0,0.15); cc=float(M.conc(M2,rng))
            md=(1-0.157)*M.nfw_enclosed(RAD*KPC,M2,cc,B[c]['z'])
            gb=P.gbar_ms2.values[idx]; gt=gb+G*md*MSUN/(RAD*KPC)**2
            y[idx]=np.log10(gt/gb)+rng.normal(0,0.10/np.log(10),len(RAD))
        s=slopes(P,y,lnlam)
        for name,S in [('primaria',PRIM),('primaria+A2319',PRIM+[STRESS]),('12',ALL)]:
            mock.append(dict(f=fg,it=it,sample=name,I1=s.loc[S,'I1'].median(),R1=s.loc[S,'R1'].median(),
                             logK=np.median(y[P.system_id.isin(S).values])))
mk=pd.DataFrame(mock); mk.to_csv(OUT/'E1_mocks_lcdm.csv',index=False)
json.dump(dict(lnlam=lnlam,x_janela=[lo,hi],n_gal=len(gal),gal_I1=gal.I1.median(),gal_R1=gal.R1.median()),open(OUT/'E1_config.json','w'),indent=1)
print(open(OUT/'E1_config.json').read())
print(ps.pivot(index='law',columns='q',values='median_primaria').round(3).to_string())
print(mk.groupby(['sample','f'])[['I1','R1','logK']].quantile([.025,.5,.975]).round(3).to_string())
print('x por raio (mediana primaria):',prim.groupby('r_kpc').x.median().round(4).tolist())
print('Mstar/Mgas mediana por raio (primaria):',(prim.Mstar/prim.Mgas).groupby(prim.r_kpc).median().round(3).tolist())
