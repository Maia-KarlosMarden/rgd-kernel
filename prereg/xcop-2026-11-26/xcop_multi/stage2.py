import json, numpy as np, pandas as pd, warnings
from common import *
from src import ingest
warnings.filterwarnings('ignore')
MASSES=['M_FORW','M_NFW','M_EIN','M_ISO','M_BUR']
cfg=json.load(open(OUT/'E1_config.json')); lnlam=cfg['lnlam']
P=pd.read_csv(OUT/'E1_barions.csv'); pred=pd.read_csv(OUT/'E1_previsoes_logK.csv'); gal=pd.read_csv(OUT/'E1_galaxias.csv')
mk=pd.read_csv(OUT/'E1_mocks_lcdm.csv'); LAWS=[c for c in pred.columns if c not in('system_id','r_kpc','x')]
tab1={x[0]:x for x in ingest.XCOP}
chk=[]; mass={m:np.empty(len(P)) for m in MASSES}; relerr=np.empty(len(P))
for c in ALL:
    with fits.open(one(c,'hydro_mass')) as h:
        t=h[1].data; r=np.array(t['RADIUS'],float); hdr=h[1].header; idx=np.where(P.system_id.values==c)[0]
        for m in MASSES: mass[m][idx]=loglog(RAD,r,np.array(t[m],float))
        relerr[idx]=loglog(RAD,r,np.array(t['EM_FORW'],float))/mass['M_FORW'][idx]
        key=[k for k in tab1 if k.upper().replace(' ','').replace('-','')[:5] in c.upper() or c.upper()[:5] in k.upper().replace(' ','')]
        row=dict(cluster=c,R500_header=float(hdr['R500']),M500_header=float(hdr['M500']))
        if len(key)==1:
            n,z,r5,m5,r2,m2,*_=tab1[key[0]]; row.update(tab=key[0],R500_tab_transcrito=r5,M500_tab_transcrito=m5,unid=(r2,m2))
            for rr,mm in [(r5*1000,m5),(r2*1000,m2)]:
                row.setdefault('dif_NFW_pct',[]).append(round(float(100*(loglog(rr,r,np.array(t['M_NFW'],float))/(mm*1e14)-1)),2))
        chk.append(row)
chk=pd.DataFrame(chk); chk.to_csv(OUT/'E2_checagem_tabela.csv',index=False); print(chk.to_string(index=False))
for m in MASSES: P[m]=mass[m]
P['relerr_forw']=relerr
for m in MASSES: P['logK_'+m]=np.log10(P[m]/P.Mb)
P.to_csv(OUT/'E2_pontos.csv',index=False)
print('K min:',{m:round(float((10**P['logK_'+m]).min()),2) for m in MASSES},' erro rel. mediano M_FORW:',round(float(np.median(relerr)),3))
samples={'primaria':PRIM,'primaria+A2319':PRIM+[STRESS],'12':ALL}
rng=np.random.default_rng(1); NB=10000
def med_ci(v,g=None):
    v=np.asarray(v); b=np.median(v[rng.integers(0,len(v),(NB,len(v)))],1)
    if g is not None: g=np.asarray(g); b=b-np.median(g[rng.integers(0,len(g),(NB,len(g)))],1); obs=np.median(v)-np.median(g)
    else: obs=np.median(v)
    return float(obs),float(np.quantile(b,.025)),float(np.quantile(b,.975))
res=[]; per=[]; laws=[]
for m in MASSES:
    y=P['logK_'+m].values; s=slopes(P,y,lnlam)
    per.append(s[['I1','R1','Dlog10K','n']].assign(mass=m).reset_index())
    for sn,S in samples.items():
        for q in ['I1','R1']:
            v=s.loc[S,q].dropna(); o,lo,hi=med_ci(v); d,dlo,dhi=med_ci(v,gal[q])
            ps={}
            for fg in (0.10,0.12,0.14):
                x=mk[(mk.f==fg)&(mk['sample']==sn)][q].values; n=len(x)
                ps[f'p_f{fg:.2f}']=min(1,2*min((1+(x<=o).sum())/(n+1),(1+(x>=o).sum())/(n+1)))
            res.append(dict(mass=m,sample=sn,q=q,n=len(v),med_cl=o,lo=lo,hi=hi,D=d,D_lo=dlo,D_hi=dhi,**ps))
        # leis transferidas
        mask=P.system_id.isin(S).values; e2={L:pd.Series((pred[L].values-y)[mask]**2).groupby(P.system_id.values[mask]).mean() for L in LAWS}
        ix=rng.integers(0,len(S),(NB,len(S)))
        oi=res[-2]
        rows=[]
        for L in LAWS:
            rec=dict(mass=m,sample=sn,law=L,rms=np.sqrt(e2[L].mean()),bias=float(np.mean((pred[L].values-y)[mask])))
            sl=slopes(P,pred[L].values,lnlam); rec['I1_prev']=sl.loc[S,'I1'].median(); rec['I1_no_IC']=bool(oi['lo']<=rec['I1_prev']<=oi['hi'])
            for ctl in ['RGD','RAR']:
                a=e2[L].values; b=e2[ctl].values; d=np.sqrt(a[ix].mean(1))-np.sqrt(b[ix].mean(1))
                gain=np.sqrt(b.mean())-np.sqrt(a.mean()); p=min(1,2*min((1+(d>=0).sum())/(NB+1),(1+(d<=0).sum())/(NB+1)))
                rec.update({f'ganho_{ctl}':gain,f'gfrac_{ctl}':gain/np.sqrt(b.mean()),f'lo_{ctl}':np.quantile(d,.025),f'hi_{ctl}':np.quantile(d,.975),f'p_{ctl}':p})
            rows.append(rec)
        t=pd.DataFrame(rows)
        cand=t.index[~t.law.isin(['RGD','RAR'])]
        allp=[(i,ctl,t.loc[i,f'p_{ctl}']) for i in cand for ctl in ['RGD','RAR']]; order=sorted(allp,key=lambda a:a[2]); run=0
        for k,(i,ctl,p) in enumerate(order): run=max(run,min(1,(len(order)-k)*p)); t.loc[i,f'pholm_{ctl}']=run
        t['P2']=False
        for i in cand:
            t.loc[i,'P2']=all((t.loc[i,f'ganho_{c}']>=0.005)&(t.loc[i,f'gfrac_{c}']>=0.05)&(t.loc[i,f'hi_{c}']<0)&(t.loc[i,f'pholm_{c}']<0.05) for c in ['RGD','RAR'])&t.loc[i,'I1_no_IC']
        laws.append(t)
R=pd.DataFrame(res); R.to_csv(OUT/'E2_contraste.csv',index=False)
LW=pd.concat(laws); LW.to_csv(OUT/'E2_leis.csv',index=False)
pd.concat(per).to_csv(OUT/'E2_por_cluster.csv',index=False)
lk=[]
for m in MASSES:
    for sn,S in samples.items():
        o=np.median(P.loc[P.system_id.isin(S),'logK_'+m])
        for fg in (0.10,0.12,0.14):
            x=mk[(mk.f==fg)&(mk['sample']==sn)].logK.values; lk.append(dict(mass=m,sample=sn,f=fg,logK_obs=o,logK_mock=np.median(x),dif=o-np.median(x)))
pd.DataFrame(lk).to_csv(OUT/'E2_nivel_mock.csv',index=False)
print(R.round(3).to_string(index=False))
print(LW[LW['mass']=='M_FORW'][['sample','law','rms','bias','ganho_RGD','ganho_RAR','pholm_RGD','pholm_RAR','I1_prev','I1_no_IC','P2']].round(3).to_string(index=False))
print('P2 verdadeiros:\n',LW[LW.P2==True][['mass','sample','law','rms']].to_string(index=False))
print(pd.DataFrame(lk).query("sample=='primaria'").round(3).to_string(index=False))
print(per[0].round(2).to_string(index=False))
y=P.logK_M_FORW; pm=P[P.system_id.isin(PRIM)]
print(pm.groupby('r_kpc')[['x','logK_M_FORW','logK_M_NFW']].median().round(3).to_string())
print(pred[pred.system_id.isin(PRIM)].groupby('r_kpc')[['RGD','RAR','F3']].median().round(3).to_string())
