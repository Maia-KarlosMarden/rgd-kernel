"""POS-HOC (nao pre-registrado): estrelas, vies hidrostatico radial, ruido do mock; e recalculo independente de I1."""
import numpy as np, pandas as pd, json, warnings, glob
from astropy.io import fits
from common import *
warnings.filterwarnings('ignore')
cfg=json.load(open(OUT/'E1_config.json')); lnlam=cfg['lnlam']; P=pd.read_csv(OUT/'E2_pontos.csv'); gal=pd.read_csv(OUT/'E1_galaxias.csv')
rows=[]
def add(name,mbfac_star=1.0,mass=None):
    mb=P.Mgas+mbfac_star*P.Mstar; Mt=P.M_FORW if mass is None else mass
    gb=G*mb*MSUN/P.radius_m**2; q=P.assign(x=gb/gdag(),gbar_ms2=gb); y=np.log10(Mt/mb).values
    s=slopes(q,y,lnlam); rows.append(dict(variante=name,I1=s.loc[PRIM,'I1'].median(),R1=s.loc[PRIM,'R1'].median(),D_I1=s.loc[PRIM,'I1'].median()-gal.I1.median(),
        logK=np.median(y[P.system_id.isin(PRIM)]),excesso_vs_RAR=np.median((y-F.predict(F.fit(SP,'RAR'),q))[P.system_id.isin(PRIM)])))
SP=sparc()
add('registrado'); add('estrelas x0',0.0); add('estrelas x0.5',0.5); add('estrelas x2 (inclui ICL hipotetica)',2.0); add('estrelas x3',3.0)
for b0,b1 in [(0.05,0.15),(0.0,0.3)]:
    b=b0+(b1-b0)*np.log(P.r_kpc/100)/np.log(10); add(f'vies HSE {b0}->{b1} de 0,1 a 1 Mpc',mass=P.M_FORW/(1-b))
t=pd.DataFrame(rows); t.to_csv(OUT/'PH_sensibilidade.csv',index=False); print(t.round(3).to_string(index=False))
# mock com ruido = erro mediano real
mk=[]
noise=float(np.median(P.relerr_forw))
exec_src=open('stage1.py').read()
B={c:baryons(c) for c in ALL}
from scipy.optimize import brentq
for fg in (0.10,0.12,0.14):
    M500={}
    for c,b in B.items():
        f=lambda lr: np.log(loglog(np.exp(lr),b['r'],b['mgas'])/(fg*4/3*np.pi*(np.exp(lr)*KPC)**3*500*rhoc(b['z'])/MSUN))
        R=np.exp(brentq(f,np.log(200),np.log(5000))); M500[c]=4/3*np.pi*(R*KPC)**3*500*rhoc(b['z'])/MSUN
    rng=np.random.default_rng(int(fg*1000)); P0=pd.read_csv(OUT/'E1_barions.csv')
    for it in range(2000):
        y=np.empty(len(P0))
        for c in ALL:
            idx=np.where(P0.system_id.values==c)[0]; M2=1.45*M500[c]*10**rng.normal(0,0.15); cc=float(M.conc(M2,rng))
            gb=P0.gbar_ms2.values[idx]; gt=gb+G*(1-0.157)*M.nfw_enclosed(RAD*KPC,M2,cc,B[c]['z'])*MSUN/(RAD*KPC)**2
            y[idx]=np.log10(gt/gb)+rng.normal(0,noise/np.log(10),len(RAD))
        s=slopes(P0,y,lnlam)
        for sn,S in [('primaria',PRIM),('primaria+A2319',PRIM+[STRESS])]: mk.append(dict(f=fg,sample=sn,I1=s.loc[S,'I1'].median()))
mk=pd.DataFrame(mk); C=pd.read_csv(OUT/'E2_contraste.csv'); out=[]
for r in C[(C.q=='I1')&C['sample'].isin(['primaria','primaria+A2319'])].itertuples():
    for fg in (0.10,0.12,0.14):
        x=mk[(mk.f==fg)&(mk['sample']==r.sample)].I1.values; o=r.med_cl
        out.append(dict(mass=r.mass,sample=r.sample,f=fg,p=min(1,2*min((1+(x<=o).sum())/(len(x)+1),(1+(x>=o).sum())/(len(x)+1)))))
o=pd.DataFrame(out).pivot_table(index=['mass','sample'],columns='f',values='p'); o.to_csv(OUT/'PH_mock_ruido_real.csv'); print('ruido',round(noise,3)); print(o.round(3))
# recalculo independente (sem common/diagnostics)
Gc=6.6743e-11; Ms=1.98847e30; kpc=3.08567758149e19; gd=299792458*(67.4e3/(kpc*1e3))/(2*np.pi); lam=np.exp(lnlam); rr=np.geomspace(100,1000,8)
I=[]
for c in PRIM:
    a=fits.getdata(f'/home/claude/xcop_in/{c}/{c}_fgas_profile.fits',1); R5=fits.getheader(f'/home/claude/xcop_in/{c}/{c}_fgas_profile.fits',1)['R500']
    s=fits.getdata(f'/home/claude/xcop_in/{c}/{c}_mstar.fits',2); h=fits.getdata(f'/home/claude/xcop_in/{c}/{c}_hydro_mass.fits',1)
    ll=lambda x,y: 10**np.interp(np.log10(rr),np.log10(np.array(x,float)),np.log10(np.array(y,float)))
    ok=np.array(s['MSTAR'])>0
    mb=ll(a['RADIUS']*R5,a['MGAS'])+ll(s['RADIUS'][ok],s['MSTAR'][ok]); mt=ll(h['RADIUS'],h['M_FORW'])
    x=Gc*mb*Ms/(rr*kpc)**2/gd; R0=np.log(mt/mb-1)+x-np.log(lam)+0.5*np.log(1-np.exp(-x)); I.append(np.polyfit(np.log(rr),R0,1)[0])
    print(c,round(I[-1],3),'logK',np.round(np.log10(mt/mb),2).tolist())
print('mediana I1 independente',round(float(np.median(I)),4),'| pipeline',round(float(C[(C['mass']=='M_FORW')&(C['sample']=='primaria')&(C.q=='I1')].med_cl.iloc[0]),4))
