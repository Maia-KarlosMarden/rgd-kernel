"""PRE-REGISTERED (written before running).
Data: XXL-365-GC (Adami+18, VizieR IX/52/xxl365gc): Mgas(<500 kpc), r500MT (from M-T relation calibrated on weak lensing, Lieu+16), z, class.
Sample: systems with Mgas AND r500MT. Primary subsample: Class 1 (C1, well-defined selection). Secondary: C1+C2.
M500 = 4pi/3 *500 rho_c(z) r500MT^3 (H0=70, Om=0.28). Total mass inside 500 kpc: NFW with c200 from Dutton&Maccio14, solved so that M(<r500)=M500.
Prediction (archived prescription: lam=.999, A=.4068, p=.5, xi_c=.115, m=4, Kravtsov stars phi=1): f_bar(<500 kpc) from predict_frozen.predict_fbar(M(<500kpc), 0.5 Mpc, M500, r500, z);
  f_gas,pred = f_bar - M*_tot(M500)/M(<500kpc).   Observed f_gas = Mgas(500kpc)/M(<500kpc).
Metric: median and bootstrap error of log10(f_gas,pred/f_gas,obs) in log-M500 bins (edges 1e13,3e13,1e14,3e14,1e15), plus overall.
Sensitivities (declared in advance): p=0.8; stars from HSC-XXL measured relation instead of Kravtsov; excluding flagged 'F'.
No parameter is tuned on these data.
"""
import sys, numpy as np
sys.path.insert(0,'.')
import blind_xxl_efeds as bt
import predict_frozen as pf
from scipy.optimize import brentq
G=pf.G; Mpc=pf.Mpc; Msun=pf.Msun
rows=[]
for ln in open('xxl365gc_vizier.txt'):
    f=ln.rstrip('\n').split('\t')
    if len(f)>=19 and f[0].strip().isdigit():
        def num(s):
            s=s.strip(); 
            try: return float(s)
            except: return np.nan
        rows.append(dict(id=int(f[1]),z=num(f[4]),cls=int(f[6]),Mg=num(f[7])*1e11,eMg=0.5*(num(f[8])+num(f[9])),r500=num(f[10]),flag=f[18].strip()))
print("rows parsed:",len(rows))
def Hz2(z,H0=70.,Om=0.28): return (H0*1e3/Mpc)*np.sqrt(Om*(1+z)**3+1-Om)
def M500_from_r(r_kpc,z):
    rho=3*Hz2(z)**2/(8*np.pi*G); r=r_kpc*1e-3*Mpc
    return 4*np.pi/3*500*rho*r**3/Msun
nfw=lambda x: np.log(1+x)-x/(1+x)
def c200_DM14(M200,z,h=.7):
    a=.520+(.905-.520)*np.exp(-.617*z**1.21); b=-.101+.026*z
    return 10**(a+b*np.log10(M200/(1e12/h)))
def M_in(M500,r500_mpc,z,r_mpc=0.5):
    rho=3*Hz2(z)**2/(8*np.pi*G)
    def f(lM200):
        M200=10**lM200; R200=(3*M200*Msun/(4*np.pi*200*rho))**(1/3)/Mpc
        c=c200_DM14(M200,z); return M200*nfw(c*r500_mpc/R200)/nfw(c)-M500
    lM=brentq(f,np.log10(M500)+0.01,np.log10(M500)+1.0)
    M200=10**lM; R200=(3*M200*Msun/(4*np.pi*200*rho))**(1/3)/Mpc; c=c200_DM14(M200,z)
    return M200*nfw(c*r_mpc/R200)/nfw(c)
out=[]
for d in rows:
    if not (np.isfinite(d['Mg']) and np.isfinite(d['r500']) and d['Mg']>0): continue
    z=d['z']; M5=M500_from_r(d['r500'],z); r5=d['r500']*1e-3
    try: Mr=M_in(M5,r5,z)
    except Exception: continue
    d.update(M5=M5,Mr=Mr,r5=r5); out.append(d)
print("usable:",len(out)," C1:",sum(d['cls']==1 for d in out))
def pred(d,p=.5,stars='krav'):
    fb=pf.predict_fbar(d['Mr'],0.5,d['M5'],d['r5'],d['z'],p=p)[0]
    st=pf.Mstar_tot(d['M5']) if stars=='krav' else bt.xxl_frac(d['M5'],d['z'])[1]*d['M5']
    return fb-st/d['Mr']
def run(sub,label,**kw):
    lr=[];M=[]
    for d in sub:
        fo=d['Mg']/d['Mr']; fp=pred(d,**kw)
        if fp>0 and fo>0: lr.append(np.log10(fp/fo)); M.append(d['M5'])
        # fp<=0 counted separately
    lr=np.array(lr);M=np.array(M); neg=sum(1 for d in sub if pred(d,**kw)<=0)
    edges=[1e13,3e13,1e14,3e14,1e15]; rng=np.random.default_rng(3)
    print("\n-- %s  %s : N=%d (pred<=0 for %d)  overall median %+.3f dex (x%.2f)"%(label,kw,len(lr),neg,np.median(lr),10**np.median(lr)))
    for lo,hi in zip(edges[:-1],edges[1:]):
        s=(M>=lo)&(M<hi)
        if s.sum()<3: print("   %.0e-%.0e  N=%d (too few)"%(lo,hi,s.sum())); continue
        bs=[np.median(rng.choice(lr[s],s.sum())) for _ in range(2000)]
        print("   %.0e-%.0e  N=%3d  median %+.3f ± %.3f dex  (x%.2f)  rms %.3f"%(lo,hi,s.sum(),np.median(lr[s]),np.std(bs),10**np.median(lr[s]),lr[s].std()))
C1=[d for d in out if d['cls']==1]
run(C1,"C1, archived")
run(out,"C1+C2, archived")
run(C1,"C1, p=0.8",p=0.8)
run(C1,"C1, HSC-XXL stars",stars='xxl')
run([d for d in C1 if 'F' not in d['flag']],"C1 without F flag, archived")
