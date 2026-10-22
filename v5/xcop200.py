import io,contextlib,numpy as np,warnings; warnings.filterwarnings("ignore")
from scipy.optimize import brentq
buf=io.StringIO()
with contextlib.redirect_stdout(buf): exec(open('stars_test.py').read().split('print("\\nALLEN')[0])
Gc=6.674e-11; cc=2.998e8; Ms=1.989e30; Mpc=3.0857e22; lam_=.999; A_=.4068
H=lambda z: 70e3/Mpc*np.sqrt(.3*(1+z)**3+.7)
# Ettori+2019 Table 2: name,z, M500,M200 (1e14, HSE), fgas500, fgas200 ; z from Eckert/our table
X=[("A1644",.0473,3.48,6.69,.128,.126),("A1795",.0622,4.63,6.53,.139,.144),("A2029",.0773,8.65,12.25,.141,.152),("A2142",.0909,8.95,13.64,.158,.168),
("A2255",.0809,5.26,10.33,.153,.146),("A2319",.0557,7.31,10.18,.189,.237),("A3158",.0597,4.26,6.63,.145,.155),("A3266",.0589,8.80,15.12,.132,.108),
("A644",.0704,5.66,7.67,.132,.139),("A85",.0555,5.65,8.50,.150,.159),("RXC1825",.0650,4.08,6.15,.133,.155),("ZW1215",.0766,7.66,13.03,.106,.092)]
def solve(m,r,xi,gd,s):
    g=Gc*m/r**2; return brentq(lambda f: f*(1+lam_*np.exp(-f*g/gd/s)/(-np.expm1(-f*g/gd/s))**(.5+A_*xi))-1,1e-4,5.0)
models={"p=0.50":lambda x:.5,"p=0.46":lambda x:.46,"p(xi)=0.46(xi/0.4)^0.37":lambda x:.46*(x/.4)**.37}
PHI=1.0
print("%-8s %5s %5s | f_bar obs: 500  200 (ratio) | pred ratio r200/r500: [p=.5 / p=.46 / p(xi)] | obs/pred at 200"%("name","xi","r2/r5"))
obs_r=[];pr_r={k:[] for k in models};o200={k:[] for k in models};o500=[];lm=[]
for n,z,M5,M2,f5,f2 in X:
    rc=3*H(z)**2/(8*np.pi*Gc); M5k=M5*1e14*Ms; M2k=M2*1e14*Ms
    R5=(M5k/(500*rc*4*np.pi/3))**(1/3); R2=(M2k/(200*rc*4*np.pi/3))**(1/3)
    gd=cc*H(z)/(2*np.pi); xi=Gc*M5k/R5**2/gd
    st=PHI*Mstar_tot(M5*1e14)*Ms          # total stellar mass, all inside r500 (Kravtsov)
    fb5=f5+st/M5k; fb2=f2+st/M2k
    pf5=solve(M5k,R5,xi,gd,1.0); o500.append(np.log(fb5/pf5)); obs_r.append(fb2/fb5)
    line=[]; 
    for k,pf in models.items():
        p=pf(xi); pf2=solve(M2k,R2,xi,gd,(0.4)**p); pr_r[k].append(pf2/pf5); o200[k].append(np.log(fb2/pf2)); line.append(pf2/pf5)
    lm.append(np.log10(M5*1e14))
    print("%-8s %5.2f %5.2f | %.3f %.3f (%.2f) | %.2f %.2f %.2f | %.2f %.2f %.2f"%(n,xi,R2/R5,fb5,fb2,fb2/fb5,*line,*[np.exp(o200[k][-1]) for k in models]))
lm=np.array(lm); oR=np.array(obs_r)
print("\nr500 control (s=1): obs/pred median %.2f rms(ln) %.3f"%(np.exp(np.median(o500)),np.sqrt(np.mean(np.square(o500)))))
for k in models:
    a=np.array(o200[k]); pr=np.array(pr_r[k])
    print("%-26s r200: obs/pred median %.2f rms(ln) %.3f trend/dex %+.2f | ratio r200/r500 pred median %.2f (obs %.2f) ; rms(ln obs/pred ratio) %.3f ; ratio-trend pred %+.2f obs %+.2f"%(k,np.exp(np.median(a)),np.sqrt(np.mean(a**2)),np.polyfit(lm,a,1)[0],np.median(pr),np.median(oR),np.sqrt(np.mean(np.log(oR/pr)**2)),np.polyfit(lm,np.log(pr),1)[0],np.polyfit(lm,np.log(oR),1)[0]))
print("\nWith A2319 (merger, f_gas200=0.237) and A3266/ZW1215 excluded:")
mask=np.array([n not in ("A2319","A3266","ZW1215") for n,*_ in X])
for k in models:
    pr=np.array(pr_r[k])[mask]; print("  %-26s pred median ratio %.2f ; obs median %.2f ; rms(ln) %.3f"%(k,np.median(pr),np.median(oR[mask]),np.sqrt(np.mean(np.log(oR[mask]/pr)**2))))
