import io,contextlib,numpy as np,warnings; warnings.filterwarnings("ignore")
from scipy.optimize import brentq
buf=io.StringIO()
with contextlib.redirect_stdout(buf):
    exec(open('stars_test.py').read().split('print("\\nALLEN')[0])
    exec(open('vik_test.py').read().split('PHI=.7')[0].split('buf=io.StringIO()')[1].split("\n",1)[1]) if False else None
Gc=6.674e-11; cc=2.998e8; Ms=1.989e30; kpc=3.0857e19; Mpc=3.0857e22; lam_=.999; A_=.4068
def solve(m,r,xi,gd,mode,w=None,dd=None,p=.5):
    g=Gc*m/r**2; al=.355*(1-np.exp(-xi/.052))
    def h(f):
        x=f*g/gd
        if mode=="v4": n=(.5+A_*xi)*w**(-al); s=1.
        elif mode=="v5": n=.5+A_*xi; s=dd**p
        elif mode=="none": n=.5+A_*xi; s=1.
        return f*(1+lam_*np.exp(-x/s)/(-np.expm1(-x/s))**n)-1
    try: return brentq(h,1e-6,5.0)
    except ValueError: return np.nan
# ---- Vikhlinin 2006
Hv=lambda z: 71e3/Mpc*np.sqrt(.3*(1+z)**3+.7)
V=[("A133",.0569,998,1.24,3.14,.065,.084),("A383",.1883,956,1.68,3.10,.090,.122),("A478",.0881,1359,4.23,7.83,.096,.118),("A907",.1603,1117,2.30,4.71,.090,.121),
("A1413",.1429,1339,3.08,7.78,.092,.104),("A1795",.0622,1283,2.66,6.57,.088,.098),("A1991",.0592,753,.60,1.28,.068,.099),("A2029",.0779,1380,4.41,8.29,.090,.121),
("A2390",.2302,1448,3.50,10.88,.127,.138),("MKW4",.0199,635,.30,.74,.044,.063)]
res={m:[] for m in ("none","v4","v5")}; rat={m:[] for m in res}; ro=[]
for n,z,r5,M25,M5,f25,f5 in V:
    rc=3*Hv(z)**2/(8*np.pi*Gc); R5=r5*kpc; M5k=M5*1e14*Ms; M25k=M25*1e14*Ms; r25=(M25k/(2500*rc*4*np.pi/3))**(1/3)
    gd=cc*Hv(z)/(2*np.pi); xi=Gc*M5k/R5**2/gd; st=Mstar_tot(M5*1e14)*Ms
    fb5=f5+st/M5k; fb25=f25+.7*st/M25k; ro.append(fb25/fb5)
    p5=solve(M5k,R5,xi,gd,"none")
    for m in res:
        p25=solve(M25k,r25,xi,gd,m,w=r25/R5,dd=M25k/M5k/(r25/R5)**3); res[m].append(np.log(fb25/p25)); rat[m].append(p25/p5)
print("INDEPENDENT test 1: Vikhlinin+2006, 10 clusters, r2500 (w=r2500/r500 = 0.39-0.48; Delta/Delta500=5)")
for m,lab in (("none","no radial term (n(xi) at w=1)"),("v4","v4 radial law n=A(xi)w^-alpha(xi)"),("v5","v5 s=(Delta/Delta500)^0.5")):
    a=np.array(res[m]); ok=np.isfinite(a); print("  %-36s no-solution %d/10 | obs/pred median %.2f  rms(ln) %.3f | f_bar(r2500)/f_bar(r500): pred %.2f (obs %.2f)"%(lab,(~ok).sum(),np.exp(np.nanmedian(a)),np.sqrt(np.nanmean(a**2)),np.nanmedian(rat[m]),np.median(ro)))
# ---- X-COP Ettori 2019 r500 -> r200
H7=lambda z: 70e3/Mpc*np.sqrt(.3*(1+z)**3+.7)
X=[("A1644",.0473,3.48,6.69,.128,.126),("A1795",.0622,4.63,6.53,.139,.144),("A2029",.0773,8.65,12.25,.141,.152),("A2142",.0909,8.95,13.64,.158,.168),("A2255",.0809,5.26,10.33,.153,.146),
("A3158",.0597,4.26,6.63,.145,.155),("A644",.0704,5.66,7.67,.132,.139),("A85",.0555,5.65,8.50,.150,.159),("RXC1825",.0650,4.08,6.15,.133,.155)]
res={m:[] for m in ("none","v4","v5")}; rat={m:[] for m in res}; ro=[]
for n,z,M5,M2,f5,f2 in X:
    rc=3*H7(z)**2/(8*np.pi*Gc); M5k=M5*1e14*Ms; M2k=M2*1e14*Ms
    R5=(M5k/(500*rc*4*np.pi/3))**(1/3); R2=(M2k/(200*rc*4*np.pi/3))**(1/3); gd=cc*H7(z)/(2*np.pi); xi=Gc*M5k/R5**2/gd; st=Mstar_tot(M5*1e14)*Ms
    fb5=f5+st/M5k; fb2=f2+st/M2k; ro.append(fb2/fb5); p5=solve(M5k,R5,xi,gd,"none")
    for m in res:
        p2=solve(M2k,R2,xi,gd,m,w=R2/R5,dd=M2k/M5k/(R2/R5)**3); res[m].append(np.log(fb2/p2)); rat[m].append(p2/p5)
print("\nINDEPENDENT test 2: X-COP (Ettori 2019 Table 2: measured f_gas at r500 and r200; 9 clusters excl. mergers/outliers), w=1.5")
for m,lab in (("none","no radial term (n(xi) at w=1)"),("v4","v4 radial law n=A(xi)w^-alpha(xi)"),("v5","v5 s=(Delta/Delta500)^0.5")):
    r=np.array(rat[m]); print("  %-36s no-solution %d/9 | f_bar(r200)/f_bar(r500): pred %.2f (obs %.2f) | rms(ln) of obs/pred ratio %.3f"%(lab,(~np.isfinite(r)).sum(),np.nanmedian(r),np.median(ro),np.sqrt(np.nanmean(np.log(np.array(ro)/r)**2))))
