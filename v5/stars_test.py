import io,contextlib,numpy as np
from scipy.optimize import minimize_scalar, brentq
buf=io.StringIO()
with contextlib.redirect_stdout(buf): exec(open('groups_test.py').read())
from allen_data import ALLEN
def Mstar_tot(M500_Msun):      # Kravtsov+2018 Table 2, 21 clusters (this work + G13): log M*,tot = 12.71 + 0.59 (log M500 - 14.5), scatter 0.11 dex, Chabrier IMF
    return 10**(12.71+0.59*(np.log10(M500_Msun)-14.5))
def Mstar_gonz(M500_Msun):     # Gonzalez+2013 Table 4: M*,3D/Mgas = 0.39 (M500/1e14)^-0.84 with Mgas = 8.8e-2 (M500/1e14)^1.26 *1e14  [r<r500]
    return 0.39*(M500_Msun/1e14)**-0.84*8.8e-2*(M500_Msun/1e14)**1.26*1e14
print("Stellar mass within r500 for a few halo masses (Msun) : Kravtsov / Gonzalez")
for M in (1.5e13,3e13,1e14,3e14,1e15): print("  M500=%.1e : M*tot = %.2e / %.2e  (M*/M500 = %.3f / %.3f)"%(M,Mstar_tot(M),Mstar_gonz(M),Mstar_tot(M)/M,Mstar_gonz(M)/M))
print("(Kravtsov relation fitted for M500 ~ 4.6e13-1.6e15; groups below ~5e13 are an extrapolation)\n")
def gp(phi,src="K"):
    st=np.array([phi*(Mstar_tot(r['M5']/Msun) if src=="K" else Mstar_gonz(r['M5']/Msun))*Msun/(r['f25']*r['M25']) for r in grows]); return st
print("GROUPS (43): stellar-to-gas ratio inside r2500 implied by the literature relation, M*(<r2500)=phi*M*tot(M500), Mgas(<r2500)=f_gas2500*M2500")
for phi in (.5,.7,1.0):
    for src in ("K","G"):
        st=gp(phi,src)
        def cost(p): rr=np.array([gresid(r,st[i],1,p=p) for i,r in enumerate(grows)]); return np.mean(rr**2)
        f=minimize_scalar(cost,bounds=(.1,1.5),method="bounded"); p=f.x; rr=np.array([gresid(r,st[i],1,p=p) for i,r in enumerate(grows)])
        r5=np.array([gresid(r,st[i],1,p=.46) for i,r in enumerate(grows)]); r5b=np.array([gresid(r,st[i],1,p=.5) for i,r in enumerate(grows)])
        xi=np.array([r['xi'] for r in grows]); lo=xi<np.median(xi)
        print(" phi=%.1f %s: M*/Mgas(<r2500) median %.2f [%.2f-%.2f] | best constant p=%.2f (rms %.3f) | at p=0.46: rms %.3f mean %+.3f | p=0.5: mean %+.3f | residual(p=.46) low-xi half %+.3f high-xi half %+.3f"%(phi,src,np.median(st),st.min(),st.max(),p,np.sqrt(f.fun),np.sqrt(np.mean(r5**2)),r5.mean(),r5b.mean(),r5[lo].mean(),r5[~lo].mean()))
print("\nALLEN+2008 clusters (42): same stars treatment, p=0.46 and 0.5")
G_=6.674e-11; c_=3e8; kpc_=3.0857e19; Mpc_=3.0857e22
Hz=lambda z: 70e3/Mpc_*np.sqrt(.3*(1+z)**3+.7)
mf=lambda x: np.log(1+x)-x/(1+x)
z=np.array([a[1] for a in ALLEN]); r25=np.array([a[2] for a in ALLEN])*kpc_; fo=np.array([a[3] for a in ALLEN])
M25=np.array([2500*3*Hz(zz)**2/(8*np.pi*G_)*4*np.pi/3*r**3 for zz,r in zip(z,r25)])
M5=[];xi=[]
for zz,r,m in zip(z,r25,M25):
    rc_=3*Hz(zz)**2/(8*np.pi*G_); f=lambda l: 500*rc_*4*np.pi/3*(10**l)**3*mf(3.0*r/10**l)/mf(3.0)-m
    R5=10**brentq(f,np.log10(r*1.2),np.log10(r*6)); M5.append(500*rc_*4*np.pi/3*R5**3); xi.append(G_*M5[-1]/R5**2/(c_*Hz(zz)/(2*np.pi)))
M5=np.array(M5); xi=np.array(xi); lM=np.log10(M25/Msun)
def fsolve(m,r,x,gd,p):
    g=G_*m/r**2; s=5.0**p; K=lambda f: 1+lam*np.exp(-f*g/gd/s)/(-np.expm1(-f*g/gd/s))**(.5+A*x); return brentq(lambda f: f*K(f)-1,1e-4,5.0)
for phi in (.5,.7,1.0):
    st=phi*Mstar_tot(M5/Msun)*Msun/(fo*M25)
    for p in (.46,.5):
        pr=np.array([fsolve(m,r,x,c_*Hz(zz)/(2*np.pi),p) for m,r,x,zz in zip(M25,r25,xi,z)]); fb=fo*(1+st)
        print(" phi=%.1f p=%.2f: M*/Mgas(<r2500) median %.2f (%.2f-%.2f) | obs/pred f_bar median %.2f | mass slope obs %+.2f RGD %+.2f | z slope obs %+.2f RGD %+.2f"%(phi,p,np.median(st),st.min(),st.max(),np.median(fb/pr),np.polyfit(lM,np.log(fb),1)[0],np.polyfit(lM,np.log(pr),1)[0],np.polyfit(np.log(1+z),np.log(fb),1)[0],np.polyfit(np.log(1+z),np.log(pr),1)[0]))
