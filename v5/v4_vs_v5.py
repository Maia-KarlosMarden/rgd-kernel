import io,contextlib,numpy as np,warnings; warnings.filterwarnings("ignore")
from scipy.optimize import minimize_scalar
buf=io.StringIO()
src=open('kernel_scenarios.py').read().split("# ---------- check NFW R500")[0]
with contextlib.redirect_stdout(buf): exec(src)
D=build(radii_mpc=[.5,1,1.5])          # paper mode: f_bar = 1.09*f_gas(M500) (Sun+2009 scaling), NFW masses, 13 clusters x 3 radii
x=D['fb']*D['g']/D['gd']; xi=D['xi']; w=D['w']
Mr=D['Mr']; M5=D['M5']; dens=(Mr/M5)/w**3            # Delta/Delta500 = rhobar(<r)/(500 rho_c)
def K_(n,s): xs=x/s; return 1+lam*np.exp(-xs)/(-np.expm1(-xs))**n
r_=lambda K: np.log10(K*D['fb'])
R=lambda a: float(np.sqrt(np.mean(np.asarray(a)**2)))
rw=lambda a,j: R(np.asarray(a)[:,j])
alpha=.355*(1-np.exp(-xi/.052))
res={}
res["v4 n(xi) fixed at w=1 (no radial term)"]=r_(K_(.5+A*xi,1.))
res["v4 n(w,xi)=1/2+A xi w^-alpha(xi)  [2 extra params]"]=r_(K_(.5+A*xi*w**(-alpha),1.))
f=lambda p: R(r_(K_(.5+A*xi,dens**p)))
o=minimize_scalar(f,bounds=(0,2),method="bounded"); p=o.x
res["v5 s=(Delta/Delta500)^p, p=%.2f fitted  [1 extra param]"%p]=r_(K_(.5+A*xi,dens**p))
res["v5 s=(Delta/Delta500)^0.5, p fixed (value used elsewhere)"]=r_(K_(.5+A*xi,dens**.5))
res["Newton (f_bar only)"]=np.log10(D['fb'])
print("X-COP 39 points (13 clusters x r=0.5,1.0,1.5 Mpc), same 'paper mode' data and units as v4 Table 4 (dex)")
print("%-62s %7s | per radius 0.5 / 1.0 / 1.5 Mpc | mean"%("model","RMS"))
for k,v in res.items(): print("%-62s %7.3f | %.3f / %.3f / %.3f | %+.3f"%(k,R(v),rw(v,0),rw(v,1),rw(v,2),np.mean(v)))
print("\nDelta/Delta500 at the three radii: median %s ; w range %.2f-%.2f"%(np.round(np.median(dens,axis=0),2),w.min(),w.max()))
# leave-one-cluster-out predictive for p (v5), and number of parameters
pl=[]
for i in range(N):
    m=np.ones(N,bool); m[i]=False
    g=lambda p: R(r_(K_(.5+A*xi[m],dens[m]**p)) if False else np.log10(Kb(m,p)))
    def Kb(m,p):
        xs=x[m]/dens[m]**p; return (1+lam*np.exp(-xs)/(-np.expm1(-xs))**(.5+A*xi[m]))*D['fb'][m]
    pl.append(minimize_scalar(lambda p: R(np.log10(Kb(m,p))),bounds=(0,2),method="bounded").x)
loo=[R(np.log10((1+lam*np.exp(-x[i]/dens[i]**pl[i])/(-np.expm1(-x[i]/dens[i]**pl[i]))**(.5+A*xi[i]))*D['fb'][i])) for i in range(N)]
print("leave-one-cluster-out p for v5: range %.2f-%.2f ; predictive RMS (dex) = %.3f"%(min(pl),max(pl),np.sqrt(np.mean(np.square(loo)))))
