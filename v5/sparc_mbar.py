import numpy as np
from scipy.optimize import least_squares
from scipy.stats import spearmanr
exec(open('sparc_test.py').read().split("print(\"N used\"")[0].replace('print("points"','#'))
G=6.674e-11; Msun=1.989e30; A=0.4068
gb=gbar(0.5,0.7); m=ok&(gb>0)
gals=sorted(set(ids)); Mb={}; xiG={}; rg={}
for g in gals:
    i=np.where((ids==g)&m)[0]
    if len(i)<5: continue
    j=i[np.argmax(R[i])]; Mb[g]=gb[j]*(R[j]*kpc)**2/G/Msun; xiG[g]=gb[j]/gd
    rg[g]=np.mean(np.log10(gb[i]*K(gb[i]/gd))-np.log10(gobs[i]))
gl=list(Mb); lM=np.log10([Mb[g] for g in gl]); res=np.array([rg[g] for g in gl]); xi=np.array([xiG[g] for g in gl])
print("galaxies",len(gl),"| logMbar range %.1f-%.1f"%(lM.min(),lM.max()))
print("Spearman(per-galaxy mean resid, log Mbar)=%+.2f ; vs xi_gal=%+.2f"%(spearmanr(res,lM)[0],spearmanr(res,xi)[0]))
for lo,hi in ((6,8),(8,9),(9,10),(10,11),(11,13)):
    k=(lM>=lo)&(lM<hi); print(" logMbar [%d,%d): N=%d mean resid %+.3f scatter %.3f"%(lo,hi,k.sum(),res[k].mean(),res[k].std()))
base=np.sqrt(np.mean((np.log10(gb[m]*K(gb[m]/gd))-np.log10(gobs[m]))**2))
xip=np.array([xiG.get(g,0) for g in ids]); Mp=np.array([Mb.get(g,np.nan) for g in ids]); m2=m&np.isfinite(Mp)
def rms(pred): return np.sqrt(np.mean((np.log10(pred[m2])-np.log10(gobs[m2]))**2))
print("baseline s=1,n=1/2: %.3f (same subset %.3f)"%(base,rms(gb*K(gb/gd))))
lm=np.log10(Mp)
# H1: n = 1/2 + A*xi_gal (paper n(xi)) with galaxy xi
print("H1 n=1/2+A*xi_gal (paper n(xi)): %.3f"%rms(gb*K(gb/gd,n=0.5+A*xip)))
# H2: amplitude lam scaled by Mbar power law
f2=least_squares(lambda q: np.log10((gb*K(gb/gd,n=0.5,lam=0.999)*0+gb*(1+q[0]*(lm-10)*0+ (K(gb/gd)-1)*10**(q[1]*(lm-10))))[m2])-np.log10(gobs[m2]),[0,-0.05]); print("H2 amplitude ~ Mbar^%.3f : %.3f"%(f2.x[1],np.sqrt(np.mean(f2.fun**2))))
# H3: g_dagger scaled by Mbar power law
f3=least_squares(lambda q: np.log10((gb*K(gb/(gd*10**(q[0]*(lm-10)))))[m2])-np.log10(gobs[m2]),[0.0]); print("H3 g_dagger ~ Mbar^%.3f : %.3f"%(f3.x[0],np.sqrt(np.mean(f3.fun**2))))
# H4: s = Mbar power law (1 param) and with x-dependence
f4=least_squares(lambda q: np.log10((gb*K(gb/gd,s=10**(q[0]*(lm-10))))[m2])-np.log10(gobs[m2]),[0.0]); print("H4 s ~ Mbar^%.3f : %.3f"%(f4.x[0],np.sqrt(np.mean(f4.fun**2))))
# dwarfs after best
for nm,pr in (("base",gb*K(gb/gd)),("H3",gb*K(gb/(gd*10**(f3.x[0]*(lm-10)))))):
    d=m2&(lm<9); b=m2&(lm>10.5); print(" %s: dwarfs(logM<9) mean %+.3f rms %.3f | massive(>10.5) mean %+.3f rms %.3f"%(nm,*(lambda r:(r.mean(),np.sqrt(np.mean(r**2))))(np.log10(pr[d])-np.log10(gobs[d])),*(lambda r:(r.mean(),np.sqrt(np.mean(r**2))))(np.log10(pr[b])-np.log10(gobs[b]))))
# cross-validation (leave-galaxy-out, 5 fold) for H3
rng=np.random.default_rng(1); fold={g:rng.integers(0,5) for g in gl}; fo=np.array([fold.get(g,-1) for g in ids]); cv=[]
for k in range(5):
    tr=m2&(fo!=k); te=m2&(fo==k)
    q=least_squares(lambda q: np.log10((gb*K(gb/(gd*10**(q[0]*(lm-10)))))[tr])-np.log10(gobs[tr]),[0.0]).x
    cv+=list(np.log10((gb*K(gb/(gd*10**(q[0]*(lm-10)))))[te])-np.log10(gobs[te]))
print("H3 5-fold CV rms %.3f vs base %.3f"%(np.sqrt(np.mean(np.square(cv))),rms(gb*K(gb/gd))))
# Upsilon sensitivity of dwarf bias
for Yd in (0.3,0.5,0.7):
    g_=gbar(Yd,0.7); mm=ok&(g_>0)&np.isfinite(Mp)&(lm<9); r=np.log10(g_[mm]*K(g_[mm]/gd))-np.log10(gobs[mm]); print(" Yd=%.1f dwarf mean resid %+.3f"%(Yd,r.mean()))
