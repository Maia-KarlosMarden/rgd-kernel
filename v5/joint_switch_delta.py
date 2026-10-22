exec(open('joint_switch.py').read().split("Mhalo_pt={}")[0])
print("SPARC Delta/Delta500: min %.1f, 1%% %.1f, 5%% %.1f ; fraction of points with Delta<20: %.3f, <10: %.4f"%(dens[mg].min(),*np.percentile(dens[mg],[1,5]),np.mean(dens[mg]<20),np.mean(dens[mg]<10)))
ld_c=np.log(densc); ld_s=np.log(dens)
print("Switch on Delta(r) itself (term active only where Delta<Delta_c), no mass/xi dependence:")
for mm in (2,4,8,32):
    best=(1e99,)
    for l0 in np.linspace(np.log(1),np.log(60),300):
        v=m2ll(rc(expit(mm*(l0-ld_c)),.5),rs(expit(mm*(l0-ld_s)),.5))
        if v<best[0]: best=(v,l0)
    l0=best[1]; print("  m=%2d: -2lnL=%7.0f (dAIC vs all-off %+.0f) Delta_c=%.1f | rms clusters %.4f SPARC %.4f"%(mm,best[0],best[0]+4-base,np.exp(l0),np.sqrt(np.mean(rc(expit(mm*(l0-ld_c)),.5)**2)),np.sqrt(np.mean(rs(expit(mm*(l0-ld_s)),.5)**2))))
print("Reference xi-switch best: -2lnL -11139 (m=2)")
