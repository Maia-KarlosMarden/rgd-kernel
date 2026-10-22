import io,contextlib,numpy as np,warnings; warnings.filterwarnings("ignore")
from scipy.optimize import brentq
buf=io.StringIO()
with contextlib.redirect_stdout(buf): exec(open('stars_test.py').read().split('print("\\nALLEN')[0])
exec(open('allen_data.py').read()); exec(open('mantz_data.py').read())
Gc=6.674e-11; cc=3e8; Ms=1.989e30; kpc=3.0857e19; Mpc=3.0857e22; lam_=.999; A_=.4068
H=lambda z: 70e3/Mpc*np.sqrt(.3*(1+z)**3+.7)
rcz=lambda z: 3*H(z)**2/(8*np.pi*Gc)
# CLASH (Umetsu+2016): name, z, M2500 (1e14), e, M500, e
CL={"A383":(.187,2.78,.63,5.88,1.73),"RXJ2129":(.234,2.08,.44,4.48,1.16),"A611":(.288,4.13,.92,10.73,2.65),"MS2137":(.313,2.47,.72,8.28,2.57),
"MACS1931":(.352,4.10,1.12,10.51,4.05),"MACS1115":(.352,3.50,.82,10.67,2.22),"MACS1532":(.363,1.83,1.01,4.17,1.71),"MACS1720":(.391,3.92,.86,9.96,2.53),
"MACS0429":(.399,2.85,.70,6.85,2.10),"MACS0329":(.450,3.29,.64,6.51,1.37),"RXJ1347":(.451,7.65,1.63,22.33,4.89),"MACS0744":(.686,4.32,1.02,11.94,2.81)}
AL={"A383":"A383","RXJ2129":"RXJ2129.6+0005","A611":"A611","MS2137":"MS2137.3-2353","MACS1931":"MACS1931.8-2635","MACS1115":"MACS1115.8+0129","MACS1532":"MACS1532.9+3021","MACS1720":"MACS1720.3+3536","MACS0429":"MACS0429.6-0253","MACS0329":"MACS0329.7-0212","RXJ1347":"RXJ1347.5-1144","MACS0744":"MACS0744.9+3927"}
Ad={a[0]:a for a in ALLEN}
def solve(m,r,xi,gd,s):
    g=Gc*m/r**2; return brentq(lambda f: f*(1+lam_*np.exp(-f*g/gd/s)/(-np.expm1(-f*g/gd/s))**(.5+A_*xi))-1,1e-4,5.0)
PHI=.7; ALPHA=1.2
out=[]
print("%-9s %5s | M2500: lens   X-ray  L/X | fgas(X)  f_bar obs(X-mass) obs(lens) | pred p=.5 | obs/pred X  lens"%("cluster","z"))
for k,(z,ML,eML,M5L,eM5) in CL.items():
    a=Ad[AL[k]]; rX=a[2]*kpc; MX=2500*rcz(z)*4*np.pi/3*rX**3/Ms/1e14; fg=a[3]
    ML_=ML*1e14*Ms; rL=(ML_/(2500*rcz(z)*4*np.pi/3))**(1/3); R5=(M5L*1e14*Ms/(500*rcz(z)*4*np.pi/3))**(1/3)
    gd=cc*H(z)/(2*np.pi); xi=Gc*M5L*1e14*Ms/R5**2/gd
    stX=PHI*Mstar_tot(M5L*1e14)*Ms/(fg*MX*1e14*Ms); fbX=fg*(1+stX)
    MgasX=fg*MX*1e14*Ms; MgasL=MgasX*(rL/rX)**ALPHA; stars=PHI*Mstar_tot(M5L*1e14)*Ms; fbL=(MgasL+stars)/ML_
    fbX=(MgasX+stars)/(MX*1e14*Ms)
    pr=solve(ML_,rL,xi,gd,5.0**.5); prX=solve(MX*1e14*Ms,rX,xi,gd,5.0**.5)
    out.append((k,ML/MX,eML/ML,np.log(fbX/prX),np.log(fbL/pr),np.log10(ML*1e14)))
    print("%-9s %5.3f |   %5.2f  %5.2f  %4.2f ±%.2f |  %.3f   %.3f   %.3f | %.3f | %.2f  %.2f"%(k,z,ML,MX,ML/MX,eML/MX,fg,fbX,fbL,pr,fbX/prX,fbL/pr))
o=np.array([x[1:] for x in out]); r=o[:,0]; w=1/o[:,1]**2/r**2
print("\nN=%d  M_L/M_X(2500): median %.2f ; mean ln-weighted %.2f ; rms scatter(ln) %.2f"%(len(r),np.median(r),np.exp(np.sum(np.log(r)/ (o[:,1])**2)/np.sum(1/o[:,1]**2)),np.std(np.log(r))))
print("Without lensing (X-ray masses, ξ from lens M500): median obs/pred = %.2f ; with lens: %.2f ; rms ln: %.2f vs %.2f"%(np.exp(np.median(o[:,2])),np.exp(np.median(o[:,3])),np.sqrt(np.mean(o[:,2]**2)),np.sqrt(np.mean(o[:,3]**2))))
lm=o[:,4]; print("trend of ln(M_L/M_X) with log M_L2500: %+.2f ; of ln(obs/pred): X %+.2f , lens %+.2f"%(np.polyfit(lm,np.log(r),1)[0],np.polyfit(lm,o[:,2],1)[0],np.polyfit(lm,o[:,3],1)[0]))
