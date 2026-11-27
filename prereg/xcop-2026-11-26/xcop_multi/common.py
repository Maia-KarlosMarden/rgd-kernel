import sys, glob, numpy as np, pandas as pd
from pathlib import Path
from scipy.integrate import cumulative_trapezoid
from astropy.io import fits
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from src import fit as F, mocks as M, diagnostics as D, features
from src.constants import G, MSUN, KPC, rhoc, gdag
DATA=Path('/home/claude/xcop_in'); OUT=ROOT/'runs'/'xcop_multi'
PRIM=['A85','A644','A1795','A2029','A2142','ZW1215']; STRESS='A2319'
ALL=['A1644','A1795','A2029','A2142','A2255','A2319','A3158','A3266','A644','A85','RXC1825','ZW1215']
RAD=np.geomspace(100,1000,8)
LAWS=['RGD','RAR','SPL','F3','G1','F3+G','N2','F3+G+N']
def one(c,pat):
    g=glob.glob(str(DATA/c/f'*{pat}*.fits')); assert len(g)==1,(c,pat,g); return g[0]
def loglog(r,rr,y): return np.exp(np.interp(np.log(r),np.log(rr),np.log(y)))
def baryons(c,fstar_missing=0.09):
    """SO barions. Le apenas RADIUS e MGAS do arquivo de gas (nunca M_NFW/FGAS) e MSTAR_SMOOTHED."""
    with fits.open(one(c,'fgas_profile')) as h:
        t=h[1].data; R500=float(h[1].header['R500']); z=float(h[1].header['REDSHIFT'])
        r=np.array(t['RADIUS'],float)*R500; mg=np.array(t['MGAS'],float)
    assert np.all(np.diff(r)>0) and np.all(mg>0)
    if glob.glob(str(DATA/c/'*_mstar.fits')):
        with fits.open(one(c,'_mstar')) as h:
            s=h['MSTAR_SMOOTHED'].data; rs=np.array(s['RADIUS'],float); ms=np.array(s['MSTAR'],float)
        ok=(ms>0)&(rs>0); star=loglog(r,rs[ok],ms[ok]); measured=True
        assert rs[ok].min()<=100 and rs[ok].max()>=1000
    else: star=fstar_missing*mg; measured=False
    assert r.min()<=100 and r.max()>=1000
    return dict(c=c,z=z,r=r,mgas=mg,mstar=star,measured=measured)
def profile(b,scaled_stars=None):
    """quantidades nos 8 raios; Z, P-, P+ calculados no perfil completo."""
    r=b['r']; star=b['mstar'] if scaled_stars is None else scaled_stars*b['mgas']
    mb=b['mgas']+star; rm=r*KPC; gb=G*mb*MSUN/rm**2
    S=np.gradient(np.log(gb),np.log(r)); I=cumulative_trapezoid(gb,rm,initial=0)
    Z=(1-S)/3; Pm=I/(rm*gb); Pp=(I[-1]-I)/(rm*gb)
    at=lambda v:np.interp(np.log(RAD),np.log(r),v)
    gbR=G*loglog(RAD,r,mb)*MSUN/(RAD*KPC)**2
    return pd.DataFrame(dict(catalog='XCOPm',system_id=b['c'],r_kpc=RAD,radius_m=RAD*KPC,gbar_ms2=gbR,x=gbR/gdag(),
        Mb=loglog(RAD,r,mb),Mgas=loglog(RAD,r,b['mgas']),Mstar=loglog(RAD,r,np.maximum(star,1e-30)),Zeq=at(Z),Pm=at(Pm),Pp=at(Pp),Rd_kpc=np.nan,Sb=at(S)))
def sparc():
    sp=M.Base().sp
    return sp
def fit_laws(sp): return {m:F.fit(sp,m) for m in LAWS}
def slopes(df,y,lnlam,min_pts=4):
    return D.slopes(df,y,lnlam,min_pts=min_pts).set_index('system')
