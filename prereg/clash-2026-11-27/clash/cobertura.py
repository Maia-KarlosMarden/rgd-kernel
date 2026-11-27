"""Cobertura e massa de gás a partir dos perfis desprojetados ACCEPT2.0 (somente bárions)."""
import glob, json, os, re, sys
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

D, XC, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
MSUN, MP, KPC, MUE = 1.98892e33, 1.6726e-24, 3.0857e21, 1.17
CLASH = {"ABELL_0209": .206, "ABELL_0383": .1871, "ABELL_0611": .288, "ABELL_2261": .224,
         "ABELL_S1063": .354, "MCXC_J2129.6+0005": .235, "MACS_J0329.6-0211": .450,
         "MACS_J0416.1-2403": .396, "MACS_J0429.6-0253": .399, "MACS_J0647.8+7015": .584,
         "MACS_J0717+3745": .548, "MACS_J0744.9+3927": .686, "MACS_J1115.8+0129": .352,
         "MACS_J1149.6+2223": .544, "MACS_J1206.2-0847": .440, "MACS_J1311.0-0311": .494,
         "MACS_J1720.2+3536": .3913, "MACS_J1931.8-2635": .352, "MACS_J2140.2-2339": .313, "LCDCS_0829": .451}
XCOP = {"ABELL_0085": "A85", "ABELL_0644": "A644", "ABELL_1795": "A1795", "ABELL_2029": "A2029",
        "ABELL_2142": "A2142", "ABELL_2319": "A2319", "NSC_J121733+033929": "ZW1215"}
cosmo = FlatLambdaCDM(71, 0.27)   # ReadMe do ACCEPT2.0: cosmologia dos perfis
ref = FlatLambdaCDM(70, 0.30)     # convencao do X-COP usada na calibracao anterior

def read(path, ncol):
    rows = []
    for line in open(path, errors="replace"):
        p = line.split()
        try:
            v = [float(x) for x in p]
        except ValueError:
            continue
        if len(v) >= ncol:
            rows.append(v)
    return np.array(rows)

def profile(name):
    base = os.path.join(D, name, "analysis", "spectral")
    f = os.path.join(base, "deproj_profiles.dat")
    if not os.path.exists(f):
        f = os.path.join(base, "deproj_profiles_interp.dat")
    d = read(f, 8)
    rin, rout, ne, elo, ehi = d[:, 0], d[:, 1], d[:, 5], d[:, 6], d[:, 7]
    # escala (kpc/arcsec) implícita no volume de casca tabulado pelo próprio arquivo
    s = np.nan
    nf = os.path.join(base, "kT_Z_norm_deproj.dat")
    if os.path.exists(nf):
        k = read(nf, 4)
        vol, a, b = k[:, -1], k[:, 0], k[:, 1]
        s = np.median((vol / (4 / 3 * np.pi * (b**3 - a**3))) ** (1 / 3)) / KPC
    return rin, rout, ne, elo, ehi, s

res, tabs = [], []
for name in list(CLASH) + list(XCOP):
    base = os.path.join(D, name, "analysis", "spectral")
    main = os.path.exists(os.path.join(base, "deproj_profiles.dat"))
    if not main and not os.path.exists(os.path.join(base, "deproj_profiles_interp.dat")):
        res.append(dict(cluster=name, status="sem arquivo")); continue
    rin, rout, ne, elo, ehi, s = profile(name)
    z = CLASH.get(name)
    if z is None:
        z = {"A85": .0558, "A644": .0698, "A1795": .0625, "A2029": .0765, "A2142": .0898, "A2319": .0562, "ZW1215": .075}[XCOP[name]]
    s70 = cosmo.angular_diameter_distance(z).to_value(u.kpc) * np.pi / 648000
    f70 = float(ref.angular_diameter_distance(z) / cosmo.angular_diameter_distance(z))
    scale = s70
    a, b = rin * scale, rout * scale
    dm = 4 / 3 * np.pi * ((b * KPC) ** 3 - (a * KPC) ** 3) * ne * MUE * MP / MSUN
    M = np.cumsum(dm)
    # conversao 71/0.27 -> 70/0.30: r ~ D_A, n_e ~ D_A^-1/2, M ~ D_A^5/2
    a, b, ne, M = a * f70, b * f70, ne * f70**-0.5, M * f70**2.5
    elo, ehi = elo * f70**-0.5, ehi * f70**-0.5
    gaps = int(np.sum(~np.isclose(rin[1:], rout[:-1])))
    r = dict(cluster=name, grupo="CLASH" if name in CLASH else "X-COP", status="ok" if main else "so_interp", z=z, n_bins=len(ne),
             r_first_kpc=round(b[0], 1), r_max_kpc=round(b[-1]), kpc_arcsec_arquivo=round(s, 3),
             kpc_arcsec_H71=round(s70, 3), fator_D70_D71=round(f70, 4), lacunas=gaps, ne_nao_positivo=int(np.sum(ne <= 0)),
             err_med_ne=round(float(np.median((elo + ehi) / 2 / np.where(ne > 0, ne, np.nan))), 3))
    if name in XCOP:
        h = fits.open(glob.glob(os.path.join(XC, XCOP[name], "*fgas_profile.fits"))[0])
        t, R500 = h[1].data, None
        for hd in h:
            if "R500" in hd.header: R500 = hd.header["R500"]
        xr, xm = np.asarray(t["RADIUS"], float) * R500, np.asarray(t["MGAS"], float)
        for rr in (50, 100, 200, 300, 500):
            if xr.min() <= rr <= min(b[-1], xr.max()) and b[0] <= rr:
                m_a = np.exp(np.interp(np.log(rr), np.log(b), np.log(M)))
                m_x = np.exp(np.interp(np.log(rr), np.log(xr), np.log(xm)))
                r[f"dex_{rr}"] = round(float(np.log10(m_a / m_x)), 3)
    res.append(r)
    tabs.append(pd.DataFrame(dict(system_id=name, r_in_kpc=a, r_kpc=b, ne_cm3=ne, ne_err_lo=elo, ne_err_hi=ehi, Mgas_msun=M)))
os.makedirs(OUT, exist_ok=True)
R = pd.DataFrame(res); R.to_csv(os.path.join(OUT, "cobertura.csv"), index=False)
pd.concat(tabs).to_csv(os.path.join(OUT, "mgas_accept2_H70.csv"), index=False)
pd.set_option("display.width", 250)
print(R.to_string(index=False))
