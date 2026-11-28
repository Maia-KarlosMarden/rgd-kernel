import sys, os, time
import numpy as np, h5py, fsspec
PIX = [9685,9558,9814,9686,9557,9813,9684,9429,9941,9559,9430,9815,9942,9428,9556,9687,9940,9812,9302,9301,
       10070,10069,9431,9683,9303,9943,9560,9816,9173,10071]
BASE = os.environ.get('DC2_BASE', 'https://portal.nersc.gov/project/lsst/cosmoDC2/')
MAG = 'SDSS_filters/magnitude:SDSS_r:observed:dustAtlas'
n = int(sys.argv[1]) if len(sys.argv) > 1 else 1
for p in PIX[:n]:
    out = 'dc2_%d.npz' % p
    if os.path.exists(out):
        print(p, 'ja feito'); continue
    t = time.time()
    for tentativa in range(3):
        try:
            with fsspec.open(BASE + 'z_0_1.step_all.healpix_%d.hdf5' % p, 'rb', block_size=2**22) as f:
                g = h5py.File(f, 'r')['galaxyProperties']
                z = g['redshiftHubble'][:]; zr = g['redshift'][:]; mag = g[MAG][:]
                ms = g['totalMassStellar'][:]; cen = g['isCentral'][:]; mh = g['hostHaloMass'][:]
                s = (z > 0.03) & (z < 0.65) & ((mag < 22.5) | (cen & (mh > 1e11)))
                d = dict(ra=g['ra_true'][:][s], dec=g['dec_true'][:][s], z=z[s], zr=zr[s], mag=mag[s],
                         ms=ms[s], cen=cen[s], mh=mh[s], tag=g['hostHaloTag'][:][s])
                rng = np.random.default_rng(p)
                src = np.where((zr > 0.4) & (zr < 1.0))[0]
                src = np.sort(rng.choice(src, len(src)//20, replace=False))
                d.update(s_z=zr[src])
                for nome, col in (('s_ra','ra'), ('s_dec','dec'), ('s_g1','shear1'), ('s_g2','shear2')):
                    d[nome] = g[col][:][src]
            break
        except Exception as e:
            print(p, 'erro, tentando de novo:', type(e).__name__, str(e)[:200]); time.sleep(10)
    else:
        sys.exit('falhou 3 vezes')
    np.savez_compressed('tmp_dc2.npz', **d); os.replace('tmp_dc2.npz', out)
    big = d['cen'] & (d['mh'] > 1e11)
    print('%d: total %d | guardadas %d | fontes %d | log M* mediano (r<20) %.2f | log Mh mediano %.2f | %.0f MB | %.0f s' % (
        p, len(z), s.sum(), len(src), np.median(np.log10(d['ms'][d['mag'] < 20])),
        np.median(np.log10(d['mh'][big])), os.path.getsize(out)/1e6, time.time()-t), flush=True)
