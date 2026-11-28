import fsspec, h5py, re, sys, time
url = sys.argv[1] if len(sys.argv) > 1 else 'https://portal.nersc.gov/project/lsst/cosmoDC2/z_0_1.step_all.healpix_8786.hdf5'
pat = re.compile(r'^(ra|dec|redshift|stellar|shear|conv|halo|is_?central|mag_(true_)?r|upid|host)', re.I)
t = time.time()
with fsspec.open(url, 'rb', block_size=2**20) as f:
    h = h5py.File(f, 'r')
    keys = list(h.keys()); print('grupos:', keys[:10], '...' if len(keys) > 10 else '')
    g = h[[k for k in keys if k.lower() != 'metadata'][0]]
    while len(g.keys()) < 30 and any(isinstance(g[k], h5py.Group) for k in g.keys()):
        g = g[[k for k in g.keys() if isinstance(g[k], h5py.Group)][0]]
    cols = list(g.keys()); print('grupo lido:', g.name, '| colunas:', len(cols))
    for c in cols:
        if pat.match(c) and isinstance(g[c], h5py.Dataset):
            print(' ', c, g[c].shape, g[c].dtype)
    d = next(g[c] for c in cols if c.lower() == 'redshift' or c.lower().startswith('redshift'))
    print('leitura de teste:', d[:3])
print(f'tempo: {time.time()-t:.0f} s')
