"""Reduz o catálogo UniverseMachine. Uso:
   python extrair_umachine.py loader.py sfr_catalog_0.760300.bin smdpl_z032.npz
O loader.py é apenas LIDO (para obter o formato do registro); nada nele é executado."""
import sys, ast, os
import numpy as np

def read_dtype(path):
    for node in ast.walk(ast.parse(open(path).read())):
        if isinstance(node, ast.Call) and getattr(node.func, 'attr', '') == 'dtype':
            kw = {k.arg: k.value for k in node.keywords}
            spec = node.args[0] if node.args else kw['dtype']
            align = ast.literal_eval(kw['align']) if 'align' in kw else False
            return np.dtype(ast.literal_eval(spec), align=align)
    sys.exit('Não achei np.dtype no loader.py')

def main(loader, cat, out):
    dt = read_dtype(loader)
    need = ['id', 'upid', 'pos', 'm', 'sm', 'obs_sm']
    miss = [f for f in need if f not in dt.names]
    if miss: sys.exit(f'Campos ausentes: {miss}\nCampos do arquivo: {dt.names}')
    size = os.path.getsize(cat)
    if size % dt.itemsize: sys.exit(f'Tamanho {size} não é múltiplo do registro {dt.itemsize}: download incompleto?')
    n = size // dt.itemsize; keep = []; done = 0
    print(f'{n} registros de {dt.itemsize} bytes')
    with open(cat, 'rb') as f:
        while done < n:
            a = np.fromfile(f, dtype=dt, count=2_000_000); done += len(a)
            s = (a['sm'] > 10**7.9) | (a['obs_sm'] > 10**7.9) | ((a['upid'] == -1) & (a['m'] > 1e11))
            a = a[s]
            host = np.where(a['upid'] == -1, a['id'], a['upid'])
            keep.append((a['pos'][:, :3].astype('f4'), a['m'].astype('f4'), a['sm'].astype('f4'),
                         a['obs_sm'].astype('f4'), a['upid'] == -1, host))
            print(f'\r{100*done/n:5.1f}%', end='', flush=True)
    cols = [np.concatenate(c) for c in zip(*keep)]
    np.savez_compressed(out, pos=cols[0], m=cols[1], sm=cols[2], obs_sm=cols[3], central=cols[4], host=cols[5])
    print(f'\nSalvo {out}: {len(cols[1])} objetos, {os.path.getsize(out)/1e6:.0f} MB')

if __name__ == '__main__':
    main(*sys.argv[1:4])
