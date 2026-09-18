"""Export candidate colormaps as (a) Origin .pal files and (b) RGB stop tables
for Origin's 'introduce other colors in mixing' dialog."""
import numpy as np, struct, os, matplotlib.pyplot as plt
from cmcrameri import cm as cmc

MAPS = {'inferno': plt.get_cmap('inferno'), 'magma': plt.get_cmap('magma'),
        'lipari': cmc.lipari, 'batlow': cmc.batlow, 'viridis': plt.get_cmap('viridis')}
OUT = 'origin_palettes'

def rgb255(cmap, x):
    return np.round(np.asarray(cmap(x))[..., :3]*255).astype(int)

def write_pal(path, cmap, n=256):                      # Microsoft RIFF PAL, what Origin imports
    cols = rgb255(cmap, np.linspace(0, 1, n))
    data = b''.join(struct.pack('BBBB', *c, 0) for c in cols)
    chunk = struct.pack('<HH', 0x0300, n) + data
    body = b'PAL ' + b'data' + struct.pack('<I', len(chunk)) + chunk
    open(path, 'wb').write(b'RIFF' + struct.pack('<I', len(body)) + body)

def stop_error(cmap, k):                               # max RGB error of a k-stop linear reproduction
    xs = np.linspace(0, 1, k); stops = np.asarray(cmap(xs))[:, :3]
    xf = np.linspace(0, 1, 256); true = np.asarray(cmap(xf))[:, :3]
    approx = np.stack([np.interp(xf, xs, stops[:, c]) for c in range(3)], 1)
    return np.abs(approx-true).max()*255

report = []
for name, cmap in MAPS.items():
    write_pal(os.path.join(OUT, f'{name}.pal'), cmap)
    errs = {k: stop_error(cmap, k) for k in (5, 7, 9, 11, 13)}
    k = next(k for k in (5, 7, 9, 11, 13) if errs[k] <= 8) if min(errs.values()) <= 8 else 13
    report.append((name, k, errs))
    xs = np.linspace(0, 1, k); cols = rgb255(cmap, xs)
    lines = [f'# {name}: {k} stops for Origin colour mixing (From = stop 1, To = stop {k})',
             '# level  R    G    B    hex']
    for i, (x, c) in enumerate(zip(xs, cols), 1):
        lines.append(f'{i:>4}  {x:4.2f}  {c[0]:>3}  {c[1]:>3}  {c[2]:>3}   #{c[0]:02X}{c[1]:02X}{c[2]:02X}')
    open(os.path.join(OUT, f'{name}_stops.txt'), 'w').write('\n'.join(lines)+'\n')

for name, k, errs in report:
    print(f'{name:8s} -> {k:2d} stops  (max RGB error {errs[k]:.1f}/255; ' +
          ', '.join(f'{kk}:{v:.0f}' for kk, v in errs.items()) + ')')

# verify the .pal files read back
for name in MAPS:
    b = open(os.path.join(OUT, f'{name}.pal'), 'rb').read()
    assert b[:4] == b'RIFF' and b[8:12] == b'PAL ' and b[12:16] == b'data'
    n = struct.unpack('<H', b[22:24])[0]
    print(f'  {name}.pal OK: {n} colours, {len(b)} bytes, first #{b[24]:02X}{b[25]:02X}{b[26]:02X} last #{b[-4]:02X}{b[-3]:02X}{b[-2]:02X}')
