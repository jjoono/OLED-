"""Real lens-free stack (ladder step D) with Al 1 nm instead of 0.2, and with ZnS >= 100 nm."""
import numpy as np, lensfree_ladder as LD
lc, wc = LD.L.spec(LD.c0, fine=False); lf, wf = LD.L.spec(LD.c0)
for al, zlo in ((0.2, 0.), (0.2, 100.), (1.0, 0.), (1.0, 100.)):
    f = lambda x: LD.real(x, lc, wc, LD.L.AGJ, al)
    e, x = LD.opt(f, [[61.5, 51.4, max(87.2, zlo + 5)], [62, 52, 220], [62, 250, 100]], [(15, 400), (15, 400), (zlo, 300)])
    print('Al %.1f ZnS>=%3.0f: eta_out %.3f (x0.96 = %.3f)  ETL/HTL/ZnS = %s' % (al, zlo, LD.real(x, lf, wf, LD.L.AGJ, al), 0.96 * LD.real(x, lf, wf, LD.L.AGJ, al), np.round(x, 1)), flush=True)
