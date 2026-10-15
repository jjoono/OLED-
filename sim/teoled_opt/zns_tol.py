import numpy as np, zns_hi2 as Z
from multiprocessing import Pool
def f(a): return a, Z.eqe(*a)
if __name__ == '__main__':
    J = []
    for b in ((130, 254, 100), (100, 252, 100)):
        J.append(b)
        for i in range(3):
            for s in (-10, 10):
                x = list(b); x[i] += s; J.append(tuple(x))
    J += [(130, 78, 0), (130, 78, 100), (130, 254, 0)]
    with Pool(4) as p:
        for a, e in p.map(f, J): print('dE %5.0f dH %5.0f ZnS %5.0f  EQE %.3f' % (*a, e))
