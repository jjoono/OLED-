import itertools, csv, numpy as np, tmodel as T
from multiprocessing import Pool
b = np.load('mla_hemi_1.65.npz'); BT, BR = b['BT'], b['BR']
DE = range(20, 301, 20); DH = range(20, 201, 30); HAT = (5, 10, 20); AGD = (8, 12, 15, 20); AGT = ('meas', 'bulk')
def f(a):
    agt, dAg, dHAT, dE, dH = a
    return (agt, dAg, dHAT, dE, dH, T.eqe(dE, dH, dHAT, dAg, BT, BR, agtop=agt))
if __name__ == '__main__':
    jobs = list(itertools.product(AGT, AGD, HAT, DE, DH))
    with Pool(4) as p: res = p.map(f, jobs, chunksize=20)
    with open('grid.csv', 'w', newline='') as fh: csv.writer(fh).writerows([['agtop', 'dAg', 'dHAT', 'dE', 'dH', 'EQE']] + res)
    print('done', len(res))
