"""ZnS 100-250 nm, ETL and HTL kept thin (first cavity order): grid d_ETL 50-130, d_HTL 30-150 (10 nm steps),
27 wavelengths. Ag 8 nominal (joint nk), Ag/Al0.2/Liq, HATCN 10, Al2O3 30 / n1.75 3 um / Al2O3 30 / n1.75 MLA."""
import csv, itertools, numpy as np, zns_hi2 as Z
from multiprocessing import Pool
ZS = (100, 125, 150, 175, 200, 225, 250); DE = range(50, 131, 10); DH = range(30, 151, 10)
def f(a): z, dE, dH = a; return (z, dE, dH, Z.eqe(dE, dH, z))
if __name__ == '__main__':
    with Pool(4) as p: res = p.map(f, list(itertools.product(ZS, DE, DH)), chunksize=6)
    with open('zns_thin_grid.csv', 'w', newline='') as fh: csv.writer(fh).writerows([['ZnS', 'd_ETL', 'd_HTL', 'EQE']] + res)
    print('done')
