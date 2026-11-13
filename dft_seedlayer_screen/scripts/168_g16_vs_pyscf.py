"""Recompute Gaussian 16 single points (workstation .out files) with PySCF at the same level.

    python scripts/168_g16_vs_pyscf.py ../benzene_t0_z0.out ...
"""
import re,sys,numpy as np
sys.path.insert(0,'scripts')
from pyscf import lib
lib.num_threads(4)
import qclocal as q
Z2S={1:'H',6:'C',7:'N',8:'O',9:'F',47:'Ag',15:'P',16:'S'}
def parse(fn):
    t=open(fn).read()
    blk=t.split('Input orientation:')[-1].split('---------------------------------------------------------------------')[2]
    S,X=[],[]
    for l in blk.strip().splitlines():
        p=l.split(); S.append(Z2S[int(p[1])]); X.append([float(v) for v in p[3:6]])
    m=re.search(r'Multiplicity = (\d+)',t); e=float(re.findall(r'SCF Done:\s+E\(\S+\)\s*=\s*(-?\d+\.\d+)',t)[-1])
    return S,np.array(X),int(m.group(1)),e
for fn in sys.argv[1:]:
    S,X,mult,eg=parse(fn)
    mf,e=q.scf(q.mol(S,X,mult))
    d3=q.d3bj(S,X)
    print(fn.split('/')[-1], len(S),'atoms', 'G16',round(eg,6),'PySCF+D3',round(e+d3,6),'diff meV',round((e+d3-eg)*27211.386,2), mf.converged, flush=True)
