import re, base64, gzip, os, numpy as np
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

D=str(DATA_DIR)

def jaw(fn, layer):
    t=open(os.path.join(D,fn),encoding='latin-1').read()
    m=re.search(r'start_Layer%d\s*\n(.*?)end_Layer%d'%(layer,layer),t,re.S)
    o={}
    for tag in ['Wvl','e1','e2']:
        mm=re.search(r"start_%s Array\s*\n\s*'([^']+)'"%tag,m.group(1))
        o[tag]=np.frombuffer(gzip.decompress(base64.b64decode(mm.group(1))),dtype='>f4').astype(float)
    return o

si=jaw('NO_on_Si.mod',0); ox=jaw('NO_on_Si.mod',1)
si_wl=si['Wvl']/10.0; ox_wl=ox['Wvl']/10.0
nSi=np.sqrt(complex(np.interp(633,si_wl,si['e1']),np.interp(633,si_wl,si['e2'])))
print('SI wl %.1f-%.1f nm  N@633=%.3f+%.4fj'%(si_wl[0],si_wl[-1],nSi.real,nSi.imag))
print('OX wl %.1f-%.1f nm  n@633=%.4f'%(ox_wl[0],ox_wl[-1],np.sqrt(np.interp(633,ox_wl,ox['e1']))))

def ce_nodes(fn):
    t=open(os.path.join(D,fn),encoding='latin-1').read()
    m=re.search(r'start_Layer2\s*\n(.*?)end_Layer2',t,re.S)
    body=m.group(1)
    nodes=re.findall(r"([\d.eE+-]+)\s+[TF]\s+[-\d.eE]+\s+[\d.eE+-]+\s+[TF]\s+'spline_e2\(([\d.]+)\)'",body)
    E=np.array([float(b) for a,b in nodes]); V=np.array([float(a) for a,b in nodes])
    einf=float(re.search(r"([\d.eE+-]+)\s+T\s+-10\.0\s+100\.0\s+F\s+'E Inf'",body).group(1))
    iramp=float(re.search(r"([\d.eE+-]+)\s+T\s+0\.0\s+1000\.0\s+F\s+'IR Amp'",body).group(1))
    o=np.argsort(E)
    return E[o],V[o],einf,iramp

E1,V1,ei1,ir1=ce_nodes('after_200C_N2_1h_IZO1_NO_on_Si.mod')
E3,V3,ei3,ir3=ce_nodes('IZO3_on_NO.mod')
print('N2_IZO1 CE: %d nodes E %.3f-%.3f  einf=%.3f irAmp=%.3f'%(len(E1),E1[0],E1[-1],ei1,ir1))
print('  e2 @0.735=%.3f @1.5=%.3f @3.0=%.3f @4.5=%.3f'%(
      V1[np.argmin(abs(E1-0.735))],V1[np.argmin(abs(E1-1.5))],V1[np.argmin(abs(E1-3.0))],V1[np.argmin(abs(E1-4.5))]))
print('IZO3pre CE: %d nodes einf=%.3f irAmp=%.3f  e2@0.735=%.3f'%(len(E3),ei3,ir3,V3[np.argmin(abs(E3-0.735))]))

d=dict(np.load(str(OUT_DIR / 'izo_data.npz')))
d.update(si_wl=si_wl,si_e1=si['e1'],si_e2=si['e2'],ox_wl=ox_wl,ox_e1=ox['e1'],ox_e2=ox['e2'],
         ceE1=E1,ceV1=V1,ce_einf1=np.array(ei1),ce_ir1=np.array(ir1),
         ceE3=E3,ceV3=V3,ce_einf3=np.array(ei3),ce_ir3=np.array(ir3))
np.savez(str(OUT_DIR / 'izo_data.npz'),**d)
print('izo_data.npz updated:',sorted(d.keys()))
