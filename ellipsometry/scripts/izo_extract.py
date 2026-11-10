import re, base64, gzip, numpy as np, openpyxl
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

D=str(DATA_DIR)

# 1) JAW substrate tables from NO_on_Si.mod (float32)
t=open(D+r'\NO_on_Si.mod',encoding='latin-1').read()
def arrs(layer):
    m=re.search(r'start_Layer%d\s*\n(.*?)end_Layer%d'%(layer,layer),t,re.S)
    out={}
    for tag in ['Wvl','e1','e2']:
        mm=re.search(r"start_%s Array\s*\n\s*'([^']+)'"%tag,m.group(1))
        out[tag]=np.frombuffer(gzip.decompress(base64.b64decode(mm.group(1))),dtype='<f4').astype(float)
    return out
si=arrs(0); ox=arrs(1)
print('SI_JAW wl %g-%g pts%d  e1@633=%.3f e2@633=%.4f'%(si['Wvl'][0],si['Wvl'][-1],len(si['Wvl']),
      np.interp(633,si['Wvl'],si['e1']),np.interp(633,si['Wvl'],si['e2'])))
print('NTVE  wl %g-%g pts%d  e1@633=%.4f'%(ox['Wvl'][0],ox['Wvl'][-1],len(ox['Wvl']),np.interp(633,ox['Wvl'],ox['e1'])))

# 2) N2_IZO1 xlsx: 4 angles interleaved, depol cols 42-46
wb=openpyxl.load_workbook(D+r'\N2_1h_IZO1.xlsx',data_only=True,read_only=True)
ws=wb['Sheet1']
rows=[r for r in ws.iter_rows(min_row=4,values_only=True)]
wb.close()
rows=[r for r in rows if isinstance(r[1],(int,float))]
n2wl=np.array([r[1] for r in rows])
n2psi=np.array([[r[2+2*j] for j in range(4)] for r in rows])
n2del=np.array([[r[3+2*j] for j in range(4)] for r in rows])
n2dep=np.array([[r[43+j] for j in range(4)] for r in rows])
print('\nN2_IZO1: %d pts %.1f-%.1f nm, angles 45/55/65/75'%(len(n2wl),n2wl[0],n2wl[-1]))
for lo,hi in [(245,600),(600,1100),(1100,1700)]:
    m=(n2wl>=lo)&(n2wl<hi)
    print('  depol %4d-%4d: mean|d|=%.2f%% max=%.2f%%'%(lo,hi,np.abs(n2dep[m]).mean(),np.abs(n2dep[m]).max()))

# 3) O2_IZO3 xlsx: 5 angles in 4-col groups; model curves in cols 21-37; depol 162..187
wb=openpyxl.load_workbook(D+r'\O2_1h_IZO3.xlsx',data_only=True,read_only=True)
ws=wb['Sheet1']
rows=[r for r in ws.iter_rows(min_row=5,values_only=True)]
wb.close()
rows=[r for r in rows if isinstance(r[1],(int,float))]
o2wl=np.array([r[1] for r in rows])
gcols=[1,5,9,13,17]      # data blocks: [wl,psi,del]
mcols=[21,25,29,33,37]   # CompleteEASE model blocks
o2psi=np.array([[r[c+1] for c in gcols] for r in rows],dtype=float)
o2del=np.array([[r[c+2] for c in gcols] for r in rows],dtype=float)
o2psi_m=np.array([[r[c+1] for c in mcols] for r in rows],dtype=float)
o2del_m=np.array([[r[c+2] for c in mcols] for r in rows],dtype=float)
depc=[163,169,175,181,187]
o2dep=np.array([[r[c] for c in depc] for r in rows],dtype=float)
print('\nO2_IZO3: %d pts %.1f-%.1f nm, angles 45/50/55/60/65'%(len(o2wl),o2wl[0],o2wl[-1]))
for lo,hi in [(245,600),(600,1100),(1100,1700)]:
    m=(o2wl>=lo)&(o2wl<hi)
    print('  depol %4d-%4d: mean|d|=%.2f%% max=%.2f%%'%(lo,hi,np.abs(o2dep[m]).mean(),np.abs(o2dep[m]).max()))
# CompleteEASE model quality vs data
dd=(o2del_m-o2del+180)%360-180
print('  CE-model residual: RMSE_psi=%.2f RMSE_del=%.2f (their failed fit)'%(
      np.sqrt(np.nanmean((o2psi_m-o2psi)**2)),np.sqrt(np.nanmean(dd**2))))

np.savez(str(OUT_DIR / 'izo_data.npz'),
         si_wl=si['Wvl'],si_e1=si['e1'],si_e2=si['e2'],
         ox_wl=ox['Wvl'],ox_e1=ox['e1'],ox_e2=ox['e2'],
         n2wl=n2wl,n2psi=n2psi,n2del=n2del,n2dep=n2dep,
         o2wl=o2wl,o2psi=o2psi,o2del=o2del,o2dep=o2dep,
         o2psi_ce=o2psi_m,o2del_ce=o2del_m)
print('\nsaved izo_data.npz')
