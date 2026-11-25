import numpy as np, json, csv, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import ag_load as L, ag_final as AFIN
class AF:
    agN = staticmethod(lambda p, wl: AFIN.agN(np.array(p), wl))

R = json.load(open('ag_final_result.json'))
HAT = [s for s in L.ORDER if L.SPEC[s][0]=='HATCN' and L.SPEC[s][2]>0]
MOO = [s for s in L.ORDER if L.SPEC[s][0]=='MoOx'  and L.SPEC[s][2]>0]
wp = np.linspace(260, 1080, 600)

fig = plt.figure(figsize=(15.5, 9.2)); gs = fig.add_gridspec(2,3,hspace=.30,wspace=.27)
for col,(grp,title) in enumerate([(HAT,'Ag on HATCN'),(MOO,'Ag on MoOx')]):
    axn = fig.add_subplot(gs[0,col]); axk = fig.add_subplot(gs[1,col])
    cm = plt.cm.viridis(np.linspace(0,.9,len(grp)))
    for c,s in zip(cm,grp):
        N = AF.agN(np.array(R[s]['p']), wp)
        lab='%d nm (fit %.1f)'%(L.SPEC[s][2], R[s]['p'][0])
        axn.plot(wp,N.real,color=c,lw=1.6,label=lab); axk.plot(wp,N.imag,color=c,lw=1.6,label=lab)
    for a,t in ((axn,'n'),(axk,'k')):
        a.set_xlabel('wavelength (nm)'); a.set_ylabel(t); a.grid(alpha=.3)
    axn.set_title('%s   n'%title,fontsize=10); axk.set_title('%s   k'%title,fontsize=10)
    axn.legend(fontsize=7.5,title='nominal Ag',title_fontsize=7.5); axn.set_ylim(0,3.2); axk.set_ylim(0,6)

ax = fig.add_subplot(gs[0,2])
for grp,mk,lab in [(HAT,'o','HATCN'),(MOO,'s','MoOx')]:
    d=[R[s]['p'][0] for s in grp]; y=[R[s]['n'][2] for s in grp]
    ax.plot(d,y,mk+'-',lw=1.6,ms=6,label=lab)
ax.axhline(0.135,ls=':',c='k',lw=1); ax.text(11.5,0.19,'bulk Ag n633 ~ 0.14',fontsize=7.5,ha='right')
ax.set_xlabel('fitted Ag thickness (nm)'); ax.set_ylabel('n at 633 nm'); ax.grid(alpha=.3); ax.legend(fontsize=8)
ax.set_title('metallisation: n(633) falls to the bulk value\nonly once the film is continuous',fontsize=10)

HB=6.582119569e-16; EPS0=8.8541878128e-12
def wp_of(p): return np.sqrt(HB**2/(EPS0*(p[3]/100.0)*(p[4]*1e-15)))
WP_REF = max(wp_of(np.array(R[s]['p'])) for s in HAT)
ax2 = fig.add_subplot(gs[1,2])
for grp,mk,lab in [(HAT,'o','HATCN'),(MOO,'s','MoOx')]:
    d=[R[s]['p'][0] for s in grp]
    y=[100*(wp_of(np.array(R[s]['p']))/WP_REF)**2 for s in grp]
    ax2.plot(d,y,mk+'-',lw=1.7,ms=7,label=lab)
    for xx,yy,s in zip(d,y,grp):
        ax2.annotate(str(L.SPEC[s][2]),(xx,yy),textcoords='offset points',xytext=(0,-14),fontsize=7,ha='center')
ax2.axhline(100,ls=':',c='k',lw=1)
ax2.set_xlabel('fitted Ag thickness (nm)   [label = nominal nm]')
ax2.set_ylabel('effective Ag density  (wp/wp_ref)^2  [%]')
ax2.grid(alpha=.3); ax2.legend(fontsize=8); ax2.set_ylim(30,115)
ax2.set_title('free-electron density from the plasma energy' + chr(10) +
              '100%% = 12 nm on HATCN (hw_p = %.2f eV; bulk Ag ~9.2 eV)'%WP_REF, fontsize=10)
fig.suptitle('260819 thin Ag on HATCN vs MoOx - Gen-Osc n,k (Einf + Drude + 3 Gaussians), 260-1080 nm',fontsize=12)
fig.savefig('Ag260819_v2_nk.png',dpi=125,bbox_inches='tight')

with open('Ag260819_v2_nk.csv','w',newline='') as f:
    w=csv.writer(f)
    w.writerow(['# 260819 thin Ag n,k - seed layer held fixed, native oxide 2.0 nm, corrected TMM'])
    w.writerow(['# '+ '; '.join('%s %s/Ag nom %d -> d_Ag %.2f nm, rough %.2f nm, MSE %.2f'
                %(s,R[s]['seed'],R[s]['ag_nom'],R[s]['p'][0],R[s]['p'][1],R[s]['mse']) for s in L.ORDER if s in R)])
    hdr=['wl_nm']
    for s in L.ORDER:
        if s in R: hdr += ['n_%s_%s%dnm'%(s,R[s]['seed'],R[s]['ag_nom']), 'k_%s_%s%dnm'%(s,R[s]['seed'],R[s]['ag_nom'])]
    w.writerow(hdr)
    Ns={s:AF.agN(np.array(R[s]['p']),wp) for s in R}
    for i,x in enumerate(wp):
        row=[round(x,2)]
        for s in L.ORDER:
            if s in R: row += [round(Ns[s].real[i],4), round(Ns[s].imag[i],4)]
        w.writerow(row)
print('saved Ag260819_v2_nk.png and Ag260819_v2_nk.csv')
