import numpy as np, json, csv, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import ce_fit as cf, ce_osc as osc, ellipsometry_fit as ef, tmm_fix

class fin:
    @staticmethod
    def NfB(p, wl, ng):
        E = 1239.841984/wl
        e1, e2 = osc.drude_rt(E, p[2], p[3])
        a, b = osc.tauc_lorentz(E, p[4], p[5], p[6], p[7]); e1 += a; e2 += b
        if ng:
            a, b = osc.gaussian(E, p[8], p[9], p[10]); e1 += a; e2 += b
        N = np.sqrt((p[1] + e1 + 1j*e2).astype(complex))
        return np.where(N.imag < 0, -N, N)

R = json.load(open('ce_v13_result.json'))
V12 = json.load(open('ce_v12_result.json'))
LAB = {'#1':'VendorA ITO-1','#2':'VendorA ITO-2','#3':'VendorB IZO-1','#4':'VendorB IZO-2','#5':'VendorA ITO 2%O2'}
OLD = {'#1':(50.0,2.093,0.0122),'#2':(50.5,2.084,0.0372),'#3':(42.0,2.05,0.0),'#4':(42.0,2.05,0.0),'#5':(42.0,2.11,0.115)}

fig = plt.figure(figsize=(15,8.4)); gs = fig.add_gridspec(2,3,hspace=.32,wspace=.26)
wp = np.linspace(340,1080,600)

# fit quality for the two extremes
for col,sh in enumerate(['#1','#5']):
    p=np.array(R[sh]['p']); d,dth=R[sh]['d'],R[sh]['dth']
    wl,Pm,Dm=cf.load(sh); m=(wl>=340)&(wl<=1080); wl,Pm,Dm=wl[m],Pm[m],Dm[m]
    ox,si=cf._mat('NTVE_JAW',wl),cf._mat('SI_JAW',wl)
    Nf=fin.NfB(p,wl,R[sh]['ng']); Nr=ef.bruggeman_ema50(Nf,np.ones_like(Nf)); amb=np.ones_like(Nf)
    ax=fig.add_subplot(gs[0,col]); ax2=ax.twinx()
    for a,an in enumerate(cf.ANG0+dth):
        rp,rs=tmm_fix.tmm(wl,[amb,Nr,Nf,ox,si],[p[0],d,cf.D_OX],an)
        rr=np.conj(rp/rs); ps=np.rad2deg(np.arctan(np.abs(rr))); dl=np.rad2deg(np.angle(rr))%360
        ax.plot(wl,Pm[:,a],color='tab:red',lw=1.3,alpha=.75)
        ax.plot(wl,ps,'k--',lw=1.0)
        ax2.plot(wl,Dm[:,a]%360,color='tab:green',lw=1.3,alpha=.75)
        ax2.plot(wl,dl,'k:',lw=1.0)
    ax.set_xlim(340,1080); ax.set_ylim(0,50); ax2.set_ylim(0,360)
    ax.set_xlabel('wavelength (nm)'); ax.set_ylabel(r'$\Psi$ (deg)',color='tab:red')
    ax2.set_ylabel(r'$\Delta$ (deg)',color='tab:green')
    ax.set_title('%s   MSE = %.2f   (was %.1f with the sign bug)'%(LAB[sh],R[sh]['mse'],
                 {'#1':21.9,'#5':51.7}[sh]),fontsize=10)

ax=fig.add_subplot(gs[0,2]); ax.axis('off')
ax.text(0,1,'what the transfer-matrix fix changed',fontsize=11,weight='bold',va='top')
ax.text(0,.90,
 'The propagation matrix had its two exponents swapped, so the\n'
 'backward wave was AMPLIFIED inside absorbing layers.  A 200 nm\n'
 'Ag film returned R = 1.033 (impossible) instead of the bulk\n'
 'Fresnel value 0.9685.  Ellipsometry only uses rp/rs, so the\n'
 'error hid completely - n and k simply absorbed it.\n\n'
 '                     before        after\n'
 '  MSE  #1            21.9          2.32\n'
 '  MSE  #5            51.7          1.80\n'
 '  thickness #1       50.0 nm       55.5 nm\n'
 '  Eo2 / Eg2 / Br2    +-64/58/22%   +-0-3%\n'
 '  Drude rho, tau     +-596%        +-3%\n\n'
 'The large error bars on Eo2/Eg2/Br2 were the BUG, not a limit\n'
 'of the data.  With it fixed every oscillator parameter is\n'
 'determined to a few percent.',
 fontsize=8.5,va='top',family='monospace')

axn=fig.add_subplot(gs[1,0]); axk=fig.add_subplot(gs[1,1])
rows=[]
for sh in ['#1','#2','#3','#4','#5']:
    N=fin.NfB(np.array(R[sh]['p']),wp,R[sh]['ng'])
    axn.plot(wp,N.real,lw=1.7,label=LAB[sh]); axk.plot(wp,N.imag,lw=1.7,label=LAB[sh])
    rows.append((sh,N))
for a,t in ((axn,'n'),(axk,'k')):
    a.set_xlabel('wavelength (nm)'); a.set_ylabel(t); a.grid(alpha=.3)
axn.legend(fontsize=8); axk.set_ylim(bottom=0)
axn.set_title('n  (corrected TMM, 340-1080 nm)',fontsize=10)
axk.set_title('k',fontsize=10)

ax=fig.add_subplot(gs[1,2]); ax.axis('off')
tab='sample            d(nm)  dth    n550   k550    MSE\n'+'-'*50+'\n'
for sh,N in rows:
    i=np.argmin(abs(wp-550))
    tab+='%-16s %5.1f %+5.2f  %5.3f  %6.4f  %5.2f\n'%(LAB[sh],R[sh]['d'],R[sh]['dth'],N.real[i],N.imag[i],R[sh]['mse'])
ax.text(0,1,tab,fontsize=9,va='top',family='monospace')
ax.text(0,.60,
 'VendorA ITO-1 and -2 are nominally the same film and land on\n'
 'd = 55.5 / 56.0 nm with MSE 2.32 / 2.33 - an independent\n'
 'consistency check that the corrected geometry is right.\n\n'
 'The 2% O2 sample is the transparent one (k550 = 0.006 vs\n'
 '0.06-0.20): oxygen fills vacancies, carriers drop, free-\n'
 'carrier absorption disappears.\n\n'
 'VendorB IZO-2 (#4) still fits poorly (MSE 16, chain residual\n'
 '0.012 vs 0.0009 for the others) - that sample needs a look.',
 fontsize=8.5,va='top')
fig.suptitle('ITO / IZO 260813 - refit with the corrected transfer matrix, 340-1080 nm',fontsize=12)
fig.savefig('ce_v13_check.png',dpi=125,bbox_inches='tight')

with open('ITO_IZO_v13_nk.csv','w',newline='') as f:
    w=csv.writer(f)
    w.writerow(['# ITO/IZO 260813 n,k - corrected transfer matrix, valid range 340-1080 nm'])
    w.writerow(['# geometry: '+'; '.join('%s d=%.1fnm dth=%+.2fdeg rough=%.2fnm MSE=%.2f'
                %(LAB[s],R[s]['d'],R[s]['dth'],R[s]['p'][0],R[s]['mse']) for s in R)])
    hdr=['wl_nm']
    for s in R: hdr += ['n_'+LAB[s].replace(' ','_'), 'k_'+LAB[s].replace(' ','_')]
    w.writerow(hdr)
    Ns={s:fin.NfB(np.array(R[s]['p']),wp,R[s]['ng']) for s in R}
    for i,x in enumerate(wp):
        row=[round(x,2)]
        for s in R: row += [round(Ns[s].real[i],4), round(Ns[s].imag[i],5)]
        w.writerow(row)
print('saved ce_v13_check.png and ITO_IZO_v13_nk.csv')
for sh,N in rows:
    i=np.argmin(abs(wp-550)); j=np.argmin(abs(wp-633))
    print('  %-16s d=%.1f  n550=%.3f k550=%.4f  n633=%.3f k633=%.4f  (was d=%.1f n550=%.3f k550=%.4f)'
          %(LAB[sh],R[sh]['d'],N.real[i],N.imag[i],N.real[j],N.imag[j],*OLD[sh]))
