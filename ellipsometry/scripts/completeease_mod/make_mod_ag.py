"""CompleteEASE .mod for the 260819 thin-Ag set.

The template has substrate + 2 layers (SI_JAW / NTVE_JAW / Gen-Osc).  The Ag
samples need substrate + 3 (NTVE / seed / Ag), so the whole Gen-Osc layer block
is duplicated and renumbered, and the layer count on the first Model Parms line
is bumped.  Layer numbering runs upward from the substrate, so Layer2 = seed and
Layer3 = Ag.
"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # measurement exports (.xlsx), CompleteEASE .mod
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fitted results, figures, intermediates

import re, os, json, numpy as np
import ag_load as L, ag_fit as AF

TEMPLATE = _os.path.join(ELLIPS_DATA, r'VendorA_ITO_GenOSC_6functions.mod')
OUTDIR   = _os.path.join(ELLIPS_OUT, r'mod_ag260819_v2')
TAIL = "\tF\tF\t0.0\t0.0\t100000.0\tF\t100.0\t"
num  = lambda x: repr(float(x))
pl   = lambda v,f,lo,hi,nm: "\t\t\t%s\t%s\t%s\t%s\tF\t'%s'%s"%(v,'T' if f else 'F',lo,hi,nm,TAIL)
gl   = lambda lo,hi,i: "\t\t\t0.0\tF\t%s\t%s\tF\t'Grade %d'%s"%(lo,hi,i,TAIL)

def blocks(einf, oscs, einf_free=True):
    """oscs = list of (type, common_eg, [(name,val,free,lo,hi,glo,ghi), ...])"""
    fit = ["\t\t\t%d\t%s\t%s\t0.5\t8.0\tF\t'Einf'%s" % (len(oscs), num(einf), 'T' if einf_free else 'F', TAIL)]
    grade = []
    for typ, ceg, lst in oscs:
        fit.append('\t\t\t' + typ + '\t')
        for nm, v, f, lo, hi, glo, ghi in lst:
            fit.append(pl(v, f, lo, hi, nm))
            grade.append('\t\t\tF\t'); grade += [gl(glo, ghi, j) for j in (1,2,3)]
        if ceg: fit.append('\t\t\tF\t')
    fit.append('\t\t\t'); grade.append('\t\t\t')
    return '\r\n'.join(fit), '\r\n'.join(grade)

def seed_oscs(p):          # p = [Einf, TL_Amp, TL_Br, TL_Eo, TL_Eg]
    return p[0], [("'Tauc-Lorentz'", True, [
        ('Amp', num(p[1]), True, '1.0',  '400.0', '0.001',  '1000.0'),
        ('Br',  num(p[2]), True, '0.10', '6.00',  '1.0E-4', '1000.0'),
        ('Eo',  num(p[3]), True, '3.50', '9.00',  '1.0E-4', '15.0'),
        ('Eg',  num(p[4]), True, '2.50', '5.50',  '0.0',    '15.0')])]

def ag_oscs(p):
    """p = [d, rough, Einf, rho, tau, (G_Amp, G_Br, G_En) x 3]"""
    GB = [('3.00','7.00'), ('1.20','4.00'), ('3.80','8.00')]   # En window per Gaussian
    oscs = [("'Drude(RT)'", False, [
        ('Resistivity (Ohm·cm)', num(p[3]), True, '1.0E-7','1.0E-1','1.0E-9','1000.0'),
        ('Scat. Time (fs)',           num(p[4]), True, '0.2',   '60.0',  '1.0E-6','1000.0')])]
    for i in range(3):
        lo_en, hi_en = GB[i]
        oscs.append(("'Gaussian'", False, [
            ('Amp',  num(p[5+3*i]), True,  '0.0','40.0',   '-100.0','1000.0'),
            ('iAmp', '0.0',         False, '-1000.0','1000.0','-1000.0','1000.0'),
            ('Br',   num(p[6+3*i]), True,  '0.10','5.00',  '0.0','100.0'),
            ('En',   num(p[7+3*i]), True,  lo_en, hi_en,   '1.0E-8','15.0')]))
    return p[2], oscs

def layer(tpl_block, idx, d_nm, fitblk, gradeblk, lo_nm, hi_nm, fit_d=True):
    b = tpl_block
    b = re.sub(r"(start_Gen-Osc Fit Parms\r\n).*?(\r\n\t\tend_Gen-Osc Fit Parms)",
               lambda m: m.group(1)+fitblk+m.group(2), b, count=1, flags=re.S)
    b = re.sub(r"(start_Gen-Osc Grade Parms\r\n).*?(\r\n\t\tend_Gen-Osc Grade Parms)",
               lambda m: m.group(1)+gradeblk+m.group(2), b, count=1, flags=re.S)
    b = re.sub(r"^\t\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tF\t'Thickness # \d+')",
               lambda m: '\t\t'+num(d_nm*10)+('\tT\t' if fit_d else '\tF\t')+num(lo_nm*10)+'\t'+num(hi_nm*10)+m.group(1),
               b, count=1, flags=re.M)
    b = b.replace('263.42946105543194\tT\t','0.0\tF\t',1).replace('11.015243042816298\tT\t','11.0\tF\t',1)
    b = b.replace('start_Layer2','start_Layer%d'%idx).replace('end_Layer2','end_Layer%d'%idx)
    b = b.replace("'Thickness # 2'", "'Thickness # %d'"%idx)
    return b

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    tpl = open(TEMPLATE,'rb').read().decode('utf-8-sig')
    Ls = tpl.split('\r\n')
    i0 = Ls.index('\tstart_Layer2'); i1 = Ls.index('\tend_Layer2')
    head = '\r\n'.join(Ls[:i0]); L2 = '\r\n'.join(Ls[i0:i1+1]); tail = '\r\n'.join(Ls[i1+1:])
    R = json.load(open('ag_final_result.json'))
    SEEDP = json.load(open('ag_seed_result.json'))
    seedp = {k:[v for v in [d for d in SEEDP[k] if abs(d['d_ox']-2.0)<1e-9][0]['p']] for k in ('HATCN','MoOx')}
    print('%-6s %-6s %4s %6s %8s %8s  %s'%('sheet','seed','Ag','layers','d_seed','d_Ag','file'))
    for sh in L.ORDER:
        name, dsn, dan = L.SPEC[sh]
        sp = seedp[name]                       # [d_seed, Einf, TLa, TLbr, TLeo, TLeg]
        e_s, o_s = seed_oscs(sp[1:])
        fb, gb = blocks(e_s, o_s)
        # seed thickness is a process constant and is degenerate with the Ag layer,
        # so it is only fittable on the bare-seed samples where it was determined
        lay = [layer(L2, 2, sp[0], fb, gb, max(sp[0]-3,0.5), sp[0]+3, fit_d=(dan==0))]
        if dan > 0:
            p = np.array(R[sh]['p'])
            e_a, o_a = ag_oscs(p)
            fb2, gb2 = blocks(e_a, o_a)
            lay.append(layer(L2, 3, p[0], fb2, gb2, max(p[0]-4,0.5), p[0]+4))
            rough = p[1]
        else:
            rough = 0.0
        t = head + '\r\n' + '\r\n\t\r\n'.join(lay) + '\r\n' + tail
        nlay = 1 + 1 + len(lay)                                   # substrate + NTVE + gen-osc layers
        t = re.sub(r"(start_Model Parms\r\n\t)\d+(\t0\t)", lambda m: m.group(1)+str(nlay)+m.group(2), t, count=1)
        t = re.sub(r"^\t[^\t\n]+(\t)[TF](\t)-5\.0\t5\.0(\tF\t'Angle Offset')",
                   lambda m: '\t0.0'+m.group(1)+'F'+m.group(2)+'-2.0\t2.0'+m.group(3), t, count=1, flags=re.M)
        t = re.sub(r"^\t[^\t\n]+(\t)[TF](\t)0\.0\t500\.0(\tF\t'Roughness')",
                   lambda m: '\t'+num(rough*10)+m.group(1)+('T' if rough>0.3 else 'F')+m.group(2)+'0.0\t150.0'+m.group(3),
                   t, count=1, flags=re.M)
        t = re.sub(r"^\t\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tT\t'Native Oxide')",
                   lambda m: '\t\t20.0\tF\t10.0\t100.0'+m.group(1), t, count=1, flags=re.M)
        fn = '%s_%s%d_Ag%d.mod'%(sh, name, dsn, dan)
        open(os.path.join(OUTDIR,fn),'w',encoding='utf-8-sig',newline='').write(t)
        print('%-6s %-6s %4d %6d %8.2f %8s  %s'%(sh,name,dan,nlay,sp[0],
              ('%.2f'%R[sh]['p'][0]) if dan>0 else '-', fn))
    print('saved ->', OUTDIR)

if __name__ == '__main__':
    main()
