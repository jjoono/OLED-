"""CompleteEASE .mod writer, corrected-TMM fit.
Parameters whose 90% confidence interval exceeds 50% of their value are written
with the fit flag OFF - they carry no information in 340-1080 nm, and leaving
them free is what produced the huge error bars on Eo2/Eg2/Br2."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # measurement exports (.xlsx), CompleteEASE .mod
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fitted results, figures, intermediates

import re, os, json

TEMPLATE = _os.path.join(ELLIPS_DATA, r'VendorA_ITO_GenOSC_6functions.mod')
OUTDIR   = _os.path.join(ELLIPS_OUT, r'mod_v13')
RESULT   = _os.path.join(ELLIPS_OUT, r'ce_v13_result.json')
NAME = {'#1':'VendorA_ITO_1','#2':'VendorA_ITO_2','#3':'VendorB_IZO_1','#4':'VendorB_IZO_2','#5':'VendorA_ITO_2pctO2'}
TAIL = "\tF\tF\t0.0\t0.0\t100000.0\tF\t100.0\t"
num  = lambda x: repr(float(x))
pl   = lambda v,f,lo,hi,nm: "\t\t\t%s\t%s\t%s\t%s\tF\t'%s'%s"%(v,'T' if f else 'F',lo,hi,nm,TAIL)
gl   = lambda lo,hi,i: "\t\t\t0.0\tF\t%s\t%s\tF\t'Grade %d'%s"%(lo,hi,i,TAIL)

def blocks(p, free):
    """p=[rough,Einf,rho,tau,TL_A,TL_Br,TL_Eo,TL_Eg,(G_A,G_Br,G_En)]  free=list[bool]"""
    ng = len(p) > 8
    osc = [("'Drude(RT)'", False, [
              ('Resistivity (Ohm\u00b7cm)', num(p[2]), free[2], '1.0E-6','1.0',   '1.0E-9','1000.0'),
              ('Scat. Time (fs)',           num(p[3]), free[3], '0.5',   '50.0',  '1.0E-6','1000.0')]),
           ("'Tauc-Lorentz'", True, [
              ('Amp', num(p[4]), free[4], '5.0',    '1500.0', '0.001',  '1000.0'),
              ('Br',  num(p[5]), free[5], '0.10',   '6.50',   '1.0E-4', '1000.0'),
              ('Eo',  num(p[6]), free[6], '3.00',   '7.00',   '1.0E-4', '15.0'),
              ('Eg',  num(p[7]), free[7], '2.50',   '4.20',   '0.0',    '15.0')])]
    if ng:
        osc.append(("'Gaussian'", False, [
              ('Amp',  num(p[8]),  free[8],  '0.0',     '15.0',   '-100.0',  '1000.0'),
              ('iAmp', '0.0',      False,    '-1000.0', '1000.0', '-1000.0', '1000.0'),
              ('Br',   num(p[9]),  free[9],  '0.20',    '3.00',   '0.0',     '100.0'),
              ('En',   num(p[10]), free[10], '1.20',    '3.60',   '1.0E-8',  '15.0')]))
    fit = ["\t\t\t%d\t%s\t%s\t0.5\t5.0\tF\t'Einf'%s" % (len(osc), num(p[1]), 'T' if free[1] else 'F', TAIL)]
    grade = []
    for typ, common_eg, lst in osc:
        fit.append('\t\t\t'+typ+'\t')
        for nm, v, f, lo, hi, glo, ghi in lst:
            fit.append(pl(v, f, lo, hi, nm))
            grade.append('\t\t\tF\t'); grade += [gl(glo, ghi, j) for j in (1,2,3)]
        if common_eg: fit.append('\t\t\tF\t')
    fit.append('\t\t\t'); grade.append('\t\t\t')
    return '\r\n'.join(fit), '\r\n'.join(grade), len(osc)

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    tpl = open(TEMPLATE,'rb').read().decode('utf-8-sig')
    R = json.load(open(RESULT))
    print('%-16s %4s %6s %7s %7s %7s  %s'%('sample','osc','MSE','d(nm)','rough','dth','locked (fit flag OFF)'))
    for sh, rec in R.items():
        p = rec['p']; rel = rec['rel']
        free = [rel[j] <= 35.0 for j in range(len(p))]
        free[1] = free[1] and True
        fb, gb, n_osc = blocks(p, free)
        t = tpl
        t = re.sub(r"(start_Gen-Osc Fit Parms\r\n).*?(\r\n\t\tend_Gen-Osc Fit Parms)",
                   lambda m: m.group(1)+fb+m.group(2), t, count=1, flags=re.S)
        t = re.sub(r"(start_Gen-Osc Grade Parms\r\n).*?(\r\n\t\tend_Gen-Osc Grade Parms)",
                   lambda m: m.group(1)+gb+m.group(2), t, count=1, flags=re.S)
        t = re.sub(r"^\t[^\t\n]+(\t)[TF](\t)-5\.0\t5\.0(\tF\t'Angle Offset')",
                   lambda m: '\t'+num(rec['dth'])+m.group(1)+'F'+m.group(2)+'-2.0\t2.0'+m.group(3),
                   t, count=1, flags=re.M)
        rgfit = 'T' if p[0] > 0.3 else 'F'
        t = re.sub(r"^\t[^\t\n]+(\t)[TF](\t)0\.0\t500\.0(\tF\t'Roughness')",
                   lambda m: '\t'+num(p[0]*10)+m.group(1)+rgfit+m.group(2)+'0.0\t150.0'+m.group(3),
                   t, count=1, flags=re.M)
        t = re.sub(r"^\t\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tF\t'Thickness # 2')",
                   lambda m: '\t\t'+num(rec['d']*10)+'\tT\t%s\t%s'%(num((rec['d']-5)*10),num((rec['d']+5)*10))+m.group(1),
                   t, count=1, flags=re.M)
        t = t.replace('263.42946105543194\tT\t','0.0\tF\t',1)
        t = t.replace('11.015243042816298\tT\t','11.0\tF\t',1)
        open(os.path.join(OUTDIR, NAME[sh]+'_GenOsc_v12.mod'),'w',encoding='utf-8-sig',newline='').write(t)
        NMS=['rough','Einf','rho','tau','TL_Amp','TL_Br','TL_Eo','TL_Eg','G_Amp','G_Br','G_En']
        print('%-16s %4d %6.2f %7.1f %7.2f %+7.2f  %s'
              %(NAME[sh], n_osc, rec['mse'], rec['d'], p[0], rec['dth'],
                ', '.join(NMS[j] for j in range(len(p)) if not free[j]) or '-'))
    print('saved ->', OUTDIR)

if __name__ == '__main__':
    main()
