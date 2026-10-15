"""Generate CompleteEASE .mod files for the 5 samples from the fitted Gen-Osc
parameters, using the user's VendorA_ITO_GenOSC_v2.mod as the structural template.

Conversions applied:
  thickness / roughness : nm -> Angstrom (x10)
  Drude  A[eV^2], Br[eV] -> Resistivity[Ohm.cm], Scattering time[fs]
        tau[fs] = hbar[eV.fs]/Br ;  rho[Ohm.m] = hbar^2/(A*eps0*tau[s])
  Gaussian Br : my sigma -> CompleteEASE Br = 2*sqrt(ln2)*sigma
  Einf : my einf goes to 'Einf'; UV/IR poles set to 0 (my model has no pole term)
"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # measurement exports (.xlsx), CompleteEASE .mod
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fitted results, figures, intermediates

import json, re, os, math

HBAR_EVS = 6.582119569e-16      # eV*s
HBAR_EVFS = 0.6582119569        # eV*fs
EPS0 = 8.8541878128e-12         # F/m

TEMPLATE = _os.path.join(ELLIPS_DATA, r'VendorA_ITO_GenOSC_v2.mod')
OUTDIR = _os.path.join(ELLIPS_OUT, r'mod_files')
PARAMS = _os.path.join(ELLIPS_OUT, r'genosc_params.json')

NAMEMAP = {'1':'VendorA_ITO_1','2':'VendorA_ITO_2','3':'VendorB_IZO_1',
           '4':'VendorB_IZO_2','5':'VendorA_ITO_2pctO2'}

def drude_to_ce(A_eV2, Br_eV):
    """A[eV^2], Br[eV] -> (resistivity Ohm*cm, scattering time fs)."""
    tau_fs = HBAR_EVFS / Br_eV
    tau_s  = tau_fs * 1e-15
    rho_m  = HBAR_EVS**2 / (A_eV2 * EPS0 * tau_s)   # Ohm*m
    return rho_m * 100.0, tau_fs                     # Ohm*cm, fs

def fmt(v):
    return repr(float(v))

def par(value, fit, lo, hi, name):
    """One CompleteEASE parameter line (without leading tabs)."""
    f = 'T' if fit else 'F'
    return "%s\t%s\t%s\t%s\tF\t'%s'\tF\tF\t0.0\t0.0\t100000.0\tF\t100.0\t" % (
        fmt(value), f, fmt(lo), fmt(hi), name)

def build_genosc_block(p, indent='\t\t\t'):
    einf   = p[0]
    A_D    = p[1]
    A_TL, E0_TL, C_TL, Eg_TL = p[2], p[3], p[4], p[5]
    A_G, Ec_G, Br_G          = p[6], p[7], p[8]
    rho, tau = drude_to_ce(A_D, 0.10)
    br_ce = Br_G * 2.0 * math.sqrt(math.log(2.0))     # sigma -> CompleteEASE Br

    L = []
    # leading "3" = number of oscillators, then Einf on the same line
    L.append("3\t" + par(einf, True, 0.0, 10.0, 'Einf'))
    L.append("'Drude(RT)'\t")
    L.append(par(rho, True, 1.0e-5, 1.0e-1, 'Resistivity (Ohm\u00b7cm)'))
    L.append(par(tau, False, 1.0, 20.0, 'Scat. Time (fs)'))
    L.append("'Tauc-Lorentz'\t")
    L.append(par(A_TL, True, 1.0, 400.0, 'Amp'))
    L.append(par(C_TL, True, 0.05, 10.0, 'Br'))
    L.append(par(E0_TL, True, 2.0, 8.0, 'Eo'))
    L.append(par(Eg_TL, True, 2.0, 5.0, 'Eg'))
    L.append("'Gaussian'\t")
    L.append(par(A_G, True, 0.0, 10.0, 'Amp'))
    L.append(par(br_ce, True, 0.05, 10.0, 'Br'))
    L.append(par(Ec_G, True, 0.3, 6.0, 'En'))
    L.append('F\t')
    L.append('')
    return '\r\n'.join(indent + x if x else '' for x in L)

def read_text(path):
    b = open(path, 'rb').read()
    for enc in ('utf-8-sig', 'utf-8', 'cp949', 'latin-1'):
        try:
            return b.decode(enc), enc
        except UnicodeDecodeError:
            continue
    raise RuntimeError('cannot decode')

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    tpl, enc = read_text(TEMPLATE)
    print('template encoding:', enc, '| lines:', tpl.count('\n')+1)
    R = json.load(open(PARAMS))

    for s in '12345':
        d = R[s]
        p = d['p']
        txt = tpl

        # 1) angle offset (global model parm)
        txt, n1 = re.subn(r"^\t[^\t\n]+(\t[TF]\t-5\.0\t5\.0\tF\t'Angle Offset')",
                          lambda m: '\t' + fmt(d['dth']) + m.group(1),
                          txt, count=1, flags=re.M)
        # 2) roughness  (nm -> Angstrom)
        txt, n2 = re.subn(r"^\t[^\t\n]+(\t[TF]\t0\.0\t500\.0\tF\t'Roughness')",
                          lambda m: '\t' + fmt(d['rough']*10.0) + m.group(1),
                          txt, count=1, flags=re.M)
        # 3) Layer2 physical thickness (nm -> Angstrom), widen bounds
        txt, n3 = re.subn(r"^\t\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tF\t'Thickness # 2')",
                          lambda m: '\t\t' + fmt(d['d']*10.0) + '\tT\t100.0\t2000.0' + m.group(1),
                          txt, count=1, flags=re.M)
        # 4) Gen-Osc oscillator block
        blk = build_genosc_block(p)
        txt, n4 = re.subn(r"(start_Gen-Osc Fit Parms\r\n).*?(\r\n\t\tend_Gen-Osc Fit Parms)",
                          lambda m: m.group(1) + blk + m.group(2),
                          txt, count=1, flags=re.S)
        # 5) poles -> 0 (this model carries no pole term; Einf absorbs it)
        txt, n5 = re.subn(r"^\t\t\t[^\t\n]+(\t[TF]\t-1000\.0\t1000\.0\tF\t'UV Pole Amp\.')",
                          lambda m: '\t\t\t0.0' + m.group(1).replace('\tT\t', '\tF\t', 1),
                          txt, count=1, flags=re.M)

        assert all([n1, n2, n3, n4, n5]), (s, n1, n2, n3, n4, n5)
        out = os.path.join(OUTDIR, NAMEMAP[s] + '_GenOsc.mod')
        open(out, 'w', encoding=enc if enc != 'utf-8-sig' else 'utf-8-sig',
             newline='').write(txt)
        rho, tau = drude_to_ce(p[1], 0.10)
        print('%-16s d=%6.1f A  rough=%5.1f A  dth=%+.2f  rho=%.4g ohm-cm  tau=%.3f fs  -> %s'
              % (NAMEMAP[s], d['d']*10, d['rough']*10, d['dth'], rho, tau,
                 os.path.basename(out)))

if __name__ == '__main__':
    main()
