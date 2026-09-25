"""Minimal OMML (Office Math) builder: equations that PowerPoint/Word open as native, editable equations."""
from xml.sax.saxutils import escape
M = 'http://schemas.openxmlformats.org/officeDocument/2006/math'
def r(t, sty=None):
    rp = f'<m:rPr><m:sty m:val="{sty}"/></m:rPr>' if sty else ''
    return f'<m:r>{rp}<m:t xml:space="preserve">{escape(t)}</m:t></m:r>'
def t(s): return r(s, 'p')                      # upright text (operators, numbers, labels)
def cat(*xs): return ''.join(xs)
def e(x): return f'<m:e>{x}</m:e>'
def sub(b, s): return f'<m:sSub>{e(b)}<m:sub>{s}</m:sub></m:sSub>'
def sup(b, s): return f'<m:sSup>{e(b)}<m:sup>{s}</m:sup></m:sSup>'
def subsup(b, s, p): return f'<m:sSubSup>{e(b)}<m:sub>{s}</m:sub><m:sup>{p}</m:sup></m:sSubSup>'
def frac(n, d): return f'<m:f><m:num>{n}</m:num><m:den>{d}</m:den></m:f>'
def par(x, l='(', rr=')'): return f'<m:d><m:dPr><m:begChr m:val="{l}"/><m:endChr m:val="{rr}"/></m:dPr>{e(x)}</m:d>'
def nary(ch, lo, hi, body):
    hide = '' if hi else '<m:supHide m:val="1"/>'
    return (f'<m:nary><m:naryPr><m:chr m:val="{ch}"/><m:limLoc m:val="undOvr"/>{hide}</m:naryPr>'
            f'<m:sub>{lo}</m:sub><m:sup>{hi or ""}</m:sup>{e(body)}</m:nary>')
def Sum(lo, body, hi=None): return nary('∑', lo, hi, body)
def Prod(lo, body, hi=None): return nary('∏', lo, hi, body)
def omath(x): return f'<m:oMath xmlns:m="{M}">{x}</m:oMath>'
def para(x): return f'<m:oMathPara xmlns:m="{M}"><m:oMath>{x}</m:oMath></m:oMathPara>'
# symbols used throughout
th = r('θ'); thi = sub(r('θ'), r('i')); thr = sub(r('θ'), r('r')); tht = sub(r('θ'), r('t'))
def Pk(k='k', arg=th): return cat(subsup(r('P'), t('sub'), par(r(k))), par(arg))
def BT(a=th): return cat(sub(r('B'), t('T')), par(a))
def BR(a, b): return cat(sub(r('B'), t('R')), par(cat(a, t(','), b)))
def RL(a=th): return cat(sub(r('R'), t('LED')), par(a))
eta = lambda s: sub(r('η'), t(s))
def etak(s, k='k'): return subsup(r('η'), t(s), par(r(k)))
pk = lambda k='k': sub(r('p'), r(k)); Ak = lambda k='k': sub(r('A′'), r(k))
cs = cat(t('cos '), th, t(' sin '), th)
