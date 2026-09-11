"""What actually determines the absorption of a CLOSED ultrathin Ag film.

Not thickness and roughness.  Matthiessen splits the resistivity -- and so the
Drude damping, and so the absorption -- into two independent terms:

    rho = rho_bulk + d_rho_surf(d, p) + d_rho_gb(D)

    Fuchs-Sondheimer   d_rho_surf/rho_bulk = 0.375 (1-p) l / d
    Mayadas-Shatzkes   d_rho_gb/rho_bulk  ~= 1.5 alpha,  alpha = (l/D) R/(1-R)

p is the surface specularity and D the lateral grain size.  Roughness enters
only through p, and rms alone does not fix p -- a surface with 1 nm rms and a
100 nm correlation length is near-specular, the same rms at a 2 nm correlation
length is fully diffuse.  D is independent of the surface entirely.

So two films with identical d and identical rms can differ in absorption by as
much as the grain-size term allows, which for these samples is a factor of 1.5.
"""
import numpy as np

RHO_B, MFP, WP, G0, LAM = 1.59, 52.0, 9.17, 0.021, 550.0
R_GB = 0.4          # grain-boundary reflection coefficient for Ag

JC = np.array([[397.4,0.05,2.07],[450.9,0.04,2.66],[495.9,0.05,3.09],
               [520.9,0.05,3.34],[548.6,0.06,3.59],[582.1,0.05,3.93],
               [616.8,0.06,4.15],[659.5,0.05,4.48],[704.5,0.041,4.84]])


def ag(lam, gamma):
    hw = 1239.84/lam
    e_ib = (np.interp(lam, JC[:,0], JC[:,1])+1j*np.interp(lam, JC[:,0], JC[:,2]))**2 \
           + WP**2/(hw**2 + 1j*G0*hw)
    return np.sqrt(e_ib - WP**2/(hw**2 + 1j*gamma*hw))


def tr(n, d, lam):
    n = [complex(x) for x in n]; k0 = 2*np.pi/lam
    M = np.eye(2, dtype=complex)
    for j in range(len(n)-1):
        r = (n[j]-n[j+1])/(n[j]+n[j+1]); t = 2*n[j]/(n[j]+n[j+1])
        I = np.array([[1,r],[r,1]], dtype=complex)/t
        if j+1 < len(n)-1:
            dl = k0*n[j+1]*d[j+1]
            I = I @ np.array([[np.exp(-1j*dl),0],[0,np.exp(1j*dl)]])
        M = M @ I
    return (n[-1].real/n[0].real)*abs(1/M[0,0])**2, abs(M[1,0]/M[0,0])**2


def rho_of(d, p, D):
    surf = 0.375*(1-p)*MFP/d
    gb = 1.5*(MFP/D)*R_GB/(1-R_GB) if D and np.isfinite(D) else 0.0
    return RHO_B*(1 + surf + gb)


def A_glass(d, p, D):
    g = G0*rho_of(d, p, D)/RHO_B
    T, R = tr([1.0, ag(LAM, g), 1.75, 1.52], [0, d, 5.0, 0], LAM)
    return 100*(1-T-R)


def A_dev(d, p, D):
    g = G0*rho_of(d, p, D)/RHO_B
    T, R = tr([1.8, ag(LAM, g), 2.10, 1.0], [0, d, 60.0, 0], LAM)
    return 100*(1-T-R)


# ---- what grain size do the measured films imply?
print("Decomposing the measured HATCN / Ag 5 nm film (Rs = 23.3 ohm/sq)\n")
rho_meas = 23.3*5*0.1
surf = RHO_B*0.375*MFP/5.0                       # p = 0, the worst case
gb = rho_meas - RHO_B - surf
alpha = (gb/RHO_B)/1.5
D_impl = MFP*(R_GB/(1-R_GB))/alpha
print(f"  measured rho          {rho_meas:6.2f} uOhm-cm")
print(f"  bulk                  {RHO_B:6.2f}")
print(f"  surfaces (p = 0)      {surf:6.2f}   <- cannot be larger than this")
print(f"  grain boundaries      {gb:6.2f}   = {100*gb/rho_meas:.0f} % of the total")
print(f"  -> implied grain size D = {D_impl:.0f} nm\n")

print("A at fixed d = 5 nm.  Rows: surface quality.  Columns: grain size.")
print("Same thickness, and the same rms could give any row.\n")
Ds = [15, 21, 30, 50, 100, np.inf]
print(f"  {'p':>5} |" + "".join(f"{('D='+(str(x) if np.isfinite(x) else 'inf')):>10}" for x in Ds))
print("  " + "-"*(7+10*len(Ds)))
for p in [0.0, 0.3, 0.6, 0.9]:
    print(f"  {p:5.1f} |" + "".join(f"{A_glass(5.0,p,D):9.2f}%" for D in Ds))

print("\n  same table, in the device (organic / Ag 5 / CPL 60 / air):")
print(f"  {'p':>5} |" + "".join(f"{('D='+(str(x) if np.isfinite(x) else 'inf')):>10}" for x in Ds))
print("  " + "-"*(7+10*len(Ds)))
for p in [0.0, 0.3, 0.6, 0.9]:
    print(f"  {p:5.1f} |" + "".join(f"{A_dev(5.0,p,D):9.2f}%" for D in Ds))

print("\nAt the measured point (p = 0, D = 21 nm), the two levers are worth:")
base = A_glass(5.0, 0.0, D_impl)
print(f"  baseline                          {base:5.2f} %")
print(f"  double the grain size, 21 -> 42   {A_glass(5.0,0.0,2*D_impl):5.2f} %"
      f"   ({A_glass(5.0,0.0,2*D_impl)-base:+.2f})")
print(f"  half-specular surfaces, p 0 -> 0.5 {A_glass(5.0,0.5,D_impl):4.2f} %"
      f"   ({A_glass(5.0,0.5,D_impl)-base:+.2f})")
print(f"  both                              {A_glass(5.0,0.5,2*D_impl):5.2f} %"
      f"   ({A_glass(5.0,0.5,2*D_impl)-base:+.2f})")
