"""The recycling series in its matrix form, Eq. (1) of the paper.

    eta_ext = sum_j  B_T . v_j ,    v_1 = P_sub (normalised),
                                    v_{j+1} = diag(R_LED) B_R v_j

B_T and B_R come from the ray trace, R_LED from the dipole model.  Unlike the
closed form eta_ext = p/[p+(1-p)A'], this does not assume that the extraction
surface randomises the angles at every bounce, which a shallow lens does not.
"""
import numpy as np


def eta_ext(BT, BR, R_LED, Psub, n_term=60):
    v = np.asarray(Psub, float).copy()
    v /= v.sum()
    M = (np.asarray(BR, float) * np.asarray(R_LED, float)[:, None])
    out, terms = 0.0, []
    for _ in range(n_term):
        t = float(BT @ v)
        terms.append(t)
        out += t
        v = M @ v
    return out, np.array(terms)
