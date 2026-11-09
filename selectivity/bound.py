"""Closed-form ceiling of multi-pass, angle-selective extraction (ideal filter, full randomization):
eta_beta = p_beta / [p_beta + (1 - p_beta) A'],  p_beta = sin^2(beta) / n_sub^2  (Lambertian, AR-ideal)
and the forward (on-axis) luminance in Lambertian-EQE units, eta_beta / sin^2(beta), relative to the
same randomizer without a filter, eta_0 = p0 / [p0 + (1 - p0) A'] with p0 = T_bar / n_sub^2."""
import numpy as np


def p_lamb(n):
    """Lambertian escape probability of a flat, uncoated n -> air interface (Yablonovitch T_bar / n^2)."""
    th = np.linspace(0, np.pi / 2, 20001)[:-1] + np.pi / 4e4
    s = n * np.sin(th); ok = s < 1; ct = np.sqrt(1 - s[ok] ** 2); c = np.cos(th[ok])
    rs = ((n * c - ct) / (n * c + ct)) ** 2; rp = ((n * ct - c) / (n * ct + c)) ** 2
    w = np.cos(th) * np.sin(th)
    return float(((1 - 0.5 * (rs + rp)) * w[ok]).sum() / w.sum())


def eta_cone(A, beta, n):
    pb = np.sin(np.radians(beta)) ** 2 / n ** 2
    return pb / (pb + (1 - pb) * A)


def eta_open(A, n):
    p0 = p_lamb(n)
    return p0 / (p0 + (1 - p0) * A)


def fwd_gain(A, beta, n):
    return eta_cone(A, beta, n) / np.sin(np.radians(beta)) ** 2 / eta_open(A, n)
