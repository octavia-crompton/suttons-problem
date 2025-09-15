import numpy as np
from .numerics import thomas, our_central_difference


def integrate_H2O_implicit(n, m, dx, dz, A, B, C, Qup, Qs, Qa):
    AA1 = -A * B
    AA2 = -C * B
    AA3 = 1 / dx
    AA4 = Qup / dx

    upd = (AA1 / (dz ** 2) + AA2 / (2 * dz))
    dia = (-2 * AA1 / (dz ** 2) + AA3)
    lod = (AA1 / (dz ** 2) - AA2 / (2 * dz))

    co = np.zeros(m)
    co[:] = AA4

    lod[0] = lod[m - 1] = 0
    dia[0] = dia[m - 1] = 1
    upd[0] = upd[m - 1] = 0

    co[0] = Qs
    co[m - 1] = Qa

    Q1 = thomas(lod, dia, upd, co)
    dQdz = our_central_difference(Q1, dz)

    Fq = np.zeros_like(Q1)
    Fq[1:] = -0.5 * (A[1:] + A[:-1]) * (0 * dQdz[:-1] + dQdz[1:]) * 0.5
    Fq[1] = Fq[2]
    Fq[0] = Fq[1]
    return Q1, Fq


essential = integrate_H2O_implicit

def integrate_T_implicit(n, m, dx, dz, A, B, C, Tup, Ts, Ta):
    AA1 = -A * B
    AA2 = -C * B
    AA3 = 1 / dx
    AA4 = Tup / dx

    upd = (AA1 / (dz ** 2) + AA2 / (2 * dz))
    dia = (-2 * AA1 / (dz ** 2) + AA3)
    lod = (AA1 / (dz ** 2) - AA2 / (2 * dz))

    co = np.zeros(m)
    co[:] = AA4

    lod[0] = lod[m - 1] = 0
    dia[0] = dia[m - 1] = 1
    upd[0] = upd[m - 1] = 0

    co[0] = Ts
    co[-1] = Ta

    T1 = thomas(lod, dia, upd, co)
    dTdz = our_central_difference(T1, dz)

    FT = np.zeros_like(T1)
    FT[1:] = -0.5 * (A[1:] + A[:-1]) * (0 * dTdz[:-1] + dTdz[1:]) * 0.5
    FT[1] = FT[2]
    FT[0] = FT[1]
    return T1, FT
