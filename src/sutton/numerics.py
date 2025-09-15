import numpy as np


def thomas(aa, bb, cc, dd):
    """
    Tridiagonal matrix solver (Thomas algorithm).
    aa: lower diag, bb: main diag, cc: upper diag, dd: RHS
    Returns solution vector q.
    """
    n = len(bb)
    bet = np.zeros(n)
    gam = np.zeros(n)
    q = np.zeros(n)

    bet[0] = bb[0]
    gam[0] = dd[0] / bb[0]

    for i in range(1, n):
        bet[i] = bb[i] - (aa[i] * cc[i - 1] / bet[i - 1])
        gam[i] = (dd[i] - aa[i] * gam[i - 1]) / bet[i]

    q[n - 1] = gam[n - 1]
    for i in range(n - 2, -1, -1):
        q[i] = gam[i] - (cc[i] * q[i + 1] / bet[i])
    return q


def thomas_alt(aa, bb, cc, dd):
    """
    Alternate Thomas algorithm implementation.
    """
    n = len(bb)
    gam = np.zeros(n)
    q = np.zeros(n)

    bet = bb[0]
    q[0] = dd[0] / bet

    for i in range(1, n):
        gam[i] = cc[i - 1] / bet
        bet = bb[i] - aa[i] * gam[i]
        if bet == 0:
            break
        q[i] = (dd[i] - aa[i] * q[i - 1]) / bet

    for i in range(n - 2, -1, -1):
        q[i] = q[i] - gam[i + 1] * q[i + 1]
    return q


def our_central_difference(s, dz):
    """
    Central differences; forward/backward at edges.
    """
    m = len(s)
    dsdz = np.zeros(m)
    dsdz[1:m - 1] = (s[2:m] - s[0:m - 2]) / (2 * dz)
    dsdz[0] = (s[1] - s[0]) / dz
    dsdz[m - 1] = (s[m - 1] - s[m - 2]) / dz
    return dsdz
