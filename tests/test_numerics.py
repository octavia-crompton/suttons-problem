import numpy as np
from sutton import thomas, our_central_difference

def test_thomas_small_system():
    # Solve tri-diagonal system for q: A q = d
    # A with lower, diag, upper: aa, bb, cc
    aa = np.array([0, -1, -1, -1], dtype=float)
    bb = np.array([2, 2, 2, 2], dtype=float)
    cc = np.array([-1, -1, -1, 0], dtype=float)
    d  = np.array([1, 0, 0, 1], dtype=float)
    q = thomas(aa, bb, cc, d)
    # Validate by reconstructing A @ q ~ d
    A = np.diag(bb) + np.diag(cc[:-1], 1) + np.diag(aa[1:], -1)
    assert np.allclose(A @ q, d)


def test_central_difference_edges():
    s = np.array([0.0, 1.0, 4.0, 9.0])
    dz = 1.0
    d = our_central_difference(s, dz)
    # Edge forward/backward; center should be (4-0)/2 = 2 at index 1
    assert d[0] == 1.0
    assert d[-1] == 5.0
    assert d[1] == 2.0
