import numpy as np


def padit(x, nz, upwind, array):
    xx = np.concatenate([np.arange(-40, 0, 2), x])
    up = upwind.reshape(1, nz)
    up = np.concatenate([up] * 10)
    array = np.concatenate([up, up, array])
    return xx, array
