import numpy as np


def stability(xi):
    """Brutsaert stability function (stable/unstable)."""
    xi = np.asarray(xi)
    phi = np.zeros_like(xi)
    inds = np.where(xi <= 1 / 16)
    phi[inds] = (1 - 16 * xi[inds]) ** -0.5
    inds = np.where(xi > 1 / 16)
    phi[inds] = 1 + 5.2 * xi[inds]
    return phi
