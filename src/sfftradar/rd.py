"""Range-Doppler processing of one RADIal frame, and CA-CFAR detection."""
import numpy as np
from scipy.signal import convolve2d

from .radial import frame_cube

N_RANGE, N_DOPPLER = 512, 256
W_RANGE   = 0.54 - 0.46 * np.cos(2 * np.pi * np.arange(N_RANGE)   / (N_RANGE - 1))    # Hamming
W_DOPPLER = 0.54 - 0.46 * np.cos(2 * np.pi * np.arange(N_DOPPLER) / (N_DOPPLER - 1))  # Hamming


def range_doppler(cube):
    """Complex (512, 256, 16) cube -> complex range-Doppler cube, same shape."""
    rng = np.fft.fft(cube * W_RANGE[:, None, None],  axis=0)
    return np.fft.fft(rng * W_DOPPLER[None, :, None], axis=1)


def frame_to_rd(seq, n):
    """Unfolded power map (512, 256): sum over the 16 channels of |X|^2. Used for sparsity."""
    cube = frame_cube(seq, n)
    cube = cube - cube.mean(axis=(0, 1), keepdims=True)        # DC offset per channel
    return np.sum(np.abs(range_doppler(cube)) ** 2, axis=2)


def frame_to_pred(seq, n, mti=False):
    """Magnitude map folded over the 12 DDM slots (512, 16). Used for detection (Phase 1)."""
    cube = frame_cube(seq, n)
    if mti:
        cube = cube - cube.mean(axis=1, keepdims=True)          # remove the mean over chirps
    else:
        cube = cube - cube.mean(axis=(0, 1), keepdims=True)
    Pm = np.sum(np.abs(range_doppler(cube)), axis=2)
    return Pm.reshape(N_RANGE, 16, 16).sum(axis=1)


def ca_cfar(P, n_rng=9, n_dop=9, g_rng=3, g_dop=3, thr_db=4):
    """Cell-averaging CFAR. True where a cell is thr_db above the mean of its training cells."""
    pw = np.abs(P) ** 2
    mask = np.ones((2 * n_rng + 1, 2 * n_dop + 1))
    mask[n_rng - g_rng: n_rng + g_rng + 1, n_dop - g_dop: n_dop + g_dop + 1] = 0

    sums   = convolve2d(pw, mask, mode="same")
    counts = convolve2d(np.ones_like(pw), mask, mode="same")
    return pw / (sums / counts) > 10 ** (thr_db / 10)