"""Phase 2 sparsity measures on a power map P (512 range x 256 Doppler)."""
import numpy as np


def k_frac(P, frac=0.9):
    """Number of bins holding frac of the total energy (largest first)."""
    p = np.sort(P.ravel())[::-1]
    c = np.cumsum(p) / p.sum()
    return int(np.searchsorted(c, frac)) + 1


def sparsity_report(P, fracs=(0.5, 0.8, 0.9, 0.99)):
    """Share of bins needed for each energy fraction, and a rough floor share."""
    p = np.sort(P.ravel())[::-1]
    c = np.cumsum(p) / p.sum()
    out = {f: (int(np.searchsorted(c, f)) + 1) / P.size for f in fracs}
    floor_share = np.median(P) * P.size / P.sum()
    return out, floor_share


def n_peaks(P, db):
    """Bins more than db above the median of the whole map."""
    return int((P > np.median(P) * 10 ** (db / 10)).sum())


def row_peaks(P, db=10):
    """Peaks per range row: bins more than db above that row's median. Returns 512 counts."""
    med = np.median(P, axis=1, keepdims=True)
    return (P > med * 10 ** (db / 10)).sum(axis=1)


def noise_map(seed=0, shape=(512, 256, 16)):
    """Noise control: complex Gaussian cube through the same windows, FFTs and channel sum."""
    from .rd import range_doppler
    rng = np.random.default_rng(seed)
    cube = rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
    return np.sum(np.abs(range_doppler(cube)) ** 2, axis=2)

def dense_peaks(Pr, thr):
    """Local maxima of one row's power that are above thr."""
    N = len(Pr)
    return [f for f in range(N) if Pr[f] > thr and Pr[f] >= Pr[(f - 1) % N] and Pr[f] >= Pr[(f + 1) % N]]