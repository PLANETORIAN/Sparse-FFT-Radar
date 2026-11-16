"""Comparing the sparse FFT with the dense FFT on real range rows."""
import numpy as np

from .sfft import sparse_fft_mc
from .sparsity import dense_peaks


def row_recall(din_row, P_row, thr, B, w, c=4.0, rng=None):
    """Run the multi-channel sparse FFT on one range row and count the dense peaks it keeps.

    A dense peak counts as kept if the sparse output is above thr within +-1 bin.
    Returns (peaks kept, peaks in the row, rounds used).
    """
    N = len(P_row)
    pk = dense_peaks(P_row, thr)
    st = {}
    est = sparse_fft_mc(din_row, B, w, c=c, rng=rng, stats=st)
    Ps = np.zeros(N)
    for f, v in est.items():
        Ps[f] = np.sum(np.abs(v) ** 2)
    got = Ps > thr
    hit = sum(got[(f - 1) % N] or got[f] or got[(f + 1) % N] for f in pk)
    return hit, len(pk), st["rounds"]