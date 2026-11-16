"""Exactly-sparse FFT (Hassanieh, Indyk, Katabi, Price, STOC 2012) for one Doppler row."""
import numpy as np


def tone(N, f, amp=1.0):
    """Test signal: one tone at bin f (f may be fractional for off-grid tests)."""
    n = np.arange(N)
    return amp * np.exp(2j * np.pi * f * n / N)


def permute(x, sigma):
    """Read x at sigma*m mod N. Moves every tone from bin f to bin sigma*f mod N."""
    N = len(x)
    m = np.arange(N)
    return x[(sigma * m) % N]

from scipy.special import erf


def flat_window(N, B, edge=1.5, thr=1e-3):
    """Window whose spectrum is ~1 over L = N/B neighbouring bins and ~0 elsewhere.

    edge: how rounded the band edges are, in bins (larger = shorter window, less flat).
    thr:  time samples smaller than thr * max are set to 0 (makes the window short).
    Returns the time window w (length N, mostly zeros) and the designed band shape G.
    """
    L = N // B
    d = np.fft.fftfreq(N, 1 / N)        # signed bin offset of each bin: 0, 1, ..., 127, -128, ..., -1
    # a box of width L centred on 0, with edges rounded over about `edge` bins
    G = 0.5 * (erf((d + L / 2) / (edge * np.sqrt(2))) - erf((d - L / 2) / (edge * np.sqrt(2))))
    w = np.real(np.fft.ifft(G)) * N   # back to time (np.fft.ifft divides by N, so undo that)
    w[np.abs(w) < thr * np.abs(w).max()] = 0
    return w, G

def hash_to_bins(x, sigma, tau, w, B):
    """Bucket values u[0..B-1]: bucket j ~ sum of permuted bins near j*N/B.

    Reads x at sigma*(m + tau) mod N (permutation + shift), multiplies by the
    flat window w, folds into B, and does one B-point FFT.
    """
    N = len(x)
    m = np.arange(N)
    y = x[(sigma * (m + tau)) % N]            # permute and shift by tau
    z = y * w                                  # window: only the nonzero w samples matter
    folded = z.reshape(N // B, B).sum(axis=0)   # cut into pieces of B, add them
    return np.fft.fft(folded)                  # one small B-point FFT


def locate(x, sigma, w, B, thresh, ratio_tol=0.1):
    """One round: find the tones that sit alone in a bucket.

    Returns {true bin f: estimated DFT value X[f]}.
    """
    N = len(x)
    L = N // B
    Gw = np.real(np.fft.fft(w)) / N              # band shape of the real window
    u0 = hash_to_bins(x, sigma, 0, w, B)
    u1 = hash_to_bins(x, sigma, 1, w, B)
    sigma_inv = pow(sigma, -1, N)

    found = {}
    for j in range(B):
        if abs(u0[j]) < thresh:                  # empty bucket
            continue
        r = u1[j] / u0[j]
        if abs(abs(r) - 1) > ratio_tol:          # |ratio| far from 1: a collision
            continue
        p = int(np.round(np.angle(r) % (2 * np.pi) * N / (2 * np.pi))) % N   # permuted bin
        off = ((p - j * L + N // 2) % N) - N // 2                            # distance from bucket centre
        if abs(off) > L // 2:                    # p is not inside bucket j's band: a collision
            continue
        f = (p * sigma_inv) % N                        # undo the permutation
        found[f] = u0[j] / Gw[(p - j * L) % N]   # undo the band shape
    return found

def sparse_fft(x, B, w, rounds=12, tol=1e-6, rng=None, stats=None):     # NEW: stats
    """Exactly-sparse FFT. Returns {bin f: X[f]} for the nonzero bins found.

    Each round: random odd sigma, locate lone tones, subtract them, repeat on the rest.
    If stats is a dict, the number of rounds used is written into it.
    """
    N = len(x)
    rng = rng if rng is not None else np.random.default_rng()
    thresh = tol * np.linalg.norm(x) * np.sqrt(N)

    est = {}
    r = x.copy()
    used = 0                                                           # NEW
    for _ in range(rounds):
        if np.linalg.norm(r) < tol * np.linalg.norm(x):
            break
        sigma = int(rng.choice(np.arange(1, N, 2)))
        found = locate(r, sigma, w, B, thresh)
        used += 1                                                      # NEW
        for f, v in found.items():
            est[f] = est.get(f, 0) + v
            r = r - tone(N, f, v / N)
    if stats is not None:                                              # NEW
        stats["rounds"] = used                                         # NEW
    return {f: v for f, v in est.items() if abs(v) > thresh}

def dense_ops(N):
    """Complex multiplications of a radix-2 N-point FFT."""
    return (N // 2) * int(np.log2(N))


def sparse_ops(N, B, w, rounds):
    """Complex multiplications of sparse_fft: 2 hashings per round + 2 per bucket in locate."""
    per_hashing = np.count_nonzero(w) + (B // 2) * int(np.log2(B))
    return rounds * (2 * per_hashing + 2 * B)