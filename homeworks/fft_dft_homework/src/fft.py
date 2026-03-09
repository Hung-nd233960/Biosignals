"""
Fast Fourier Transform Module
==============================
Contains a recursive Cooley–Tukey FFT implementation.
"""

import numpy as np


def fft_recursive(x: np.ndarray) -> np.ndarray:
    """Compute the FFT using the recursive Cooley–Tukey algorithm.

    The input length **must** be a power of 2.  The algorithm splits the
    signal into even- and odd-indexed samples, recursively computes their
    FFTs, and combines the results using twiddle factors.

    Parameters
    ----------
    x : np.ndarray
        Input signal (1-D, real or complex).  Length must be a power of 2.

    Returns
    -------
    X : np.ndarray
        FFT of the input signal (complex-valued).

    Raises
    ------
    ValueError
        If the length of *x* is not a power of 2.
    """
    N = len(x)

    if N == 0:
        return np.array([], dtype=complex)

    # Check that N is a power of 2
    if N & (N - 1) != 0:
        raise ValueError(f"Input length must be a power of 2, got {N}.")

    # Base case
    if N == 1:
        return np.array(x, dtype=complex)

    # Recursive split
    X_even = fft_recursive(x[0::2])
    X_odd = fft_recursive(x[1::2])

    # Twiddle factors
    half = N // 2
    twiddle = np.exp(-1j * 2 * np.pi * np.arange(half) / N)

    # Combine (butterfly)
    X = np.empty(N, dtype=complex)
    X[:half] = X_even + twiddle * X_odd
    X[half:] = X_even - twiddle * X_odd

    return X
