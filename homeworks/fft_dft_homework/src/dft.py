"""
Discrete Fourier Transform Module
==================================
Contains loop-based and vectorized DFT implementations.
"""

import numpy as np


def dft_loop(x: np.ndarray) -> np.ndarray:
    """Compute the Discrete Fourier Transform using nested loops.

    Implements the DFT formula directly:
        X[k] = Σ_{n=0}^{N-1} x[n] * exp(-j 2π k n / N)

    Parameters
    ----------
    x : np.ndarray
        Input signal (1-D, real or complex).

    Returns
    -------
    X : np.ndarray
        DFT of the input signal (complex-valued).
    """
    N = len(x)
    X = np.zeros(N, dtype=complex)

    for k in range(N):
        for n in range(N):
            X[k] += x[n] * np.exp(-1j * 2 * np.pi * k * n / N)

    return X


def dft_vector(x: np.ndarray) -> np.ndarray:
    """Compute the Discrete Fourier Transform using matrix multiplication.

    Constructs the DFT matrix W where W[k, n] = exp(-j 2π k n / N)
    and computes X = W @ x.

    Parameters
    ----------
    x : np.ndarray
        Input signal (1-D, real or complex).

    Returns
    -------
    X : np.ndarray
        DFT of the input signal (complex-valued).
    """
    N = len(x)
    n = np.arange(N)
    k = n.reshape(-1, 1)  # column vector

    # DFT matrix
    W = np.exp(-1j * 2 * np.pi * k * n / N)

    X = W @ x
    return X
