"""
Square Wave – Fourier Series Reconstruction
============================================
Reference signal via scipy.signal.square and harmonic synthesis
using odd-harmonic sine summation.
"""

import numpy as np
from scipy import signal


def square_reference(t: np.ndarray, f0: float = 5.0, A: float = 1.0) -> np.ndarray:
    """Generate a reference square wave using SciPy.

    Parameters
    ----------
    t : np.ndarray
        Time vector.
    f0 : float
        Fundamental frequency in Hz.
    A : float
        Amplitude.

    Returns
    -------
    np.ndarray
        Square wave samples.
    """
    return A * signal.square(2 * np.pi * f0 * t)


def square_fourier(
    t: np.ndarray, N: int, f0: float = 5.0, A: float = 1.0
) -> np.ndarray:
    """Reconstruct a square wave using its Fourier series.

    Formula:
        x(t) = (4A/π) Σ_{k=1,3,5,...,N} (1/k) sin(2π k f0 t)

    Parameters
    ----------
    t : np.ndarray
        Time vector.
    N : int
        Maximum harmonic index (odd harmonics up to N are used).
    f0 : float
        Fundamental frequency in Hz.
    A : float
        Amplitude.

    Returns
    -------
    np.ndarray
        Reconstructed waveform.
    """
    x = np.zeros_like(t, dtype=float)
    for k in range(1, N + 1, 2):  # odd harmonics: 1, 3, 5, ...
        x += (1.0 / k) * np.sin(2 * np.pi * k * f0 * t)
    x *= (4.0 * A) / np.pi
    return x
