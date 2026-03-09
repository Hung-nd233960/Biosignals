"""
Sawtooth Wave – Fourier Series Reconstruction
==============================================
Reference signal via scipy.signal.sawtooth and harmonic synthesis
using all-harmonic sine summation with 1/k decay.
"""

import numpy as np
from scipy import signal


def sawtooth_reference(t: np.ndarray, f0: float = 5.0, A: float = 1.0) -> np.ndarray:
    """Generate a reference sawtooth wave using SciPy.

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
        Sawtooth wave samples.
    """
    return A * signal.sawtooth(2 * np.pi * f0 * t)


def sawtooth_fourier(
    t: np.ndarray, N: int, f0: float = 5.0, A: float = 1.0
) -> np.ndarray:
    """Reconstruct a sawtooth wave using its Fourier series.

    The SciPy reference ``sawtooth(2πf₀t)`` rises from −1 to +1. The
    matching Fourier series is:

        x(t) = −(2A/π) Σ_{k=1}^{N} (1/k) sin(2π k f₀ t)

    Parameters
    ----------
    t : np.ndarray
        Time vector.
    N : int
        Maximum harmonic index.
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
    for k in range(1, N + 1):  # all harmonics: 1, 2, 3, ...
        x += (1.0 / k) * np.sin(2 * np.pi * k * f0 * t)
    x *= -(2.0 * A) / np.pi
    return x
