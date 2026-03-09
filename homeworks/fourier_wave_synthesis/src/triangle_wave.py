"""
Triangle Wave – Fourier Series Reconstruction
==============================================
Reference signal via scipy.signal.sawtooth(width=0.5) and harmonic
synthesis using odd-harmonic sine summation with 1/k² decay.
"""

import numpy as np
from scipy import signal


def triangle_reference(t: np.ndarray, f0: float = 5.0, A: float = 1.0) -> np.ndarray:
    """Generate a reference triangle wave using SciPy.

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
        Triangle wave samples.
    """
    return A * signal.sawtooth(2 * np.pi * f0 * t, width=0.5)


def triangle_fourier(
    t: np.ndarray, N: int, f0: float = 5.0, A: float = 1.0
) -> np.ndarray:
    """Reconstruct a triangle wave using its Fourier series.

    The SciPy reference ``sawtooth(2πf₀t, width=0.5)`` produces a
    cosine-phased triangle (value −1 at t = 0). The matching Fourier
    series uses **cosine** terms:

        x(t) = −(8A/π²) Σ_{k=1,3,5,…,N} (1/k²) cos(2π k f₀ t)

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
        x += (1.0 / k**2) * np.cos(2 * np.pi * k * f0 * t)
    x *= -(8.0 * A) / (np.pi**2)
    return x
