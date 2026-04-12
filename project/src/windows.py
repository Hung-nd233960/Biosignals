"""Window (taper) functions implemented from first principles.

All windows return a length-N numpy array. These functions avoid using
scipy.signal.windows or any other high-level window generators.
"""

from __future__ import annotations

import numpy as np


def _validate_length(n: int) -> int:
    """Validate and normalize window length.

    Parameters
    ----------
    n
        Window length in samples.

    Returns
    -------
    int
        Validated window length.
    """
    if n <= 0:
        raise ValueError("Window length must be positive.")
    return int(n)


def rectangular(n: int) -> np.ndarray:
    """Return a rectangular (square) window.

    w[n] = 1
    """
    n = _validate_length(n)
    return np.ones(n, dtype=float)


def hann(n: int) -> np.ndarray:
    """Return a Hann window.

    w[n] = 0.5 * (1 - cos(2*pi*n/(N-1)))
    """
    n = _validate_length(n)
    if n == 1:
        return np.ones(1, dtype=float)
    idx = np.arange(n, dtype=float)
    return 0.5 * (1.0 - np.cos(2.0 * np.pi * idx / (n - 1.0)))


def hamming(n: int, a0: float = 0.54, a1: float = 0.46) -> np.ndarray:
    """Return a Hamming window with configurable coefficients.

    w[n] = a0 - a1 * cos(2*pi*n/(N-1))
    """
    n = _validate_length(n)
    if n == 1:
        return np.ones(1, dtype=float)
    idx = np.arange(n, dtype=float)
    return a0 - a1 * np.cos(2.0 * np.pi * idx / (n - 1.0))


def blackman(n: int, a0: float = 0.42, a1: float = 0.5, a2: float = 0.08) -> np.ndarray:
    """Return a Blackman window with configurable coefficients.

    w[n] = a0 - a1 * cos(2*pi*n/(N-1)) + a2 * cos(4*pi*n/(N-1))
    """
    n = _validate_length(n)
    if n == 1:
        return np.ones(1, dtype=float)
    idx = np.arange(n, dtype=float)
    angle = 2.0 * np.pi * idx / (n - 1.0)
    return a0 - a1 * np.cos(angle) + a2 * np.cos(2.0 * angle)


def gaussian(n: int, sigma: float = 0.4) -> np.ndarray:
    """Return a Gaussian window with configurable width.

    The standard form is:
        w[n] = exp(-0.5 * ((n - (N-1)/2) / (sigma*(N-1)/2))^2)

    Parameters
    ----------
    n
        Window length in samples.
    sigma
        Normalized width in [0, 0.5]. Use larger values for wider windows.
        If sigma == 0, this function returns a rectangular window.
    """
    n = _validate_length(n)
    if sigma < 0 or sigma > 0.5:
        raise ValueError("sigma must be in [0, 0.5].")
    if n == 1:
        return np.ones(1, dtype=float)
    if sigma == 0:
        return np.ones(n, dtype=float)
    idx = np.arange(n, dtype=float)
    center = (n - 1.0) / 2.0
    denom = sigma * center
    return np.exp(-0.5 * ((idx - center) / denom) ** 2)


def get_window(name: str, n: int, **kwargs: float) -> np.ndarray:
    """Return a window by name.

    Parameters
    ----------
    name
        Window name: "rectangular", "hann", "hamming", "blackman", "gaussian".
    n
        Window length in samples.
    **kwargs
        Optional window coefficients.
    """
    name = name.lower()
    if name in {"rectangular", "square", "boxcar"}:
        return rectangular(n)
    if name == "hann":
        return hann(n)
    if name == "hamming":
        return hamming(n, **kwargs)
    if name == "blackman":
        return blackman(n, **kwargs)
    if name == "gaussian":
        return gaussian(n, **kwargs)
    raise ValueError(f"Unsupported window name: {name}")
