"""Signal generation utilities used by the STFT analysis suite."""

from __future__ import annotations

import numpy as np


def linear_chirp(
    duration: float, fs: float, f0: float, f1: float
) -> tuple[np.ndarray, np.ndarray]:
    """Generate a linear chirp signal using an analytic phase model.

    The instantaneous frequency is:
        f(t) = f0 + k * t,  k = (f1 - f0) / duration

    The phase is:
        phi(t) = 2*pi*(f0*t + 0.5*k*t^2)

    Returns
    -------
    t, x
        Time vector (s) and chirp signal.
    """
    if duration <= 0:
        raise ValueError("duration must be positive.")
    if fs <= 0:
        raise ValueError("fs must be positive.")
    t = np.arange(0.0, duration, 1.0 / fs, dtype=float)
    k = (f1 - f0) / duration
    phase = 2.0 * np.pi * (f0 * t + 0.5 * k * t**2)
    x = np.sin(phase)
    return t, x
