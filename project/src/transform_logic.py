"""DFT and PSD computations from explicit summations (no FFT)."""

from __future__ import annotations

import math
from typing import Tuple

import numpy as np


def dft(x: np.ndarray) -> np.ndarray:
    """Compute the full DFT of x using nested loops."""
    x = np.asarray(x, dtype=float)
    N = x.size
    X = np.empty(N, dtype=complex)
    for k in range(N):
        real = 0.0
        imag = 0.0
        for n in range(N):
            angle = -2.0 * math.pi * k * n / N
            real += x[n] * math.cos(angle)
            imag += x[n] * math.sin(angle)
        X[k] = complex(real, imag)
    return X


def dft_one_sided(x: np.ndarray) -> np.ndarray:
    """Return one-sided DFT for real-valued x (bins 0..N//2)."""
    X = dft(x)
    N = X.size
    return X[: N // 2 + 1]


def frequency_bins(N: int, fs: float, one_sided: bool = True) -> np.ndarray:
    if N <= 0 or fs <= 0:
        raise ValueError("N and fs must be positive.")
    if one_sided:
        n_bins = N // 2 + 1
    else:
        n_bins = N
    freqs = np.empty(n_bins, dtype=float)
    for k in range(n_bins):
        freqs[k] = k * fs / N
    return freqs


def periodogram(
    frame: np.ndarray, fs: float, window: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute one-sided PSD using the DFT definition.

    PSD(f_k) = |X_w[k]|^2 / (fs * U)
    where U = sum(w[n]^2)
    """
    if frame.size != window.size:
        raise ValueError("window and frame must have same length")
    xw = frame * window
    U = float(np.sum(window * window))
    X = dft_one_sided(xw)
    psd = (np.abs(X) ** 2) / (fs * U)
    freqs = frequency_bins(frame.size, fs, one_sided=True)
    return psd, freqs
