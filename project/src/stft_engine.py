"""STFT from first principles (explicit DFT, explicit segmentation)."""

from __future__ import annotations

import numpy as np

from window_math import get_window, window_energy
from transform_logic import dft_one_sided, frequency_bins


def stft_psd(
    x: np.ndarray,
    fs: float,
    *,
    nperseg: int,
    hop: int,
    window_name: str = "hann",
    sigma: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute PSD per STFT frame.

    Returns (psd_frames, freqs, times) where psd_frames shape is (frames, n_freqs).
    """
    if nperseg <= 0:
        raise ValueError("nperseg must be positive.")
    if hop <= 0 or hop > nperseg:
        raise ValueError("hop must satisfy 1 <= hop <= nperseg.")

    x = np.asarray(x, dtype=float)
    if x.size < nperseg:
        raise ValueError("Signal is shorter than nperseg.")

    window = get_window(window_name, nperseg, sigma=sigma)
    U = window_energy(window)

    starts = np.arange(0, x.size - nperseg + 1, hop)
    n_frames = starts.size
    freqs = frequency_bins(nperseg, fs, one_sided=True)
    psd_frames = np.empty((n_frames, freqs.size), dtype=float)
    times = np.empty(n_frames, dtype=float)

    for i, start in enumerate(starts):
        frame = x[start : start + nperseg]
        Xw = dft_one_sided(frame * window)
        psd_frames[i, :] = (np.abs(Xw) ** 2) / (fs * U)
        times[i] = (start + (nperseg - 1) / 2.0) / fs

    return psd_frames, freqs, times
