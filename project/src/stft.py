"""STFT utilities that avoid high-level signal processing helpers."""

from __future__ import annotations

import numpy as np

from windows import get_window


def compute_window_length(fs: float, f_min: float, cycles: int = 3) -> int:
    """Return a window length that spans at least the requested cycles.

    N = ceil(cycles * fs / f_min)
    """
    if fs <= 0 or f_min <= 0:
        raise ValueError("fs and f_min must be positive.")
    if cycles <= 0:
        raise ValueError("cycles must be positive.")
    return int(np.ceil(cycles * fs / f_min))


def stft_power(
    x: np.ndarray,
    fs: float,
    *,
    window_name: str,
    nperseg: int,
    noverlap: int,
    fmin: float = 4.0,
    fmax: float = 100.0,
    window_kwargs: dict[str, float] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute a magnitude-squared STFT (power) with explicit steps.

    Steps
    -----
    1) Segment signal into length-N frames with hop = N - noverlap.
    2) Apply a window (taper) from first principles.
    3) Compute FFT using numpy.fft.rfft.
    4) Keep frequencies in [fmin, fmax].
    5) Return power = |FFT|^2 at frame center times.
    """
    if nperseg <= 0:
        raise ValueError("nperseg must be positive.")
    if noverlap < 0 or noverlap >= nperseg:
        raise ValueError("noverlap must satisfy 0 <= noverlap < nperseg.")
    if fmin < 0 or fmax <= fmin:
        raise ValueError("fmax must be greater than fmin.")

    x = np.asarray(x, dtype=float)
    hop = nperseg - noverlap
    if x.size < nperseg:
        raise ValueError("Signal is shorter than nperseg.")

    window_kwargs = window_kwargs or {}
    window = get_window(window_name, nperseg, **window_kwargs)

    n_frames = 1 + (x.size - nperseg) // hop
    times = np.empty(n_frames, dtype=float)

    freqs = np.fft.rfftfreq(nperseg, d=1.0 / fs)
    mask = (freqs >= fmin) & (freqs <= fmax)
    freqs = freqs[mask]

    power = np.empty((freqs.size, n_frames), dtype=float)

    for frame in range(n_frames):
        start = frame * hop
        segment = x[start : start + nperseg]
        tapered = segment * window
        spectrum = np.fft.rfft(tapered)
        power[:, frame] = np.abs(spectrum[mask]) ** 2
        center = start + (nperseg - 1) / 2.0
        times[frame] = center / fs

    return freqs, times, power
