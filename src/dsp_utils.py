"""DSP utilities implementing the SOP's contracts: DFT/STFT/PSD with correct window normalization.

Functions:
- window_*: window generators
- periodogram: computes PSD for one frame with U-normalization
- welch_psd: computes Welch PSD by averaging periodograms
- compute_stft: simple STFT (returns complex spectra, freqs, times)
- band_power: integrates PSD over a frequency band

All PSD outputs are in power per Hz (physical units) following:
PSD(f_k) = |X_w[k]|^2 / (f_s * U)
where U = sum_n w[n]^2

"""

from typing import Tuple, Literal
import numpy as np

WindowName = Literal["rect", "hann", "gaussian"]


def rect_window(N: int) -> np.ndarray:
    return np.ones(N, dtype=float)


def hann_window(N: int) -> np.ndarray:
    n = np.arange(N)
    return 0.5 * (1 - np.cos(2 * np.pi * n / (N - 1)))


def gaussian_window(N: int, std: float | None = None) -> np.ndarray:
    n = np.arange(N)
    mu = (N - 1) / 2.0
    if std is None:
        std = N / 8.0
    return np.exp(-0.5 * ((n - mu) / std) ** 2)


def get_window(name: WindowName, N: int, **kwargs) -> np.ndarray:
    if name == "rect":
        return rect_window(N)
    if name == "hann":
        return hann_window(N)
    if name == "gaussian":
        return gaussian_window(N, kwargs.get("std", None))
    raise ValueError(f"Unknown window: {name}")


def periodogram(
    frame: np.ndarray, fs: float, window: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute one-sided PSD (power/Hz) of `frame` using window `window`.

    Returns (psd, freqs) where psd length is N//2+1 for even N (one-sided).
    Uses normalization U = sum(window**2) and freq bin width implicitly via `fs`.
    Formula: PSD(f_k) = |X_w[k]|^2 / (fs * U)
    """
    N = len(frame)
    if len(window) != N:
        raise ValueError("window and frame must have same length")
    w = window.astype(float)
    U = np.sum(w * w)
    xw = frame * w
    X = np.fft.rfft(xw, n=N)
    psd = (np.abs(X) ** 2) / (fs * U)
    freqs = np.fft.rfftfreq(N, 1.0 / fs)
    return psd, freqs


def welch_psd(
    x: np.ndarray, fs: float, N: int, H: int, window_name: WindowName = "hann"
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute Welch PSD according to SOP: frame length N, hop H, window `window_name`.

    Returns (psd_avg, freqs).
    """
    w = get_window(window_name, N)
    L = len(x)
    if L < N:
        raise ValueError("Signal length must be >= segment length N")
    starts = np.arange(0, L - N + 1, H)
    psd_acc = None
    for s in starts:
        frame = x[s : s + N]
        p_seg, freqs = periodogram(frame, fs, w)
        if psd_acc is None:
            psd_acc = np.zeros_like(p_seg)
        psd_acc += p_seg
    psd_avg = psd_acc / len(starts)
    return psd_avg, freqs


def compute_stft(
    x: np.ndarray, fs: float, N: int, H: int, window_name: WindowName = "hann"
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Simple STFT: returns (STFT, freqs, times).

    STFT is shape (frames, n_freqs) with complex values (one-sided rfft results).
    times are frame centers in seconds.
    """
    w = get_window(window_name, N)
    L = len(x)
    starts = np.arange(0, L - N + 1, H)
    frames = []
    for s in starts:
        frame = x[s : s + N]
        X = np.fft.rfft(frame * w, n=N)
        frames.append(X)
    STFT = np.vstack(frames)
    freqs = np.fft.rfftfreq(N, 1.0 / fs)
    times = (starts + (N - 1) / 2.0) / fs
    return STFT, freqs, times


def band_power(psd: np.ndarray, freqs: np.ndarray, band: Tuple[float, float]) -> float:
    """Integrate PSD over `band` using trapezoidal integration to produce power (not power/Hz).

    The PSD is assumed to be power per Hz; integrating over Hz returns power in that band.
    """
    f0, f1 = band
    mask = (freqs >= f0) & (freqs <= f1)
    if not np.any(mask):
        return 0.0
    return float(np.trapz(psd[mask], freqs[mask]))
