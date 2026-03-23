"""STFT analysis utilities for chirp time-frequency experiments.

This module provides reusable functions for:
1) computing the short-time Fourier transform (STFT),
2) plotting spectrograms,
3) generating all required figures for the homework.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import stft

from signal_generation import generate_chirp


def compute_stft(
    x: np.ndarray,
    fs: int,
    *,
    window: str | tuple[str, float] = "hann",
    nperseg: int = 128,
    noverlap: int = 64,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute STFT of a 1D signal.

    STFT formula applied by ``scipy.signal.stft``:
        Zxx[k, m] = sum_{n=0}^{N-1} x[n] w[n - mR] exp(-j 2pi k n / N)

    where:
    - ``w`` is the chosen analysis window,
    - ``R`` is the hop size (R = nperseg - noverlap),
    - ``m`` indexes short-time frames,
    - ``k`` indexes frequency bins.
    """
    # SciPy returns:
    # f   : frequency bins in Hz
    # tau : frame-center times in seconds
    # Zxx : complex STFT matrix with shape (n_freq_bins, n_time_frames)
    f, tau, zxx = stft(
        x,
        fs=fs,
        window=window,
        nperseg=nperseg,
        noverlap=noverlap,
        boundary=None,
        padded=False,
    )
    return f, tau, zxx


def _half_power_width(axis_values: np.ndarray, profile: np.ndarray) -> float:
    """Return -3 dB (half-power) contiguous width around the peak.

    The half-power threshold is:
        A_half = A_peak / sqrt(2)

    Width is measured on the provided axis between nearest left/right points
    around the peak that remain above the threshold.
    """
    peak_idx = int(np.argmax(profile))
    peak_val = float(profile[peak_idx])
    if peak_val <= 0:
        return 0.0

    threshold = peak_val / np.sqrt(2.0)

    left = peak_idx
    while left > 0 and profile[left - 1] >= threshold:
        left -= 1

    right = peak_idx
    while right < profile.size - 1 and profile[right + 1] >= threshold:
        right += 1

    return float(axis_values[right] - axis_values[left])


def _energy_weighted_width(axis_values: np.ndarray, profile: np.ndarray) -> float:
    """Return robust width estimate using second central moment of energy.

    Energy weighting uses:
        w_i = |A_i|^2

    and width is converted from standard deviation to a diameter-like measure:
        width = 2 * sigma

    This avoids zero-width artifacts on coarse STFT grids.
    """
    energy = np.square(np.maximum(profile, 0.0))
    total = float(np.sum(energy))
    if total <= 0:
        return 0.0

    mu = float(np.sum(axis_values * energy) / total)
    variance = float(np.sum(np.square(axis_values - mu) * energy) / total)
    sigma = np.sqrt(max(variance, 0.0))
    return float(2.0 * sigma)


def quantify_resolution_blur(
    f: np.ndarray,
    tau: np.ndarray,
    zxx: np.ndarray,
    *,
    f0: float,
    f1: float,
    duration: float,
) -> dict[str, float]:
    """Quantify spectrogram blur in frequency and time.

    Frequency blur metric:
        Median frame-wise energy-weighted ridge thickness in Hz.

    Time blur metric:
        For several target frequencies, take the STFT magnitude along time at
        that frequency and measure energy-weighted temporal spread in seconds.

    Returned values are provided both in physical units and in bin/frame units
    to make the "blurry" concept directly measurable.
    """
    magnitude = np.abs(zxx)

    # Frequency blur: how thick the chirp ridge is in Hz at each time frame.
    freq_widths_hz: list[float] = []
    for frame_idx in range(magnitude.shape[1]):
        frame_profile = magnitude[:, frame_idx]
        peak = float(np.max(frame_profile))
        focus_mask = frame_profile >= (0.1 * peak)
        if np.any(focus_mask):
            freq_widths_hz.append(
                _energy_weighted_width(f[focus_mask], frame_profile[focus_mask])
            )
        else:
            freq_widths_hz.append(_energy_weighted_width(f, frame_profile))

    # Time blur: how wide the energy is in time for selected frequency slices.
    # Select frequencies away from boundaries for stable measurements.
    target_freqs = np.linspace(f0 + 0.1 * (f1 - f0), f1 - 0.1 * (f1 - f0), 5)
    time_widths_s: list[float] = []
    for target_freq in target_freqs:
        freq_idx = int(np.argmin(np.abs(f - target_freq)))
        time_profile = magnitude[freq_idx, :]
        # Use the full temporal profile at each chosen frequency slice to
        # avoid zero-width values when thresholded masks collapse to one frame.
        time_widths_s.append(_energy_weighted_width(tau, time_profile))

    delta_f = float(f[1] - f[0]) if f.size > 1 else 1.0
    delta_tau = float(tau[1] - tau[0]) if tau.size > 1 else duration

    freq_blur_hz = float(np.median(freq_widths_hz))
    time_blur_s = float(np.median(time_widths_s))

    return {
        "freq_blur_hz": freq_blur_hz,
        "time_blur_s": time_blur_s,
        "freq_blur_bins": freq_blur_hz / delta_f,
        "time_blur_frames": time_blur_s / delta_tau,
    }


def save_resolution_metrics_csv(
    rows: list[dict[str, float | int | str]],
    *,
    save_path: str | Path,
) -> None:
    """Save quantitative resolution metrics to CSV."""
    fieldnames = [
        "window",
        "nperseg",
        "noverlap",
        "freq_blur_hz",
        "time_blur_s",
        "freq_blur_bins",
        "time_blur_frames",
    ]
    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)

    with save_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def plot_spectrogram(
    f: np.ndarray,
    tau: np.ndarray,
    zxx: np.ndarray,
    *,
    title: str,
    save_path: str | Path,
    cmap: str = "magma",
) -> None:
    """Plot and save STFT magnitude spectrogram.

    The displayed magnitude uses:
        |Zxx(f, tau)| = sqrt(Re(Zxx)^2 + Im(Zxx)^2)

    and is converted to dB scale for readability:
        20 * log10(|Zxx| + eps)
    """
    # Magnitude of the complex STFT coefficients.
    magnitude = np.abs(zxx)
    eps = 1e-12
    magnitude_db = 20.0 * np.log10(magnitude + eps)

    fig, ax = plt.subplots(figsize=(10, 4.5))
    mesh = ax.pcolormesh(
        tau,
        f,
        magnitude_db,
        shading="gouraud",
        cmap=cmap,
    )
    ax.set_title(title)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Frequency [Hz]")
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label("Magnitude [dB]")
    ax.set_ylim(0, np.max(f))
    fig.tight_layout()

    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=180)
    plt.close(fig)


def plot_chirp_signal(
    t: np.ndarray,
    x: np.ndarray,
    *,
    save_path: str | Path,
) -> None:
    """Plot chirp signal in the time domain and save figure."""
    fig, ax = plt.subplots(figsize=(10, 3.8))
    ax.plot(t, x, color="tab:blue", linewidth=0.9)
    ax.set_title("Linear Chirp Signal in Time Domain")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=180)
    plt.close(fig)


def run_full_analysis() -> None:
    """Generate all required STFT figures for the homework."""
    root = Path(__file__).resolve().parents[1]
    figures_dir = root / "figures"
    data_dir = root / "data"

    # Signal setup used across all experiments.
    fs = 1000
    duration = 1.0
    f0 = 5.0
    f1 = 200.0
    t, x = generate_chirp(fs=fs, duration=duration, f0=f0, f1=f1)

    # 1) Time-domain chirp plot.
    plot_chirp_signal(t, x, save_path=figures_dir / "chirp_signal.png")

    # 2) Baseline STFT (Hann, nperseg=128, noverlap=64).
    f, tau, zxx = compute_stft(x, fs, window="hann", nperseg=128, noverlap=64)
    plot_spectrogram(
        f,
        tau,
        zxx,
        title="Baseline STFT Spectrogram (window=hann, nperseg=128, noverlap=64)",
        save_path=figures_dir / "stft_hann_baseline.png",
    )

    # 3) Window comparison with the same time/frequency tiling.
    # Gaussian window uses tuple format: ("gaussian", std_samples).
    window_specs: list[tuple[str, str | tuple[str, float]]] = [
        ("hann", "hann"),
        ("hamming", "hamming"),
        ("blackman", "blackman"),
        ("gaussian", ("gaussian", 16.0)),
    ]
    for win_name, win_spec in window_specs:
        f, tau, zxx = compute_stft(x, fs, window=win_spec, nperseg=128, noverlap=64)
        plot_spectrogram(
            f,
            tau,
            zxx,
            title=f"STFT Spectrogram ({win_name} window, nperseg=128, noverlap=64)",
            save_path=figures_dir / f"stft_{win_name}.png",
        )

    # 4) Resolution analysis by varying segment length with Hann window.
    # To keep overlap ratio consistent at 50%, noverlap = nperseg // 2.
    metrics_rows: list[dict[str, float | int | str]] = []
    for nperseg in [64, 128, 256]:
        f, tau, zxx = compute_stft(
            x,
            fs,
            window="hann",
            nperseg=nperseg,
            noverlap=nperseg // 2,
        )
        plot_spectrogram(
            f,
            tau,
            zxx,
            title=f"STFT Resolution Study (hann, nperseg={nperseg}, noverlap={nperseg // 2})",
            save_path=figures_dir / f"stft_nperseg_{nperseg}.png",
        )

        blur_metrics = quantify_resolution_blur(
            f,
            tau,
            zxx,
            f0=f0,
            f1=f1,
            duration=duration,
        )
        metrics_rows.append(
            {
                "window": "hann",
                "nperseg": nperseg,
                "noverlap": nperseg // 2,
                **blur_metrics,
            }
        )

    save_resolution_metrics_csv(
        metrics_rows,
        save_path=data_dir / "resolution_blur_metrics.csv",
    )


if __name__ == "__main__":
    run_full_analysis()
