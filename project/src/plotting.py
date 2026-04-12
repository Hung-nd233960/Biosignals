"""Plotting helpers for time-frequency analysis."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def plot_spectrogram(
    freqs: np.ndarray,
    times: np.ndarray,
    power: np.ndarray,
    *,
    title: str,
    save_path: str | Path,
    vmin_db: float | None = None,
    vmax_db: float | None = None,
) -> None:
    """Plot a power spectrogram in dB and save it to disk."""
    eps = 1e-12
    power_db = 10.0 * np.log10(power + eps)

    fig, ax = plt.subplots(figsize=(10, 4.5))
    mesh = ax.pcolormesh(times, freqs, power_db, shading="auto")
    ax.set_title(title)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Frequency (Hz)")
    if vmin_db is not None or vmax_db is not None:
        mesh.set_clim(vmin_db, vmax_db)
    fig.colorbar(mesh, ax=ax, label="Power (dB)")
    fig.tight_layout()

    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=180)
    plt.close(fig)
