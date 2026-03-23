"""
Signal Generation Module
========================
Generates a chirp signal for DFT / FFT analysis.
"""

import numpy as np
from scipy.signal import chirp
import matplotlib.pyplot as plt
from pathlib import Path


def generate_chirp(
    A: float = 1.0,
    duration: float = 1.0,
    fs: int = 1000,
    f0: float = 5.0,
    f1: float = 200.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate a linear chirp signal.

    Parameters
    ----------
    A : float
        Amplitude of the chirp signal.
    duration : float
        Duration in seconds.
    fs : int
        Sampling frequency in Hz.
    f0 : float
        Start frequency in Hz.
    f1 : float
        End frequency in Hz.

    Returns
    -------
    t : np.ndarray
        Time vector.
    x : np.ndarray
        Chirp signal samples.
    """
    t = np.arange(0, duration, 1 / fs)
    x = A * chirp(t, f0=f0, f1=f1, t1=duration, method="linear")
    return t, x


def plot_chirp(
    t: np.ndarray,
    x: np.ndarray,
    save_path: str | Path | None = None,
) -> None:
    """Plot a chirp signal in the time domain.

    Parameters
    ----------
    t : np.ndarray
        Time vector.
    x : np.ndarray
        Chirp signal samples.
    save_path : str or Path, optional
        If provided, save the figure to this path.
    """
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(t, x, linewidth=0.6)
    ax.set_title("Chirp Signal (f₀ = 5 Hz → f₁ = 200 Hz)")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150)
        print(f"Figure saved to {save_path}")

    return fig, ax
