"""Text helpers for report generation inside notebooks."""

from __future__ import annotations


def theoretical_background_markdown() -> str:
    """Return a markdown block for the theoretical background section."""
    return (
        "## 1. Theoretical Background\n"
        "Short-Time Fourier Transform (STFT) analyzes nonstationary signals by \n"
        "windowing the signal and computing a local Fourier transform.\n\n"
        "The discrete-time STFT is:\n\n"
        "$$\\n"
        "X(m, k) = \\sum_{n=0}^{N-1} x[n + mR] \, w[n] \, e^{-j 2\\pi k n / N}\\n"
        "$$\n\n"
        "where $N$ is the window length, $R$ is hop size, and $w[n]$ is a taper.\n\n"
        "Power is computed as:\n\n"
        "$$\\n"
        "P(m, k) = |X(m, k)|^2\\n"
        "$$\n"
    )
