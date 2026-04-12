"""Analysis helpers for overlap and window sensitivity studies."""

from __future__ import annotations

import numpy as np

from stft import stft_power


def overlap_sweep(
    x: np.ndarray,
    fs: float,
    *,
    window_name: str,
    nperseg: int,
    overlap_fracs: list[float],
    fmin: float = 4.0,
    fmax: float = 100.0,
    window_kwargs: dict[str, float] | None = None,
) -> dict[float, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Compute STFT power for multiple overlap ratios.

    Returns a mapping of overlap fraction -> (freqs, times, power).
    """
    results: dict[float, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for frac in overlap_fracs:
        if frac < 0 or frac >= 1:
            raise ValueError("overlap fraction must be in [0, 1).")
        noverlap = int(round(frac * nperseg))
        freqs, times, power = stft_power(
            x,
            fs,
            window_name=window_name,
            nperseg=nperseg,
            noverlap=noverlap,
            fmin=fmin,
            fmax=fmax,
            window_kwargs=window_kwargs,
        )
        results[frac] = (freqs, times, power)
    return results
