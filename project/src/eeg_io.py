"""EEG data access helpers."""

from __future__ import annotations

from pathlib import Path

import mne
import numpy as np


def load_channel_edf(
    path: str | Path,
    *,
    channel_name: str | None = None,
    preferred_channels: list[str] | None = None,
    duration_s: float | None = None,
) -> tuple[np.ndarray, float, str]:
    """Load one EEG channel from an EDF file.

    Parameters
    ----------
    path
        Path to EDF file.
    channel_name
        Explicit channel name to load if present.
    preferred_channels
        Ordered list of channel names to try (first match wins).
    duration_s
        Optional duration in seconds to load from the start of the file.

    Returns
    -------
    data, fs, channel_name
        Signal samples, sampling rate (Hz), and channel name.
    """
    path = Path(path)
    raw = mne.io.read_raw_edf(path.as_posix(), preload=True, verbose="ERROR")
    fs = float(raw.info["sfreq"])
    ch_names = list(raw.ch_names)

    if channel_name and channel_name in ch_names:
        pick = ch_names.index(channel_name)
    elif preferred_channels:
        pick = 0
        for name in preferred_channels:
            if name in ch_names:
                pick = ch_names.index(name)
                break
    else:
        pick = 0

    stop = None
    if duration_s is not None:
        if duration_s <= 0:
            raise ValueError("duration_s must be positive.")
        stop = int(round(duration_s * fs))

    data = raw.get_data(picks=[pick], start=0, stop=stop).squeeze().astype(float)
    channel_name = ch_names[pick]
    return data, fs, channel_name


def load_first_channel_edf(path: str | Path) -> tuple[np.ndarray, float, str]:
    """Load the first EEG channel from an EDF file (legacy helper)."""
    return load_channel_edf(path)
