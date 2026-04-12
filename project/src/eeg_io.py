"""EEG data access helpers."""

from __future__ import annotations

from pathlib import Path

import mne
import numpy as np


def load_first_channel_edf(path: str | Path) -> tuple[np.ndarray, float, str]:
    """Load the first EEG channel from an EDF file.

    Parameters
    ----------
    path
        Path to EDF file.

    Returns
    -------
    data, fs, channel_name
        Signal samples, sampling rate (Hz), and channel name.
    """
    path = Path(path)
    raw = mne.io.read_raw_edf(path.as_posix(), preload=True, verbose="ERROR")
    data = raw.get_data(picks=[0]).squeeze().astype(float)
    fs = float(raw.info["sfreq"])
    channel_name = raw.ch_names[0]
    return data, fs, channel_name
