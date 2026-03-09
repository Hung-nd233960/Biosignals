"""
Metrics Module
==============
Error metrics for waveform reconstruction accuracy.
"""

import numpy as np


def rms_error(signal_a: np.ndarray, signal_b: np.ndarray) -> float:
    """Compute the Root Mean Square error between two signals.

    Formula:
        RMS = sqrt( mean( (signal_a - signal_b)² ) )

    Parameters
    ----------
    signal_a : np.ndarray
        First signal (e.g. reference).
    signal_b : np.ndarray
        Second signal (e.g. reconstruction).

    Returns
    -------
    float
        RMS error value.
    """
    return float(np.sqrt(np.mean((signal_a - signal_b) ** 2)))
