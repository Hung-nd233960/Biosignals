import numpy as np
from src.dsp_utils import welch_psd, band_power


def test_single_sine_band_power():
    fs = 200.0
    duration = 4.0
    t = np.arange(0, duration, 1.0 / fs)
    A = 1.0
    f0 = 10.0
    x = A * np.sin(2 * np.pi * f0 * t)

    N = 256
    H = 64
    psd, freqs = welch_psd(x, fs, N, H, window_name="hann")

    # Theoretical total power of a pure sine of amplitude A is A^2/2
    # Integrate PSD around f0 (±1 Hz) to capture the peak.
    p = band_power(psd, freqs, (f0 - 1.0, f0 + 1.0))
    assert 0.3 < p < 0.7, f"band power {p} out of expected range"
