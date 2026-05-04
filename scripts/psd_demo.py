"""Demo script: generate synthetic EEG-like signal and compute Welch PSD using `src.dsp_utils`.

Run:
python scripts/psd_demo.py

Produces a PNG `psd_demo.png` in repo root.
"""

import numpy as np
import matplotlib.pyplot as plt
from src.dsp_utils import welch_psd, band_power


def main():
    fs = 250.0
    t = np.arange(0, 10.0, 1.0 / fs)
    # two sinusoids + broadband noise
    s1 = 0.8 * np.sin(2 * np.pi * 10.0 * t)  # 10 Hz (alpha)
    s2 = 0.5 * np.sin(2 * np.pi * 20.0 * t)  # 20 Hz (beta)
    noise = 0.4 * np.random.randn(len(t))
    x = s1 + s2 + noise

    N = 512
    H = 128
    psd, freqs = welch_psd(x, fs, N, H, window_name="hann")

    plt.figure(figsize=(8, 4))
    plt.semilogy(freqs, psd)
    plt.title("Welch PSD (demo)")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Power / Hz")
    plt.xlim(0, 60)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("psd_demo.png")
    print("Saved psd_demo.png")

    # band powers
    alpha = band_power(psd, freqs, (8, 13))
    beta = band_power(psd, freqs, (13, 30))
    print(f"Alpha power (8-13 Hz): {alpha:.4f}")
    print(f"Beta power (13-30 Hz): {beta:.4f}")


if __name__ == "__main__":
    main()
