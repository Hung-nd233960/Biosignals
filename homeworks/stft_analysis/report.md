# STFT Time-Frequency Analysis of a Chirp Signal

**Course:** Biomedical Signal Processing  
**Topic:** Short-Time Fourier Transform (STFT) of a non-stationary chirp

---

## 1. Introduction

A chirp signal is a signal whose frequency changes over time. In this assignment, we use a **linear chirp** that sweeps from low to high frequency.

A standard Fourier transform gives global frequency content, but it does not indicate **when** each frequency occurs. For non-stationary signals like chirps, we need a time-frequency method. STFT solves this by computing Fourier transforms over short, shifted windows.

---

## 2. Method

### 2.1 STFT idea

STFT divides the signal into short overlapping segments. Each segment is multiplied by a window and transformed into frequency domain:

$$
Z_{xx}(f_k, \tau_m) = \sum_{n=0}^{N-1} x[n] \, w[n-mR] \, e^{-j2\pi kn/N}
$$

Where:

- $x[n]$: input signal
- $w[\cdot]$: analysis window
- $N$: window length (`nperseg`)
- $R = N - \text{noverlap}$: hop size
- $\tau_m$: time position of frame $m$
- $f_k$: frequency bin $k$
- $Z_{xx}(f_k, \tau_m)$: complex STFT coefficient

The plotted spectrogram uses magnitude:

$$
|Z_{xx}| = \sqrt{\Re(Z_{xx})^2 + \Im(Z_{xx})^2}
$$

In figures, magnitude is shown in decibels:

$$
20\log_{10}(|Z_{xx}| + \varepsilon)
$$

with small $\varepsilon$ to avoid $\log(0)$.

### 2.2 Baseline parameters

- Sampling frequency: $f_s = 1000$ Hz
- Baseline STFT: `window='hann'`, `nperseg=128`, `noverlap=64`

---

## 3. Spectral Windows

Windowing reduces spectral leakage by tapering segment edges before Fourier transform.

- **Hann**: smooth and widely used general-purpose window.
- **Hamming**: similar to Hann with slightly different sidelobe behavior; often slightly better leakage suppression in some cases.
- **Blackman**: stronger sidelobe attenuation, usually less leakage but wider main lobe.
- **Gaussian**: smooth bell-shaped taper, often providing a clean compromise between leakage suppression and localization.

In this implementation, the comparison figures use Hann, Hamming, Blackman, and Gaussian with the same `nperseg=128`, `noverlap=64`.

---

## 4. Results

### 4.1 Chirp Signal

Time-domain plot (`figures/chirp_signal.png`) shows an oscillation that becomes denser over time, indicating increasing instantaneous frequency.

### 4.2 Window Comparison

Compared figures:

- `figures/stft_hann.png`
- `figures/stft_hamming.png`
- `figures/stft_blackman.png`
- `figures/stft_gaussian.png`

Observations:

- All windows capture the upward chirp trend clearly.
- Hann and Hamming give similar overall behavior.
- Blackman often appears smoother with reduced leakage energy around the main chirp trace, at the cost of slight broadening.
- Gaussian gives a smooth ridge with good leakage control and balanced visual localization.

### 4.3 Resolution Analysis (Window Size Trade-off)

Compared figures:

- `figures/stft_nperseg_64.png`
- `figures/stft_nperseg_128.png`
- `figures/stft_nperseg_256.png`

Key trade-off:

- **Small window (`nperseg=64`)**:
  - better time localization (rapid changes easier to localize in time)
  - poorer frequency resolution (thicker/blurry frequency trace)
- **Medium window (`nperseg=128`)**:
  - balanced time and frequency resolution
- **Large window (`nperseg=256`)**:
  - better frequency resolution (sharper frequency localization)
  - poorer time resolution (events spread over time)

This is the standard STFT uncertainty trade-off.

### 4.4 Quantitative Meaning of "Blurry"

To quantify blur (instead of only visual judgment), the script computes two metrics and saves them in:

- `data/resolution_blur_metrics.csv`

Definitions:

1. **Frequency blur** (ridge thickness in Hz)

For each time frame, take STFT magnitude along frequency and compute an energy-weighted spread:

$$
w_i = |Z_{xx}(f_i, \tau)|^2,
\quad
\mu_f = \frac{\sum_i f_i w_i}{\sum_i w_i},
\quad
\sigma_f = \sqrt{\frac{\sum_i (f_i-\mu_f)^2 w_i}{\sum_i w_i}}
$$

Then define width as:

$$
W_f = 2\sigma_f
$$

Report value = median of $W_f$ over frames.

1. **Time blur** (temporal spread in seconds)

For several fixed frequency slices, take magnitude along time and compute:

$$
w_m = |Z_{xx}(f, \tau_m)|^2,
\quad
\mu_\tau = \frac{\sum_m \tau_m w_m}{\sum_m w_m},
\quad
\sigma_\tau = \sqrt{\frac{\sum_m (\tau_m-\mu_\tau)^2 w_m}{\sum_m w_m}}
$$

Then define:

$$
W_\tau = 2\sigma_\tau
$$

Report value = median of $W_\tau$ over selected frequency slices.

Measured values (Hann window):

|nperseg|noverlap|freq_blur_hz|time_blur_s|freq_blur_bins|time_blur_frames|
|---|---:|---:|---:|---:|---:|
|64|32|18.150|0.0923|1.162|2.883|
|128|64|11.228|0.0573|1.437|0.895|
|256|128|14.610|0.0304|3.740|0.237|

Interpretation:

- Time blur decreases clearly as `nperseg` increases (from 0.0923 s to 0.0304 s in this metric setting).
- Frequency blur in bins increases strongly for larger `nperseg`, showing broader ridge occupancy in discrete-bin terms.
- These numeric widths provide a concrete definition of "blurry": larger width means blurrier representation along that axis.

---

## 5. Image Resolution Interpretation

A spectrogram can be read as an image:

- x-axis: time ($\tau$)
- y-axis: frequency ($f$)
- color: magnitude (energy level)

For this chirp, the main energy appears as a rising diagonal ridge. Interpretation:

- A clear thin ridge means better frequency concentration.
- A thick or blurry ridge suggests lower frequency precision or stronger leakage.
- Temporal sharpness (how quickly features move along x-axis) depends on window length.

---

## 6. Conclusion

STFT is effective for analyzing non-stationary signals such as chirps because it preserves both time and frequency information.

Main conclusions:

- Window choice affects spectral leakage and visual smoothness.
- Window length controls the time-frequency resolution trade-off.
- The chirp behavior is clearly tracked as an increasing-frequency ridge in all spectrograms.

---

## 7. Deliverables

- Time-domain chirp plot: `figures/chirp_signal.png`
- Baseline STFT: `figures/stft_hann_baseline.png`
- Window comparison: `figures/stft_hann.png`, `figures/stft_hamming.png`, `figures/stft_blackman.png`, `figures/stft_gaussian.png`
- Resolution study: `figures/stft_nperseg_64.png`, `figures/stft_nperseg_128.png`, `figures/stft_nperseg_256.png`
- Quantitative metrics CSV: `data/resolution_blur_metrics.csv`
- Reusable code module: `src/stft_analysis.py`
- Notebook: `notebook/stft_analysis.ipynb`
- Report: `report.md`
