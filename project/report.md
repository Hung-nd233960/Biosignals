# EEG Preprocessing and Fourier Analysis Report

## 1. Overview

This report summarizes two connected analyses:

1. EEG preprocessing of real EDF data to produce analysis-ready signals.
2. Frequency-domain analysis using manual DTFT, FFT, IFFT, and stationarity comparison.

Primary notebooks:

- project/notebooks/preprocessing.ipynb
- project/notebooks/FFT.ipynb

Primary dataset:

- project/data/sub-NORB00055_ses-1_task-EEG_eeg.edf

---

## 2. Data and Experimental Setup

### 2.1 Raw data

- File: sub-NORB00055_ses-1_task-EEG_eeg.edf
- Sampling rate: 200 Hz
- Duration: 1140.00 s
- Total samples per channel: 228000
- Channels after cleaning: 21
- Dropped channels: 25+, 26+, 27+

### 2.2 Exported preprocessed outputs

The preprocessing notebook exports:

- project/data/preprocessed/epochs_rereferenced.npy
- project/data/preprocessed/epochs.npy
- project/data/preprocessed/epochs_with_buffer.npy
- project/data/preprocessed/preprocessing_metadata.npz

Shapes from exported arrays:

- epochs: (1, 21, 220000)
- epochs_rereferenced: (1, 21, 220000)
- epochs_with_buffer: (1, 21, 220600)

---

## 3. Preprocessing Pipeline

## 3.1 Step 1: Data loading and channel cleanup

The EDF recording is loaded with MNE, converted from volts to microvolts, and cleaned by dropping channels 25+, 26+, and 27+.

## 3.2 Step 2: Continuous filtering

Filtering is performed on continuous data before epoching to avoid edge artifacts contaminating short trial windows.

Filters:

1. High-pass Butterworth filter

- Cutoff: 0.5 Hz
- Order: 4
- Application: zero-phase filtering via filtfilt

1. Notch filter

- Center: 60 Hz
- Width: 2 Hz (Q = 30)
- Application: zero-phase filtering via filtfilt

Observed global effect:

- Mean voltage before filtering: 6.6077 uV
- Mean voltage after filtering: 0.0012 uV
- RMS before filtering: 37.8612 uV
- RMS after filtering: 31.8578 uV
- RMS reduction: 15.86%

## 3.3 Step 3: Epoching with buffer zones

Epoching is done after filtering.

Configured values:

- Trial duration: 1100.0 s
- Buffer duration: 1.5 s on each side
- Buffer rationale: preserve low-frequency validity and keep boundary effects away from trial core

## 3.4 Step 4: Common Average Reference (CAR)

CAR is applied on epoched data by subtracting the cross-channel mean at each time sample.

Per-channel mean RMS across trials:

- Before CAR: 0.005485 uV
- After CAR: 0.004203 uV

This confirms reduction of global common-mode offset.

## 3.5 Step 5: Export for downstream analysis

The final preprocessed arrays and metadata are saved in the preprocessed folder for reproducible downstream analysis.

---

## 4. DTFT and FFT Analysis

The Fourier notebook analyzes two signals:

1. Simulated EEG-like signal (1 s at 1000 Hz).
2. Real preprocessed EEG segment (1 s at 200 Hz).

### 4.1 Manual DTFT formulation

Manual DTFT/DFT is implemented as:

$$
X[k] = \sum_{n=0}^{N-1} x[n] e^{-j 2\pi kn/N}
$$

This is a dot product between the signal and complex sinusoidal basis functions. The DC term corresponds to k = 0.

### 4.2 FFT and IFFT validation

For each signal:

1. Compute manual DTFT coefficients.
2. Compute FFT with NumPy.
3. Compare maximum coefficient difference.
4. Reconstruct using IFFT and measure max reconstruction error.

### 4.3 One-sided amplitude scaling

A one-sided spectrum is computed from 0 Hz to Nyquist:

- Normalize by N.
- Double all non-DC, non-Nyquist bins.

This preserves physically interpretable amplitudes for real-valued signals.

---

## 5. Results

## 5.1 Simulated EEG-like signal (7, 15, 30 Hz + noise)

Configuration:

- fs = 1000 Hz
- duration = 1.0 s
- components: 7 Hz (2.5), 15 Hz (1.2), 30 Hz (0.8)

Numerical checks:

- Max |DTFT - FFT| = 1.885e-10
- Max IFFT reconstruction error = 2.220e-15

Top spectral peaks:

- 7.0 Hz, amplitude 2.4936
- 15.0 Hz, amplitude 1.1478
- 30.0 Hz, amplitude 0.7813

These recovered amplitudes align closely with the simulated ground truth.

## 5.2 Preprocessed real EEG segment (1 s)

Configuration:

- fs = 200 Hz
- segment length = 200 samples
- signal from trial 1, channel 1 (after mean removal)

Numerical checks:

- Max |DTFT - FFT| = 8.210e-11
- Max IFFT reconstruction error = 2.132e-14

Top spectral peaks:

- 1.0 Hz, amplitude 28.9870
- 3.0 Hz, amplitude 10.9058
- 4.0 Hz, amplitude 6.9667
- 2.0 Hz, amplitude 5.2188
- 5.0 Hz, amplitude 4.9251

Dominant low-frequency content is expected in slow cortical drift and residual low-frequency neural/non-neural components.

---

## 6. Stationarity Demonstration

A stationarity comparison was performed between:

1. A 1-second EEG window.
2. The full trial segment from the same channel.

If EEG were strictly stationary, their normalized spectra should remain similar.

Observed metrics:

- Relative L2 spectral difference: 33.812
- Low/high bandpower ratio (1 to 8 Hz) / (8 to 30 Hz), 1-second window: 110.227
- Low/high bandpower ratio (1 to 8 Hz) / (8 to 30 Hz), full trial: 44.740

Interpretation:

The large spectral mismatch indicates clear nonstationarity across the full trial duration. This supports using localized time-frequency methods (for example STFT or wavelets) for temporal dynamics.

---

## 7. Theoretical Mapping to Core Concepts

## 7.1 Euler representation and complex basis

The transform uses complex exponentials to capture both cosine and sine projections:

$$
e^{-j\theta} = \cos(\theta) - j\sin(\theta)
$$

## 7.2 Nyquist theorem

For sampling rate fs, the highest resolvable frequency is fs/2. Unique bins for real signals are N/2 + 1.

## 7.3 DC component

The k = 0 coefficient measures constant offset (signal mean contribution).

## 7.4 Convolution theorem relevance

FFT enables efficient filtering because time-domain convolution corresponds to frequency-domain multiplication.

---

## 8. Limitations and Notes

1. Manual DTFT is O(N^2), so practical only for short windows.
2. The real EEG Fourier demonstration uses one channel and one selected 1-second segment for tractable manual DTFT.
3. Spectral interpretation can change by channel, epoch selection, and referencing strategy.
4. Full stationarity conclusions should ideally be extended with sliding-window statistics across channels.

---

## 9. Reproducibility

To reproduce results:

1. Run all cells in preprocessing notebook to regenerate filtered, epoched, re-referenced exports.
2. Run all cells in FFT notebook for simulated-signal validation, real EEG Fourier analysis, and stationarity comparison.

Notebook order:

1. project/notebooks/preprocessing.ipynb
2. project/notebooks/FFT.ipynb

---

## 10. Conclusion

The preprocessing pipeline successfully reduced drift/noise and produced reusable analysis-ready EEG arrays. Manual DTFT and FFT produced numerically identical spectra up to floating-point precision, and IFFT reconstructed signals with near-zero error. The added 1-second versus full-trial comparison demonstrated strong nonstationarity in real EEG, motivating localized time-frequency methods for comprehensive neural dynamics analysis.
