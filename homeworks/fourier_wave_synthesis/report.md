# Waveform Construction Using Fourier Series – Report

**Course:** Biomedical Signal Processing  
**Topic:** Fourier Series Waveform Synthesis

---

## 1. Introduction

This report documents the reconstruction of three periodic waveforms — square, triangle, and sawtooth — using their Fourier series representations. The accuracy of the reconstruction is evaluated via Root Mean Square (RMS) error as the number of harmonics increases from 1 to 400.

---

## 2. Signal Parameters

| Parameter              | Value   |
|------------------------|---------|
| Sampling frequency     | 1000 Hz |
| Duration               | 1 s     |
| Fundamental frequency  | 5 Hz    |
| Amplitude              | 1       |
| Maximum harmonics      | 400     |

---

## 3. Fourier Series Formulas

### 3.1 Square Wave

$$
x(t) = \frac{4A}{\pi} \sum_{\substack{k=1,3,5,\ldots}}^{N} \frac{1}{k} \sin(2\pi k f_0 t)
$$

Uses **odd harmonics** only; coefficients decay as $1/k$.

### 3.2 Triangle Wave

$$
x(t) = \frac{8A}{\pi^2} \sum_{\substack{k=1,3,5,\ldots}}^{N} \frac{(-1)^{(k-1)/2}}{k^2} \sin(2\pi k f_0 t)
$$

Uses **odd harmonics** only; coefficients decay as $1/k^2$.

### 3.3 Sawtooth Wave

$$
x(t) = \frac{2A}{\pi} \sum_{k=1}^{N} \frac{(-1)^{k+1}}{k} \sin(2\pi k f_0 t)
$$

Uses **all harmonics**; coefficients decay as $1/k$.

---

## 4. Waveform Reconstruction

Each waveform was reconstructed at N = 5, 20, and 100 harmonics and visually compared against the SciPy reference signal.

| Figure | File |
|--------|------|
| Square wave reconstruction  | `figures/square_reconstruction.png`  |
| Triangle wave reconstruction | `figures/triangle_reconstruction.png` |
| Sawtooth wave reconstruction | `figures/sawtooth_reconstruction.png` |

At N = 5 harmonics, the general shape is recognizable but significant deviations exist. By N = 100, the reconstruction closely matches the reference for all waveforms, though overshoot persists at discontinuities for the square and sawtooth waves.

---

## 5. RMS Error Analysis

RMS error was computed for each harmonic count from 1 to 400:

$$
\text{RMS} = \sqrt{\frac{1}{N_s} \sum_{i=1}^{N_s} \left( x_{\text{ref}}[i] - x_{\text{recon}}[i] \right)^2 }
$$

Results are stored in:

| File | Waveform |
|------|----------|
| `data/square_rms.csv`   | Square   |
| `data/triangle_rms.csv` | Triangle |
| `data/sawtooth_rms.csv` | Sawtooth |

### Convergence Summary

| Waveform | Coefficient Decay | Convergence Speed |
|----------|------------------|-------------------|
| Square   | 1/k              | Slow              |
| Sawtooth | 1/k              | Slow              |
| Triangle | 1/k²             | Fast              |

The triangle wave reaches negligible RMS error with far fewer harmonics than the square or sawtooth waves. This is directly explained by the faster spectral roll-off (1/k² vs 1/k).

---

## 6. Gibbs Phenomenon

The **Gibbs phenomenon** manifests as persistent overshoots (≈ 9 % of the step height) near the jump discontinuities of the square and sawtooth waves. These overshoots do not vanish as the number of harmonics increases — they narrow but maintain the same amplitude. This is a fundamental property of Fourier series convergence at points of discontinuity.

The triangle wave, being continuous, does **not** exhibit Gibbs phenomenon, and its Fourier reconstruction converges uniformly.

---

## 7. Why Triangle Waves Converge Faster

The triangle wave is a **continuous** function — it has no jump discontinuities. Its first derivative is discontinuous (piecewise constant), but the function itself is smooth enough for the Fourier coefficients to decay as $1/k^2$.

In contrast, both the square and sawtooth waves have **jump discontinuities**, which limit the Fourier coefficient decay to $1/k$. A slower decay means more harmonics are needed for the partial sums to approximate the function, resulting in higher RMS error at any given harmonic count.

In general, smoother signals have faster-decaying Fourier coefficients and therefore converge more rapidly in the Fourier series sense.

---

## 8. Conclusion

This assignment confirmed that:

- Fourier series can reconstruct periodic waveforms with increasing accuracy as harmonics are added.
- RMS error monotonically decreases with harmonic count for all three waveforms.
- The rate of convergence depends on the smoothness of the signal (spectral decay rate).
- Gibbs phenomenon is an inherent limitation of Fourier truncation at discontinuities.

---

## 9. Deliverables

| Deliverable              | Location                                  |
|--------------------------|-------------------------------------------|
| Jupyter Notebook         | `notebook/wave_synthesis.ipynb`           |
| Square RMS CSV           | `data/square_rms.csv`                     |
| Triangle RMS CSV         | `data/triangle_rms.csv`                   |
| Sawtooth RMS CSV         | `data/sawtooth_rms.csv`                   |
| Square reconstruction    | `figures/square_reconstruction.png`       |
| Triangle reconstruction  | `figures/triangle_reconstruction.png`     |
| Sawtooth reconstruction  | `figures/sawtooth_reconstruction.png`     |
| RMS plot – Square        | `figures/rms_square.png`                  |
| RMS plot – Triangle      | `figures/rms_triangle.png`                |
| RMS plot – Sawtooth      | `figures/rms_sawtooth.png`                |
| Report                   | `report.md`                               |
