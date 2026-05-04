% DSP + EEG: Complete SOP (Theory → Implementation)

Purpose
-------

This document is a single, self-contained specification that explains from first principles why each DSP step is required for EEG analysis, what each formula means, and how to implement transforms so their outputs are physically interpretable. It is written as a theory-to-implementation contract rather than only code examples.

Environment
-----------

Use the project's `environment.yml` to create a reproducible Python environment. Typical required packages: `numpy`, `scipy`, `matplotlib`, `mne` (optional), `pandas`.

Big picture
-----------

We want a system that answers: "What frequencies exist in a noisy biological signal, how strong are they, and how do they evolve over time?"

Pipeline (summary)

RAW EEG
  ↓
DFT (structure inspection: DC, noise, anomalies)
  ↓
PREPROCESSING (DC removal, filtering)
  ↓
WINDOWING (leakage correction)
  ↓
STFT (time-local frequency analysis)
  ↓
PSD (energy normalization)
  ↓
BAND ANALYSIS (alpha, beta, gamma)
  ↓
BIOLOGICAL INTERPRETATION

1. Signal model (foundation)

--------------------------------
Model discrete-time EEG as
$$x[n] = s[n] + a[n] + \eta[n]$$
where:

- $s[n]$ — neural oscillations (true brain activity)
- $a[n]$ — artifacts (muscle, eye, powerline)
- $\eta[n]$ — noise (thermal, instrumentation)

This decomposition is conceptual: analysis algorithms must separate these components using frequency/time structure and heuristics.

1. Why frequency-domain exists

--------------------------------
Time-domain waveforms hide overlapping oscillations and can obscure rhythmic structure. The Fourier approach projects the signal onto orthogonal complex exponentials
$$e^{-j2\pi kn/N}$$
so that the decomposition isolates energy per frequency index $k$. The basis functions are orthogonal on the finite-length discrete-time vector space, letting us treat components independently (subject to windowing caveats below).

1. From continuous Fourier transform → DFT (derivation sketch)

------------------------------------------------------------
Start from continuous-time Fourier transform (CTFT) for band-limited signals and then sample at $f_s$ to obtain a discrete-time sequence. Truncating a finite time interval and evaluating its discrete Fourier series gives the DFT. The DFT is therefore the sampled Fourier transform of a finite-length signal segment, with implicit periodic extension.

1. Discrete Fourier Transform (DFT)

---------------------------------
Definition (length-$N$ DFT):
$$X[k] = \sum_{n=0}^{N-1} x[n]\,e^{-j\frac{2\pi}{N}kn}\qquad k=0,\dots,N-1$$

Meaning:

- $X[k]$ measures complex amplitude of the $k$-th Fourier basis over the entire window.
- $k=0$ (DC) equals the windowed-sum: $X[0]=\sum_n x[n]$ and represents offset/drift.

Important assumptions and consequences:

- Implicit periodic extension of the $N$-sample segment — if the segment does not align with an integer number of periods of underlying sinusoids, spectral leakage occurs.

1. The core problem: truncation and leakage

-------------------------------------------
Truncation multiplies the infinite-time signal by a rectangle in time. In frequency this is convolution with a sinc (the DFT of a rectangular window), which spreads energy across frequencies. This is spectral leakage and explains why a pure tone sampled into a non-harmonic-length window results in smeared energy.

1. Windowing as correction

--------------------------
Window the segment by $w[n]$ before computing the DFT:
$$x_w[n]=x[n]w[n],\\
X_w[k]=\sum_{n=0}^{N-1} x[n]w[n]e^{-j\frac{2\pi}{N}kn}$$

Time-multiplication → frequency-convolution. A smooth $w[n]$ has a frequency response with less sidelobe energy than a rectangle, thus suppressing leakage at the cost of widening the main lobe (poorer frequency resolution).

6. Window energy and normalization
----------------------------------
Windowing changes the signal's total energy. Define window energy
$$U=\sum_{n=0}^{N-1} w[n]^2.$$ The periodogram estimate of power spectral density (single segment) must be normalized by both sampling frequency $f_s$ and $U$ to obtain physical units (power per Hz):
$$\mathrm{PSD}(f_k) \approx \frac{|X_w[k]|^2}{f_s\,U}$$
where $f_k = k\frac{f_s}{N}$.

Rationale: $|X_w[k]|^2$ scales with squared amplitude and with the squared window amplitude. Dividing by $U$ returns an estimate proportional to the true signal power in physical units.

7. Common windows and trade-offs
--------------------------------
- Rectangular: $w[n]=1$. Narrowest main lobe, highest sidelobes (maximum leakage).
- Hann (Hanning):
$$w[n]=0.5\left(1-\cos\frac{2\pi n}{N-1}\right).$$
  Smooth edges, reduced sidelobes, modestly wider main lobe.
- Gaussian:
$$w[n]=\exp\left(-\frac{1}{2}\left(\frac{n-\mu}{\sigma}\right)^2\right)$$
  Very low sidelobes, wide main lobe — best leakage suppression, worst frequency sharpness.

Trade-off principle: reducing leakage increases main-lobe width, decreasing frequency resolution. Choose a window based on whether leakage suppression or sharp frequency discrimination is more important.

8. Frequency resolution and main-lobe width
------------------------------------------
Theoretical frequency resolution is proportional to main-lobe width in Hz. For length $N$ and sampling $f_s$, the rectangular main-lobe width (zero-to-zero) is approximately $2f_s/N$. For other windows multiply by a window-specific factor.

9. STFT: time-local spectra
---------------------------
DFT assumes stationarity across the whole segment. For non-stationary signals we use a sliding window with hop $H$:
$$X(m,k)=\sum_{n=0}^{N-1} x[n+mH]w[n]e^{-j\frac{2\pi}{N}kn}$$
where $m$ indexes time frames. STFT yields a time-frequency representation; resolution in time is proportional to $H$ and $N$, and in frequency to main-lobe width.

Overlapping windows
-------------------
To reduce variance and avoid missing transient energy, frames usually overlap (common choices: 50%–75% overlap). Overlap plus an appropriate synthesis window enables perfect or near-perfect reconstruction when needed.

10. PSD and averaging (Welch)
----------------------------
Single-segment periodogram is noisy. Welch's method averages multiple windowed periodograms (possibly overlapped) to reduce variance. For each segment $r$ compute $P_r(f_k)=|X_{w,r}[k]|^2/(f_s U)$ and average:
$$\mathrm{PSD}_{\mathrm{Welch}}(f_k)=\frac{1}{R}\sum_{r=1}^R P_r(f_k).$$

This averaging reduces variance (noise) at the cost of frequency resolution depending on segment length.

11. From PSD to band power (EEG interpretation)
----------------------------------------------
Band power in frequency band $B$ is a Riemann-sum over discrete frequency bins:
$$P_{\mathrm{band}} = \sum_{k:\,f_k\in B} \mathrm{PSD}(f_k)\Delta f\\ = \frac{f_s}{N}\sum_{k:\,f_k\in B} \mathrm{PSD}(f_k)$$
where $\Delta f=f_s/N$ is the bin width. Band definitions (typical): delta (0.5–4 Hz), theta (4–8 Hz), alpha (8–13 Hz), beta (13–30 Hz), gamma (>30 Hz).

12. Artifact decision rules (operational)
---------------------------------------
Any strong narrowband peak should be classified as one of: neural, muscular (EMG), environmental (50/60 Hz), or unknown. Operational heuristics:
- 50/60 Hz narrowband + harmonics → environmental
- Broad high-frequency increase (>30 Hz) across channels → EMG
- Channel-specific transient high amplitude → electrode artifact

13. Self-check questions (use as tests and unit checks)
-----------------------------------------------------
1. What distortion does a rectangular window introduce? (sinc sidelobes → leakage)
2. Why does windowing change amplitude but not true signal energy after normalization? (window multiplies samples, changing raw sum-of-squares; dividing by $U$ restores scale)
3. Why does STFT require overlapping windows? (variance reduction, transient capture, and reconstructability)
4. Why is PSD necessary even if DFT exists? (PSD expresses energy per Hz in physical units, comparable across sampling rates and window choices after normalization)
5. Why is EEG bandpower more meaningful than single frequency peaks? (physiological rhythms are broad and distributed; single-bin peaks may be artifacts)

14. Implementation notes (practical contract)
-----------------------------------------
- Always store sampling rate $f_s$ with recorded data.
- Remove DC (high-pass or detrend) before spectral estimates if drift exists.
- Choose window length $N$ considering stationarity: longer $N$ → better frequency precision, worse time localization.
- Use Welch averaging for stable PSD estimates; choose overlap for smoother time continuity.
- Normalize periodograms by $f_s$ and $U$ to obtain PSD in power/Hz.

15. Quick reference formulas
---------------------------
- DFT: $X[k]=\sum_{n=0}^{N-1}x[n]e^{-j2\pi kn/N}$
- STFT frame: $X(m,k)=\sum_{n=0}^{N-1}x[n+mH]w[n]e^{-j2\pi kn/N}$
- Window energy: $U=\sum_{n=0}^{N-1}w[n]^2$
- PSD (single frame): $\mathrm{PSD}(f_k)=\dfrac{|X_w[k]|^2}{f_s\,U}$

References & further reading
---------------------------
- Oppenheim & Schafer — Discrete-Time Signal Processing
- Welch, P.D. (1967) The use of fast Fourier transform for the estimation of power spectra: A method based on time averaging over short, modified periodograms.
- Brigham — The Fast Fourier Transform and Its Applications

Appendix: Example implementation checklist
-----------------------------------------
1. Load raw data, record `f_s`.
2. Detrend or high-pass to remove DC if needed.
3. Choose `N`, window `w[n]`, and hop `H`.
4. Compute STFT using chosen window and hop.
5. For each frame compute PSD via $|X_w|^2/(f_s U)$.
6. Average frames (Welch) or compute time-frequency visualization.
7. Compute band powers by summing PSD bins in defined bands.
8. Apply artifact rules and annotate.

-- end
