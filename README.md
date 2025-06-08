# Sampling-Theorem-and-Spectral-Analysis-with-Basic-ASK-Modulation

## Project Overview
This project explores fundamental digital signal processing (DSP) concepts that are crucial for understanding modern wireless communication systems. It includes:
-   A demonstration of the **Nyquist-Shannon Sampling Theorem**, illustrating the effects of sampling at various rates, including **aliasing**.
-   **Spectral analysis using the Fast Fourier Transform (FFT)** to visualize the frequency content of time-domain signals.
-   Implementation of a **basic Amplitude Shift Keying (ASK)** modulation scheme and analysis of its spectral characteristics.

All functionalities are implemented using **core MATLAB functions**, ensuring compatibility with basic MATLAB Online environments where specialized toolboxes are unavailable.

This project showcases:
-   **Core Digital Signal Processing (DSP) Principles:** Practical understanding of sampling, aliasing, FFT, and elementary modulation.
-   **MATLAB for Engineering Simulation:** Proficiency in using MATLAB for signal generation, manipulation, and visualization.
-   **Resourcefulness & Problem-Solving:** Ability to implement complex DSP concepts from first principles, highlighting adaptability to tool limitations.

## DSP & Communication Concepts Demonstrated
-   **Nyquist-Shannon Sampling Theorem:** Visualized continuous-time signals converting to discrete-time samples and the conditions for perfect reconstruction.
-   **Aliasing:** Explicitly showing how undersampling a signal causes higher frequencies to appear as lower frequencies (aliased components).
-   **Fast Fourier Transform (FFT):** Transforming signals from the time domain to the frequency domain to reveal their spectral content.
-   **Single-Sided Amplitude Spectrum:** Interpreting frequency domain plots to identify signal components.
-   **Basic Amplitude Shift Keying (ASK) Modulation:** A simple digital modulation technique, demonstrating how binary data can modulate a carrier signal.
-   **ASK Signal Spectrum:** Analyzing the frequency characteristics of an ASK modulated signal, illustrating bandwidth occupation.

4.  **Output:** The script will execute and automatically generate **three separate figure windows:** 
    * A figure demonstrating the **Sampling Theorem** with examples of sampling above and below the Nyquist rate.
    * A figure showing a time-domain signal and its corresponding **frequency spectrum** obtained via FFT.
    * A figure illustrating a **basic ASK modulated signal** and its spectrum.

## Expected Results
-   The importance of sampling rate for accurate signal representation and reconstruction.
-   How aliasing distorts the signal when the sampling rate is insufficient.
-   The distinct frequency components of signals clearly resolved by the FFT.
-   The spectral characteristics of a simple modulated signal (ASK), showing how modulation shifts and spreads the signal's energy in frequency.

