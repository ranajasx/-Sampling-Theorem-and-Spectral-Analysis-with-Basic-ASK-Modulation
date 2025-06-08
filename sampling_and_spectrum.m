% sampling_and_spectrum.m
% Demonstrates Sampling Theorem and Spectral Analysis with Basic ASK.

clear;
close all;
clc;

% --- 1. System Parameters Definition ---
Fs_cont = 10000; % High sampling rate to approximate a continuous-time signal (Hz)
T_sim = 1;       % Simulation duration (seconds)
t_cont = 0:1/Fs_cont:T_sim-(1/Fs_cont); % Time vector for continuous signal

% --- Part A: Sampling Theorem Demonstration ---
fprintf('--- Part A: Sampling Theorem Demonstration ---\n');

F_signal = 5; % Frequency of sine wave (Hz)
analog_signal = cos(2*pi*F_signal*t_cont); % A simple analog signal

% Scenario 1: Sampling ABOVE Nyquist Rate
Fs_above_nyquist = 2 * F_signal * 2; % Example: 4x Nyquist rate
Ts_above = 1 / Fs_above_nyquist;
t_sample_above = 0:Ts_above:T_sim-(Ts_above);
sampled_signal_above = cos(2*pi*F_signal*t_sample_above);

% Scenario 2: Sampling BELOW Nyquist Rate (demonstrates aliasing)
F_alias = Fs_above_nyquist - F_signal; % Expected aliased frequency
Fs_below_nyquist = 2 * F_signal * 0.8; % Example: below Nyquist rate
Ts_below = 1 / Fs_below_nyquist;
t_sample_below = 0:Ts_below:T_sim-(Ts_below);
sampled_signal_below = cos(2*pi*F_signal*t_sample_below); % The same analog signal

% Plotting Sampling Theorem
figure;
subplot(3,1,1);
plot(t_cont, analog_signal);
title('Original Analog Signal (5 Hz)');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;

subplot(3,1,2);
plot(t_cont, analog_signal, 'b-', t_sample_above, sampled_signal_above, 'ro');
title(sprintf('Sampled Above Nyquist (Fs = %g Hz, Recoverable)', Fs_above_nyquist));
xlabel('Time (s)'); ylabel('Amplitude');
grid on;
legend('Analog', 'Samples');

subplot(3,1,3);
plot(t_cont, analog_signal, 'b-', t_sample_below, sampled_signal_below, 'ro');
hold on;
% Plot a lower frequency sine wave to show the aliased signal
plot(t_cont, cos(2*pi*F_alias*t_cont), 'g--'); % Shows the aliased component
title(sprintf('Sampled Below Nyquist (Fs = %g Hz, Aliasing to %g Hz)', Fs_below_nyquist, F_alias));
xlabel('Time (s)'); ylabel('Amplitude');
grid on;
legend('Analog', 'Samples', 'Aliased Signal');

sgtitle('Demonstration of Sampling Theorem and Aliasing'); % Super title for the figure

% --- Part B: Spectral Analysis with FFT ---
fprintf('\n--- Part B: Spectral Analysis with FFT ---\n');

Fs_fft = 100; % Sampling rate for FFT example (Hz)
T_fft = 1;    % Duration for FFT example
t_fft = 0:1/Fs_fft:T_fft-(1/Fs_fft);

F_freq1 = 10; % Component 1 frequency
F_freq2 = 25; % Component 2 frequency
signal_for_fft = sin(2*pi*F_freq1*t_fft) + 0.5*cos(2*pi*F_freq2*t_fft);

N_fft = length(signal_for_fft);
Y = fft(signal_for_fft); % Compute FFT
P2 = abs(Y/N_fft);      % Two-sided spectrum
P1 = P2(1:N_fft/2+1);   % One-sided spectrum
P1(2:end-1) = 2*P1(2:end-1); % Multiply by 2 for power in one side

f_axis = Fs_fft*(0:(N_fft/2))/N_fft; % Frequency axis

figure;
subplot(2,1,1);
plot(t_fft, signal_for_fft);
title(sprintf('Time Domain Signal (F1=%g Hz, F2=%g Hz)', F_freq1, F_freq2));
xlabel('Time (s)'); ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(f_axis, P1);
title('Single-Sided Amplitude Spectrum');
xlabel('Frequency (Hz)'); ylabel('|Y(f)|');
grid on;
sgtitle('FFT Spectrum Analysis');

% --- Part C: Basic Amplitude Shift Keying (ASK) and its Spectrum ---
fprintf('\n--- Part C: Basic ASK Modulation and Spectrum ---\n');

% Message signal (binary data)
num_bits_ask = 20;
bits_ask = randi([0 1], num_bits_ask, 1);

% Carrier parameters
Fc = 50; % Carrier frequency (Hz)
Fs_carrier = 500; % Sampling rate for carrier (must be > 2*Fc)
t_carrier_res = 0:1/Fs_carrier: (num_bits_ask/10)-(1/Fs_carrier); % Time vector for a few bits

% Create the ASK modulated signal
samples_per_bit = Fs_carrier / (Fc * 2); % roughly 10 samples per cycle if bit rate = Fc/2
% Let's fix samples_per_bit to get distinct pulses
bit_duration = 0.05; % Duration of one bit (s)
samples_per_bit_fixed = Fs_carrier * bit_duration; % Number of samples per bit interval

ask_signal = zeros(1, num_bits_ask * samples_per_bit_fixed);
t_ask = 0:1/Fs_carrier:(length(ask_signal)/Fs_carrier)-(1/Fs_carrier);

for i = 1:num_bits_ask
    bit_val = bits_ask(i);
    start_sample = (i-1)*samples_per_bit_fixed + 1;
    end_sample = i*samples_per_bit_fixed;

    carrier_segment = cos(2*pi*Fc*t_ask(start_sample:end_sample));

    if bit_val == 1 % Transmit carrier for '1'
        ask_signal(start_sample:end_sample) = carrier_segment;
    else % Transmit nothing (or a lower amplitude carrier) for '0'
        ask_signal(start_sample:end_sample) = zeros(size(carrier_segment)); % On-Off Keying (OOK)
    end
end

% Plotting ASK Signal
figure;
subplot(2,1,1);
plot(t_ask, ask_signal);
title(sprintf('Basic ASK (On-Off Keying) Signal (Carrier: %g Hz)', Fc));
xlabel('Time (s)'); ylabel('Amplitude');
grid on;
xlim([0 num_bits_ask/10]); % Show a few bits

% Compute and plot spectrum of ASK signal
N_ask_fft = length(ask_signal);
Y_ask = fft(ask_signal);
P2_ask = abs(Y_ask/N_ask_fft);
P1_ask = P2_ask(1:N_ask_fft/2+1);
P1_ask(2:end-1) = 2*P1_ask(2:end-1);
f_ask_axis = Fs_carrier*(0:(N_ask_fft/2))/N_ask_fft;

subplot(2,1,2);
plot(f_ask_axis, P1_ask);
title('ASK Signal Spectrum');
xlabel('Frequency (Hz)'); ylabel('|Y(f)|');
grid on;
xlim([0 2*Fc]); % Zoom in on relevant frequencies
sgtitle('Basic ASK Modulation and Spectrum');