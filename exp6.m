alpha= 1 + mod(282,3);
%%%Q1
fc = 10;
fs = 720;
Wp = 10/360;
Ws = 20/360;

[n,Wn] = buttord(Wp,Ws,3,40);

[b,a] = butter(n,Wn);
sys= tf(b,a,1/720)
figure
pzmap(sys)

figure
pzmap(sys)
grid on
title('Pole-Zero Plot');
figure
bode(sys)
title('Bode Plot');



% Define the time vector for 1 second
t = 0:1/fs:1;
impulse_response = impulse(sys, t);
step_response = step(sys, t);

figure;
plot(t, impulse_response, 'b', 'LineWidth', 2, 'DisplayName', 'Impulse Response');
hold on;
plot(t, step_response, 'r', 'LineWidth', 2, 'DisplayName', 'Step Response');
xlabel('Time (s)');
ylabel('Amplitude');
title('Impulse Response and Step Response');
legend('show');
grid on;

%%%% Q2

ecg_signal = load('ECG_Data.txt');
filtered_signal = filter(b,a,ecg_signal);

t = (0:length(ecg_signal)-1) / fs;

%Plot the original and filtered ECG signals
subplot(2,1,1);
plot(t, ecg_signal);
title('Original ECG Data');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t, filtered_signal);
title('Filtered ECG Data');
xlabel('Time (s)');
ylabel('Amplitude');

%%%%Q3
% Load the audio file
filename = 'instru3.wav';
[y, Fs] = audioread(filename);

% Compute the spectrogram of the audio signal
window_size = 1024;  % You can adjust this based on your preference
overlap = 512;       % Adjust overlap as needed
nfft = 2048;         % Adjust FFT size as needed

spectrogram(y, hamming(window_size), overlap, nfft, Fs, 'yaxis');
title('Original Spectrogram');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

% Design a digital Butterworth bandpass filter
f_passband = [100 1000];  % Passband frequency range (adjust as needed)
order = 7;                % Filter order (adjust as needed)

[b, a] = butter(order, f_passband / (Fs/2), 'bandpass');

% Apply the filter to the audio signal
filtered_signal = filter(b, a, y);

% Save the filtered audio to a new WAV file
output_filename = 'filtered_instru3.wav';
audiowrite(output_filename, filtered_signal, Fs);

% Plot the spectrogram of the filtered audio
figure;
spectrogram(filtered_signal, hamming(window_size), overlap, nfft, Fs, 'yaxis');
title('Filtered Spectrogram');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

% Play the filtered audio
sound(filtered_signal, Fs);
% 
% 
%%%%Q4
stopband_attenuation = 40; % Minimum stopband attenuation in dB
Fs = 720;            % Sampling frequency in samples/sec
Fp = 10;             % Passband edge frequency in Hz
Fstop = 20;          % Stopband edge frequency in Hz

% Calculate normalized cutoff frequencies (Wc and Ws)
Wc = 2 * Fp / Fs;
Ws = 2 * Fstop / Fs;

% Determine filter order (N) using the formula for Chebyshev Type I filter
[n, Wn] = cheb1ord(Wc, Ws, alpha, stopband_attenuation);

% Design the Chebyshev Type I filter
[z, p] = cheby1(n, alpha, Wn);
sos_chebyshev = tf(z, p, 1/Fs);

% Design the Butterworth filter as in the previous example
[n_butterworth, Wn_butterworth] = buttord(Wc, Ws, alpha, stopband_attenuation);
[z_butterworth, p_butterworth] = butter(n_butterworth, Wn_butterworth);
sos_butterworth = tf(z_butterworth, p_butterworth, 1/Fs);

% Compare the filter orders
fprintf('Chebyshev Type I Filter Order: %d\n', n);
fprintf('Butterworth Filter Order: %d\n', n_butterworth);

% Plot the Bode plots of Chebyshev Type I and Butterworth filters
figure;
bode(sos_chebyshev, 'r', sos_butterworth, 'b');
title('Bode Plot Comparison - Chebyshev Type I vs. Butterworth');
legend('Chebyshev Type I', 'Butterworth');
grid on;

% Define the time vector for 1 second
t = 0:1/Fs:1;

% Compute the impulse response of the Chebyshev Type I filter
impulse_response_chebyshev = impulse(sos_chebyshev, t);

% Compute the step response of the Chebyshev Type I filter
step_response_chebyshev = step(sos_chebyshev, t);

% Compute the impulse response and step response of the Butterworth filter
impulse_response_butterworth = impulse(sos_butterworth, t);
step_response_butterworth = step(sos_butterworth, t);

% Plot the impulse response and step response of both filters on the same graph
figure;
plot(t, impulse_response_chebyshev, 'r', 'LineWidth', 2, 'DisplayName', 'Chebyshev Type I');
hold on;
plot(t, impulse_response_butterworth, 'b', 'LineWidth', 2, 'DisplayName', 'Butterworth');
xlabel('Time (s)');
ylabel('Amplitude');
title('Impulse Response Comparison - Chebyshev Type I vs. Butterworth');
legend('Location', 'northeast');
grid on;

figure;
plot(t, step_response_chebyshev, 'r', 'LineWidth', 2, 'DisplayName', 'Chebyshev Type I');
hold on;
plot(t, step_response_butterworth, 'b', 'LineWidth', 2, 'DisplayName', 'Butterworth');
xlabel('Time (s)');
ylabel('Amplitude');
title('Step Response Comparison - Chebyshev Type I vs. Butterworth');
legend('Location', 'southeast');
grid on;
