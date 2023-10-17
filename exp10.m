% Parameters for the simulation
N = 1000; 
% Number ofwaves (adjust as needed)
E0 = 1; 
% Amplitudeof local average E-field (adjust asneeded)
fc = 2e9; 
% Carrier frequency in Hz (adjust as needed)
% Generate random amplitudes for Cn(Equation 4.58)
Cn = rand(1, N);
% Normalize amplitudes (Equation4.62)
Cn = Cn / sum(Cn);
% Initialize time parameters
t = linspace(0, 1, 1000); % Timevector (adjust as needed)
% Initialize arrays to store Tc(t)and Ts(t)
% Probability density function of
abs(z)
zBin = [0:0.01:7];
sigma2 = 1;
pzTheory = (zBin / sigma2) .* exp(-
(zBin.^2) / (2 * sigma2)); % Theory
[nzSim, zBinSim] = hist(abs(z),
zBin); % Simulation
% Probability density of theta
thetaBin = [-pi:0.01:pi];
pThetaTheory = 1 / (2 * pi) *
ones(size(thetaBin));
[nThetaSim, thetaBinSim] =
hist(angle(z), thetaBin); %
Simulation
% Plot the PDFs
figure;
subplot(2, 1, 1);
plot(zBinSim, nzSim / (N * 0.01),
'm', 'LineWidth', 2);
hold on;
plot(zBin, pzTheory, 'b.-');
xlabel('z');
ylabel('Probability Density, p(z)');
legend('Simulation', 'Theory');
title('Probability Density Function
of |z|');
axis([0 7 0 0.7]);
grid on;
subplot(2, 1, 2);
plot(thetaBinSim, nThetaSim / (N *
0.01), 'm', 'LineWidth', 2);
hold on;
plot(thetaBin, pThetaTheory, 'b.-');
xlabel('?');
ylabel('Probability Density, p(?)');
legend('Simulation', 'Theory');
title('Probability Density Function
of ?');
axis([-pi pi 0 0.2]);
grid on;
Tc = zeros(size(t));
Ts = zeros(size(t));
% Calculate Tc(t) and Ts(t) over time(Equations 4.64 and 4.65)
for n = 1:N
 % Random phase for the nthcomponent (Equation 4.61)
 phase_n = rand * 2 * pi;
 % Calculate Tc(t) and Ts(t)components
 Tc = Tc + E0 * Cn(n) * cos(2 * pi* fc * t + phase_n);
 Ts = Ts + E0 * Cn(n) * sin(2 * pi* fc * t + phase_n);
end
% Calculate E_z field component(Equation 4.63)
Ez_field = Tc .* cos(2 * pi * fc * t)- Ts .* sin(2 * pi * fc * t);
% Calculate the Doppler shift
v = 30; 
% Velocity of the mobilereceiver in m/s (adjust as needed)
angle_of_arrival_deg = 30; 
% Angle ofarrival in degrees (adjust as needed)
% Convert angle of arrival fromdegrees to radians
angle_of_arrival_rad =deg2rad(angle_of_arrival_deg);
% Calculate the Doppler shift in Hertz
% Apply the filter to the received
signal (Ez_field)
Fc = 0; % Center frequency (adjust as
needed)
Fm = doppler_shift; % Maximum Doppler
shift (adjust as needed)
% Define the frequency axis
fs = 1 / (t(2) - t(1)); % Sampling
frequency
f_axis = linspace(-fs / 2, fs / 2,
length(t));
% Calculate the filter response
frequency_response = (1.5 / (pi *
Fm)) * sqrt(1 - ((f_axis - Fc) /
Fm).^2);
% Apply the filter to the received
signal (Ez_field)
filtered_signal = ifft(fft(Ez_field)
.* fftshift(frequency_response));
% Plot the magnitude of Ez_field as
the received signal (r)
figure;
plot(t, abs(filtered_signal));
xlabel('Time (s)');
ylabel('Received Signal Amplitude
(r)');
title('Received Signal (Magnitude of
Ez\_field)');
% Plot the Doppler shift as before
figure;
plot(t, doppler_shift *
ones(size(t)), 'r--');
xlabel('Time (s)');
ylabel('Doppler Shift (Hz)');
title('Doppler Shift Over Time');
% Plot the power spectrum of the
filtered signal
psd_filtered = (1 / (fs *
length(filtered_signal))) *
abs(fft(filtered_signal)).^2;
f_axis_filtered = linspace(-fs / 2,
fs / 2, length(psd_filtered));
figure;
plot(f_axis_filtered, psd_filtered);
% Power spectrum
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectrum of Filtered
Received Signal');% Parameters for
the simulation
N = 1000; % Number of
waves (adjust as needed)
E0 = 1; % Amplitude
of local average E-field (adjust as
needed)
fc = 2e9; % Carrier
frequency in Hz (adjust as needed)
fs = 1000; % Sampling
frequency in Hz (adjust as needed)
% Initialize time parameters
t = linspace(0, 1, fs); % Time
vector (adjust as needed)
% Number of waveforms
num_waveforms = 3;
% Initialize cross-correlation matrix
cross_corr_matrix =
zeros(num_waveforms, num_waveforms);
% Generate three r(t) waveforms
rt_waveforms = cell(num_waveforms,
1);
for waveform_idx = 1:num_waveforms
 % Initialize arrays to store
Tc(t) and Ts(t)
 Tc = zeros(size(t));
 Ts = zeros(size(t));
 % Calculate Tc(t) and Ts(t) over
time (Equations 4.64 and 4.65)
 for n = 1:N
 % Random phase for the nth
component (Equation 4.61)
 phase_n = rand * 2 * pi;
 % Calculate Tc(t) and Ts(t)
components
 Tc = Tc + E0 * rand * cos(2 *
pi * fc * t + phase_n);
 Ts = Ts + E0 * rand * sin(2 *
pi * fc * t + phase_n);
 end
 % Calculate E_z field component
(Equation 4.63)
 Ez_field = Tc .* cos(2 * pi * fc
* t) - Ts .* sin(2 * pi * fc * t);
 
 % Apply the filter to the
received signal (Ez_field)
 Fc = 0; % Center frequency
(adjust as needed)
 Fm = doppler_shift; % Maximum
Doppler shift (adjust as needed)
 % Define the frequency axis
 f_axis = linspace(-fs / 2, fs /
2, length(t));
% Calculate the filter response
 frequency_response = (1.5 / (pi *
Fm)) * sqrt(1 - ((f_axis - Fc) /
Fm).^2);
 % Apply the filter to the
received signal (Ez_field)
 filtered_signal =
ifft(fft(Ez_field) .*
fftshift(frequency_response));
 % Store the generated r(t)
waveform
 rt_waveforms{waveform_idx} =
abs(filtered_signal);
end
% Calculate cross-correlation values
among r(t) waveforms
for i = 1:num_waveforms
 for j = 1:num_waveforms
 cross_corr =
xcorr(rt_waveforms{i},
rt_waveforms{j});
 cross_corr_matrix(i, j) =
max(cross_corr);
 end
end
% Display the cross-correlation
matrix
disp('Cross-Correlation Matrix:');
disp(cross_corr_matrix);