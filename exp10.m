N = 1000;
E0 = 1;
fc = 2e9;
Cn = rand(1, N);
Cn = Cn / sum(Cn);
t = linspace(0, 1, 1000);
Tc = zeros(size(t));
Ts = zeros(size(t));
for n = 1:N
phase_n = rand * 2 * pi;
Tc = Tc + E0 * Cn(n) * cos(2 * pi * fc * t + phase_n);
Ts = Ts + E0 * Cn(n) * sin(2 * pi * fc * t + phase_n);
end
Ez_field = Tc .* cos(2 * pi * fc * t) - Ts .* sin(2 * pi * fc * t);
v = 30;
angle_of_arrival_deg = 30;
angle_of_arrival_rad = deg2rad(angle_of_arrival_deg);
c = 3e8;
doppler_shift = (v / c) * fc * cos(angle_of_arrival_rad);
N = 10^6;
x = randn(1, N);
y = randn(1, N);
z = (x + 1i * y);
zBin = [0:0.01:7];
sigma2 = 1;
pzTheory = (zBin / sigma2) .* exp(-(zBin.^2) / (2 * sigma2));
[nzSim, zBinSim] = hist(abs(z),zBin);
thetaBin = [-pi:0.01:pi];
pThetaTheory = 1 / (2 * pi) * ones(size(thetaBin));
[nThetaSim, thetaBinSim] = hist(angle(z), thetaBin);
figure;
subplot(2, 1, 1);
plot(zBinSim, nzSim / (N * 0.01), 'm', 'LineWidth', 2);
hold on;
plot(zBin, pzTheory, 'b.-');
xlabel('z');
ylabel('Probability Density, p(z)');
legend('Simulation', 'Theory');
title('Probability Density Function of |z|');
axis([0 7 0 0.7]);
grid on;
subplot(2, 1, 2);
plot(thetaBinSim, nThetaSim / (N *0.01), 'm', 'LineWidth', 2);
hold on;
plot(thetaBin, pThetaTheory, 'b.-');
xlabel('θ');
ylabel('Probability Density, p(θ)');
legend('Simulation', 'Theory');
title('Probability Density Function of θ');
axis([-pi pi 0 0.2]);
grid on;
Fc = 0;
Fm = doppler_shift;
fs = 1 / (t(2) - t(1));
f_axis = linspace(-fs / 2, fs / 2, length(t));
frequency_response = (1.5 / (pi * Fm)) * sqrt(1 - ((f_axis - Fc) / Fm).^2);
filtered_signal = ifft(fft(Ez_field).* fftshift(frequency_response));
figure;
plot(t, abs(filtered_signal));
xlabel('Time (s)');
ylabel('Received Signal Amplitude(r)');
title('Received Signal (Magnitude of Ez\_field)');
figure;
plot(t, doppler_shift * ones(size(t)), 'r--');
xlabel('Time (s)');
ylabel('Doppler Shift (Hz)');
title('Doppler Shift Over Time');
psd_filtered = (1 / (fs * length(filtered_signal))) * abs(fft(filtered_signal)).^2;
f_axis_filtered = linspace(-fs / 2, fs / 2, length(psd_filtered));
figure;
plot(f_axis_filtered, psd_filtered);
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectrum of Filtered Received Signal');
N = 1000;
E0 = 1;
fc = 2e9;
fs = 1000;
t = linspace(0, 1, fs);
num_waveforms = 3;
cross_corr_matrix = zeros(num_waveforms, num_waveforms);
rt_waveforms = cell(num_waveforms, 1);
for waveform_idx = 1:num_waveforms
Tc = zeros(size(t));
Ts = zeros(size(t));
for n = 1:N
phase_n = rand * 2 * pi;
Tc = Tc + E0 * rand * cos(2 * pi * fc * t + phase_n);
Ts = Ts + E0 * rand * sin(2 * pi * fc * t + phase_n);
end
Ez_field = Tc .* cos(2 * pi * fc * t) - Ts .* sin(2 * pi * fc * t);
Fc = 0;
Fm = doppler_shift;
f_axis = linspace(-fs / 2, fs / 2, length(t));
frequency_response = (1.5 / (pi * Fm)) * sqrt(1 - ((f_axis - Fc) / Fm).^2);
filtered_signal = ifft(fft(Ez_field) .* fftshift(frequency_response));
rt_waveforms{waveform_idx} = abs(filtered_signal);
end
for i = 1:num_waveforms
for j = 1:num_waveforms
cross_corr = xcorr(rt_waveforms{i}, rt_waveforms{j});
cross_corr_matrix(i, j) = max(cross_corr);
end
end
disp('Cross-Correlation Matrix:');
disp(cross_corr_matrix)
