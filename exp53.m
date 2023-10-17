clc;
clear all;
close all;
b = 10;
bs = randi([0, 1], 1, b); % Generate random binary data
t = 0:0.001:b; % Time vector
snr_values = -10:2:100; % Range of SNR values to test
BER = zeros(size(snr_values)); %% Initialize BER array
for snr_idx = 1:length(snr_values)
 snr = snr_values(snr_idx);
 % Modulation parameters
 Ac = 1; % Carrier amplitude
 f = 1; % Carrier frequency
 fs = 1000; % Sampling frequency
 phi_0 = 0; % Phase for binary 0
 phi_1 = pi; % Phase for binary 1
 % Modulating signal
 mt = zeros(1, length(t));
 for i = 1:b
 mt((i-1)*fs+1:i*fs) = bs(i);
 end
 % PSK Modulation
 phi = zeros(1, length(t));
 for i = 1:b
 phi((i-1)*fs+1:i*fs) = bs(i) * phi_1 + (1 - bs(i)) *phi_0;
 end
 st = Ac * cos(2 * pi * f * t + phi);
 % AWGN channel
 snr = 10; % Signal-to-noise ratio
 received_signal = awgn(st, snr);
 % PSK Demodulation
 d = zeros(1, length(t));
 for i = 1:b
 x = sum(received_signal((i-1)*fs+1:i*fs) .* cos(2 * pi *f * t((i-1)*fs+1:i*fs) + phi_0)) / fs;
 if x < 0
 d((i-1)*fs+1:i*fs) = 1;
 else
 d((i-1)*fs+1:i*fs) = 0;
 end
 end
end
% Plotting
figure;
subplot(4,1,1);
plot(t, mt);
title('Modulating Signal');
xlabel('Time (sec)');
ylabel('Amplitude (volts)');
subplot(4,1,2);
plot(t, st);
title('PSK Modulated Signal');
xlabel('Time (sec)');
ylabel('Amplitude (volts)');
received_t = 0:0.001:length(received_signal) / fs - 0.001;
subplot(4,1,3);
plot(received_t, received_signal);
title('Received Signal (with AWGN)');
xlabel('Time (sec)');
ylabel('Amplitude (volts)');
subplot(4,1,4);
plot(t, d);
title('Demodulated Signal');
xlabel('Time (sec)');
ylabel('Amplitude');
num_bit = 1000;
BER_iter = 20;
Eb=1;
SNRdB=0:0.2:10;
SNR= 10.^(SNRdB/10);
for count = 1:length(SNR)
 avgError = 0;
 No = Eb/SNR(count);
 
 for run_time = 1:BER_iter
 Error=0;
 data = randi([0,1],1,num_bit);
 s=2*data-1;
 
 N=sqrt(No/2)*randn(1,num_bit);
 Y=s+N;
 
 for k=1:num_bit
 if ((Y(k)>0 && data(k)==0)||(Y(k)<0 && data(k)==1))
 Error=Error+1;
 end
 end
 Error=Error/num_bit;
 avgError=avgError+Error;
 end
 BER_sim(count)=avgError/BER_iter;
end
figure (3)
semilogy(SNRdB,BER_sim,'g', 'linewidth', 2.5);
grid on;
hold on;
BER_th = (1/2)*erfc(sqrt(SNR));
semilogy(SNRdB,BER_th,'r','linewidth',2.5);
grid on;
hold on;
title('Curve for Bit Error Rate verses SNR for PSK modulation');
xlabel('SNR(dB)');
ylabel('BER');
legend('Simulation','Theoretical');
