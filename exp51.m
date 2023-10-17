clc
clear all;
close all;
b = 10;
bs = randi([0, 1], 1, b); % Generate random binary data
t = 0:0.001:b; % Time vector
snr_values = -10:1:10000; % Range of SNR values to test
BER = zeros(size(snr_values)); %% Initialize BER array
%Modulation parameters
Ac = 1; % Carrier amplitude
fe = 2; % Carrier frequency
fs = 1000; % Sampling frequency
for snr_idx = 1:length(snr_values)
 snr = snr_values(snr_idx);
 %Modulating signal
 mt = zeros(1, length(t));
 for i = 1:b
 mt((i-1)*fs+1:i*fs) = bs(i);
 end
 %ASK Modulation
 ct = Ac * cos(2 * pi * fe * t);
 st = mt .* ct;
 %AWGN channel

 t = 0:0.001:b;
 received_signal = awgn(st, snr);
 %ASK Demodulation
 d = zeros(1, length(t));
 for i = 1:b
 x = sum(received_signal((i-1)*fs+1:i*fs).*ct((i-1)*fs+1:i*fs));
 if (x / fs) > 0
 d((i-1)*fs+1:i*fs) = 1;
 else
 d((i-1)*fs+1:i*fs) = 0;
 end
 end
end
%Plotting
figure;
subplot(4,1,1);
plot(t, mt);
title('Modulating Signal');
xlabel('Time (sec)');
ylabel('Amplitude (volts)');
subplot(4,1,2);
plot(t, st);
title('ASK Modulated Signal');
xlabel('Time (sec)');
ylabel('Amplitude (volts)');
subplot(4,1,3);
plot(t, received_signal);
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
for count=1:length(SNR)
 avgError=0;
 No=Eb/SNR(count);
 for run_time=1:BER_iter
 Error=0;
 data = randi([0 1],1,num_bit);
 Y= awgn(complex(data),SNRdB(count));
 for k=1:num_bit
 if ((Y(k)>0.5 && data(k)==0)||(Y(k)<0.5 && data(k)==1))
 Error=Error+1;
 end
 end
 Error=Error/num_bit;
 avgError=avgError+Error;
 end
 BER_sim(count)=avgError/BER_iter;
end
figure (1)
semilogy(SNRdB,BER_sim,'g','linewidth', 2.5);
grid on;
hold on;
BER_th = (1/2)*erfc(0.5*sqrt(SNR));
semilogy(SNRdB,BER_th,'r','linewidth',2.5);
grid on;
hold on;
title('Curve for Bit Error Rate verses SNR for ASK modulation');
xlabel('SNR(dB)');
ylabel('BER');
legend('Simulation','Theoretical');