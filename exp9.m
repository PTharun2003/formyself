numTaps = 5; % Number of taps in the equalizer
channelDelay = 5; % Delay of the channel (for simulation purposes)
SNR = 20; % Signal-to-noise ratio (dB)
symbolRate = 10e3; 
sincDuration = 2; 
numSamples = symbolRate * sincDuration; 
t = linspace(-2, sincDuration, numSamples); 
sincSignal = sinc(10 * t);
channel = zeros(1, channelDelay + 1);
channel(channelDelay + 1) = 1; 
receivedSignal = filter(channel, 1, sincSignal); 
receivedSignal = awgn(receivedSignal, SNR, 'measured');
equalizerOutput = zeros(1, numSamples);
for i = numTaps + 1:numSamples
 equalizerOutput(i) = receivedSignal(i - numTaps:i) * sincSignal(i - numTaps:i).';
end
figure;
subplot(2, 1, 1);
plot(t, receivedSignal, 'b', 'LineWidth', 1.5);
title('Received Sinc Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
subplot(2, 1, 2);
plot(t, equalizerOutput, 'r', 'LineWidth', 1.5);
title('Signal after Tapped Delay Equalization');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;