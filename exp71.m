clc;
clear all;
close all;
%input sequence.
input_bits =input('Entera 4 bit input:');
t=0:0.01:28;
%PN Sequence
PN_Out=[];
s1o=1;
s2o=0;
s3o=0;
for i=1:7
 PN_Out=[PN_Out s3o];
 s3n =s2o;
 s2n=s1o;
 s1n=xor(s1o,s3o);
 slo=s1n;
 s2o=s2n;
 s3o=s3n;
end
%disp(PN_Out);
%PN to bipolar
for i=1:7
 if (PN_Out(i)==0)
 PN_Out(i)= -1;
 end
end
%disp[PN_Out);
PN_Outx= repmat(repelem(PN_Out,100), 1, 4);
figure(1);
subplot(3, 1, 1);
plot(t(1:2800), PN_Outx);
axis([0 28 -1.3 1.3]);
grid on;
box on;
xlabel('Time (in Tb)');
ylabel('Amplitude (in V)');
title('Pseudo Noise Sequence');
%input to bipolar and extension
input_ext=[];
for i=1:4
 if(input_bits(i)==0)
 input_ext = [input_ext repmat(-1,1,7)];
 else
 input_ext = [input_ext repmat(1,1,7)];
 end
end
%disp(input_ext);
input_extx = repelem(input_ext,100);
figure(1);
subplot(3, 1, 2);
plot(t(1:2800),input_extx);
axis([0 28 -1.3 1.3]);
grid on;
box on;
xlabel('Time (in Tb)');
ylabel('Amplitude (in V)');
title('Input Sequence');
%multiplying Input and PN.
final_out_bits = [];
for i=1:7:28
 for j=1:7
 final_out_bits = [final_out_bits input_ext(i+j-1)*PN_Out(j)];
 end
end
%disp(final_out_bits);
final_out_bitsx = repelem(final_out_bits,100);
figure(1);
subplot(3, 1, 3);
plot(t(1:2800), final_out_bitsx);
axis([0 28 -1.3 1.3]);
grid on;
box on;
xlabel('Time (in Tb)');
ylabel('Amplitude (in V)');
title('Input Sequence multiplied by PN Sequence');
t=0:0.01:56;
BPSK_in = repelem(final_out_bitsx,2);
figure(2);
subplot(3, 1, 1);
plot(t(1:5600),BPSK_in);
axis([0 56 -1.3 1.3]);
grid on;
box on;
xlabel('Time (in secs)');
ylabel('Amplitude (in V)');
title('Input bits');
sine_wave = sin(2*pi*t);
figure(2);
subplot(3, 1, 2);
plot(t,sine_wave);
axis([0 56 -1.3 1.3]);
grid on;
box on;
xlabel('Time (in secs)');
ylabel('Amplitude (in V)');
title('Carrier Wave');
DSSS_wave = sine_wave(1:5600).*BPSK_in;
figure(2);
subplot(3, 1, 3);
plot(t(1:5600),DSSS_wave);
axis([0 56 -1.3 1.3]);
grid on;
box on;
xlabel('Time (in secs)');
ylabel('Amplitude (in V)');
title('Direct Sequence Spread Spectrum Wave');