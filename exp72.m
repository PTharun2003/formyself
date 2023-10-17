clc;
close all;
clear all;
input_bits = input('Enter the input bits: ');
t=0:0.01:10;
t=t(1:end-1);
t1 =0:0.01:100;
t1=t1(1:end-1);
t2 =0:0.01:200;
t2 =t2(1:end-1);
%PN Sequence
PN_Out=[];
s1o = 1;
s2o = 0;
s3o = 0;
for i=1:7
 PN_Out=[PN_Out s3o];
  s3n=s2o;
 s2n=s1o;
 s1n= xor(s1o,s3o);
 slo=s1n;
 s2o=s2n;
 s3o= s3n;
end
disp(PN_Out);
PN_Outx = repmat(PN_Out, 1, 10);
disp(PN_Outx);
input_dec= [];
for i = 1:2:20
 if input_bits(i) == 0
   if input_bits(i+1) == 0
     input_dec = [input_dec 1];
   else
     input_dec = [input_dec 2];
   end
 else
 if input_bits(i+1) == 0
 input_dec= [input_dec 3];
 else
 input_dec= [input_dec 4];
 end
 end
end
PN_dec=[];
for i=1:3:60
 if PN_Outx(i)==0
if PN_Outx(i+1)==0
 if PN_Outx(i+2) == 0
 PN_dec=[PN_dec 0];
 else
 PN_dec=[PN_dec 1];
 end
else
 if PN_Outx(i+2) == 0
 PN_dec=[PN_dec 2];
 else
 PN_dec= [PN_dec 3];
 end
end
else
 if PN_Outx(i+1)==0
 if PN_Outx(i+2) ==0
 PN_dec= [PN_dec 4];
 else
 PN_dec=[PN_dec 5];
 end
 else
 if PN_Outx(i+2) == 0
 PN_dec= [PN_dec 6];
 else
 PN_dec= [PN_dec 7];
 end
 end
 end
end
input_dec_slow = input_dec;
sine_wave_input=[];
for i=1:10
 sine_temp = sin(2*pi^input_dec_slow(i)*t);
 sine_wave_input= [sine_wave_input sine_temp];
end
figure(1);
subplot(4, 1, 1);
plot(t1,sine_wave_input);
axis([0 100 -1.3 1.3]);
grid on;
box on;
xlabel('Time (in secs)');
ylabel('Amplitude (in V)');
title('Input after FSK Modulation');

figure(1)
subplot(4, 1, 2);
PN_map= repmat(repelem(PN_Out,660), 1, 50);
plot(t1,PN_map(1:10000));
xlabel('Time (in secs)');
ylabel('Amplitude (in V)');
title('PN sequence');

PN_dec_slow = repelem(PN_dec,2);
sine_wave_PN_slow =[];
for i=1:10;
 sine_temp = sin(2*pi*PN_dec_slow(i)*t);
 sine_wave_PN_slow = [sine_wave_PN_slow sine_temp];
end
figure(1);
subplot(4, 1, 3);
plot(t1,sine_wave_PN_slow);
axis([0 100 -1.3 1.3]);
grid on;
box on;
xlabel('Time (in secs)');
ylabel('Amplitude (in V)');
title('PN Sequence from Frequency Synthesizer');
Slow_freq_hop =[];
freq_slow = [];
for i=1:10
 freq = PN_dec_slow(i)*input_dec_slow(i);
 freq_slow = [freq_slow freq];
 sine_temp = sin(pi*freq*t);
 Slow_freq_hop = [Slow_freq_hop sine_temp];
end
figure(1);
subplot(4, 1,4);
plot(t1,Slow_freq_hop);
axis([0 100 -1.3 1.3]);
grid on;
box on;
xlabel('Time (in secs)');
ylabel('Amplitude (in V)');
title('Slow Frequency Hopping');
input_dec_fast = repelem(input_dec,2);
sine_wave_input =[];
for i=1:20
 sine_temp = sin(2*pi*input_dec_fast(i)*t);
 sine_wave_input=[sine_wave_input sine_temp];
end
figure(2);
subplot(4, 1, 1);
plot(t2,sine_wave_input);
axis([0 200 -1.3 1.3]);
grid on;
box on,
xlabel('Time (in secs)');
ylabel('Amplitude (in V)');
title('inputafter FSK Modulation');
figure(2)
subplot(4, 1, 2);
PN_map= repmat(repelem(PN_Out,333), 1, 100);
plot(t2,PN_map(1:20000));
xlabel('Time (in secs)');
ylabel('Amplitude (in V)');
title('PN sequence');
PN_dec_fast = PN_dec;
sine_wave_PN_fast = [];
for i=1:20
 sine_temp = sin(2*pi*PN_dec_fast(i)*t);
 sine_wave_PN_fast =[sine_wave_PN_fast sine_temp];
end
figure(2);
subplot(4, 1, 3);
plot(t2,sine_wave_PN_fast);
axis([0 200 -1.3 1.3]);
grid on;
box on;
xlabel('Time (in secs)');
ylabel('Amplitude (in V)');
title('PN Sequence from Frequency Synthesizer');
Fast_freq_hop = [];
freq_fast= [];
for i=1:20
 freq = PN_dec_fast(i)* input_dec_fast(i);
 freq_fast = [freq_fast freq];
 sine_temp = sin(pi*freq*t);
 Fast_freq_hop = [Fast_freq_hop sine_temp];
end
figure(2);
subplot(4, 1, 4);
plot(t2,Fast_freq_hop);
axis([0 200 -1.3 1.3]);
grid on;
box on;
xlabel('Time (in secs)');
ylabel('Amplitude (in V)');
title('Fast Frequency Hopping');