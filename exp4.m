
bits = input('enter the sequence :');
bitrate = 1;
T = length(bits)/bitrate;
n = 200;
N = n*length(bits);
dt = T/N;
t = 0:dt:T;
x = zeros(1,length(t));
for i = 0:length(bits)-1
    if bits(i+1) == 1
        x(i*n+1:(i+1)*n) = 1;
    else
        x(i*n+1:(i+1)*n) = 0;
    end
end
figure(1)
plot(t,x,'linewidth',3);
grid on;
title('unipolar nrz');
xlabel('time')
ylabel('amplitude')
T = length(bits)/bitrate;
n = 200;
N = n*length(bits);
dt = T/N;
t = 0:dt:T;
x = zeros(1,length(t));
for i = 0:length(bits)-1
    if bits(i+1) == 1
        x(i*n+1:(i+1)*n) = 1;
    else
        x(i*n+1:(i+1)*n) = -1;
    end
end
figure(2);
plot(t,x,'linewidth',3);
title('polar nrx');
xlabel('time');
ylabel('amp');

T = length(bits)/bitrate;
n = 200;
N = n*length(bits);
dt = T/N;
t = 0:dt:T;
x = zeros(1,length(t));
for i = 0:length(bits)-1
    if bits(i+1) == 1
        x(i*n+1:(i+0.5)*n) = 1;
        x((i+0.5)*n+1:(i+1)*n) = 0;
    else
        x(i*n+1:(i+1)*n) = 0;
    end
end
figure(3);
plot(t,x,'linewidth',3);
title('unipolar rx');
xlabel('time');
ylabel('amp');
T = length(bits)/bitrate;
n = 200;
N = n*length(bits);
dt = T/N;
t = 0:dt:T;
x = zeros(1,length(t));
for i = 0:length(bits)-1
    if bits(i+1) == 1
        x(i*n+1:(i+0.5)*n) = 1;
        x((i+0.5)*n+1:(i+1)*n) = 0;
    else
        x(i*n+1:(i+0.5)*n) = -1;
        x((i+0.5)*n+1:(i+1)*n) = 0;
    end
end
figure(4);
plot(t,x,'linewidth',3);
title('polar rx');
xlabel('time');
ylabel('amp');
T = length(bits)/bitrate;
n = 200;
N = n*length(bits);
dt = T/N;
t = 0:dt:T;
x = zeros(1,length(t));
p = 1;
for i = 0:length(bits)-1
    if bits(i+1) == 1
        x(i*n+1:(i+1)*n) = p;
        p = -1*p;
    else
         x(i*n+1:(i+1)*n) = 0;
    end
end
figure(5);
plot(t,x,'linewidth',3);
title('bipolar nrx');
xlabel('time');
ylabel('amp');
T = length(bits)/bitrate;
n = 200;
N = n*length(bits);
dt = T/N;
t = 0:dt:T;
x = zeros(1,length(t));
p = 1;
for i = 0:length(bits)-1
    if bits(i+1) == 1
        x(i*n+1:(i+0.5)*n) = p;
        x((i+0.5)*n+1:(i+1)*n) = 0;
        p = -1*p;
    else
        x(i*n+1:(i+1)*n) = 0;
    end
end
figure(6);
plot(t,x,'linewidth',3);
title('bipolar rx');
xlabel('time');
ylabel('amp');
T = length(bits)/bitrate;
n = 200;
N = n*length(bits);
dt = T/N;
t = 0:dt:T;
x = zeros(1,length(t));
for i = 0:length(bits)-1
    if bits(i+1) == 1
        x(i*n+1:(i+0.5)*n) = 1;
        x((i+0.5)*n+1:(i+1)*n) = -1;
    else
        x(i*n+1:(i+0.5)*n) = -1;
        x((i+0.5)*n+1:(i+1)*n) = 1;
    end
end
figure(7);
plot(t,x,'linewidth',3);
title('man');
xlabel('time');
ylabel('amp');
Rb=1;
Tb=1/Rb;
f=0:0.05*Rb:2*Rb;
x=f*Tb;
P=Tb*(sinc(x).*sinc(x)); %Polar
P1=0.5*Tb*(sinc(x).*sinc(x))+ 0.5 *dirac(f); %Unipolar
P2=Tb*(sinc(x/2)).^2.*(sin(pi*x/2)).^2; %Manchester
P3=Tb*(sinc(x/2)).^2.*(sin(pi*x)).^2; %Bipolar
hold on
figure(9)
plot(f,P,'r')
hold on
plot(f,P1,'g')
plot(f,P2,'b')
plot(f,P3,'m')
grid on
box on
xlabel('f ---->')
ylabel('Power Spectral Density ---->')
title('PSD for Various Binary Line Codes')