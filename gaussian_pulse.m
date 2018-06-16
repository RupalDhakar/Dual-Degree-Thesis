clear all;
close all;
%load('Momega.mat'); %propogation matrix

% section 1
% generating gaussian pulse 
fs = 100e+9; %sampling frequency 100GHz
L = 500; % Length of signal
sigma = 10e-10; % Pulse duration 1ns
variance=sigma^2;
NFFT = 1024;
T = 1/fs;
t = (-L/2:L/2)*T; % Time base
t0 = max(t)/2; % Used to centering the pulse
Egauss = 1/(sqrt(2*pi*variance))*(exp(-2*log(2)*(t).^2/(sigma)^2));
figure(1);
plot(t,Egauss/norm(Egauss),'b');
xlabel('Time (fs)');
ylabel('Amplitude');
grid on


% section 2
% fft of gaussian pulse
freq = fs*linspace(-NFFT/2,NFFT/2,NFFT);
for i = 1:NFFT
    omega(i) = 2*pi*freq(i)/NFFT;
    delay(i) = exp(-1j*omega(i)*(1e-10));
end 

save('delay', 'delay');

X = (fft(Egauss,NFFT));

for i = 1:NFFT/2
    X_wrp(NFFT/2+i) = (X(i));
end
for i = NFFT/2+1:NFFT
    X_wrp(i-NFFT/2) = (X(i));
end


figure(2);
plot(freq,abs(X_wrp)/norm(abs(X_wrp)));
title('Magnitude of FFT');
xlabel('Frequency (THz)')
ylabel('Magnitude |X(f)|');
grid on

sign_power_in = 0;
for i = 1:L
    sign_power_in = sign_power_in + Egauss(i)^2;
end
signp = sign_power_in;

noiseEner = signp/(10^(0/10));        % energy of noise to be added
noiseVar = noiseEner/2;     % variance of noise to be added
noiseStd = sqrt(noiseVar); 
noise = noiseStd*(randn(1,1024)+1j*randn(1,1024));
out = noise + X_wrp;


figure(2);
hold on;
plot(freq,(out)/norm(abs(out)));
title('Magnitude of FFT');
xlabel('Frequency (GHz)')
ylabel('Magnitude |X(f)|');
grid on
% df = fs/NFFT;
% dw = 2*pi*df;
% dt = 1./fs;
% t1 = 501*dt;
% figure(1);
% hold on
% t_axis = (-t1/2+dt:dt:t1/2);
% plot(t_axis,((out))/norm(abs(out)), 'r');
% % freq_1 = linspace(-NFFT/2,NFFT/2,NFFT);
% % plot(freq_1, abs(X_wrp)/norm(abs(X_wrp)));