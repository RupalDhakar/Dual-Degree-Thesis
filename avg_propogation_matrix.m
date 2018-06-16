clear all;
close all;

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
% figure(1);
% plot(t,Egauss/norm(Egauss),'b');
% xlabel('Time (fs)');
% ylabel('Amplitude');
% grid on
count = 10000;
snr_c = 21;
sign_power_in = 0;
for i = 1:L+1
    sign_power_in = sign_power_in + Egauss(i)^2;
end
signp = sign_power_in;
noiseEner = signp/(10^(snr_c/10));        % energy of noise to be added
noiseVar = noiseEner/2;     % variance of noise to be added
noiseStd = sqrt(noiseVar); 

% fft of gaussian pulse
freq = fs*linspace(-NFFT/2,NFFT/2,NFFT);
for i = 1:NFFT
    omega(i) = 2*pi*freq(i)/NFFT;
    delay(i) = exp(-1j*omega(i)*(5e-10));
end 
X = fft(Egauss,NFFT)/L;

% passing gaussian through fiber
%taking inverse fourier tranform of output

load('Momega.mat');
for i = 1:NFFT
    X_out(:,:,i) = X(i).*M_omega(:,:,i);
    %X_out_1(:,:,i) = awgn(X_out(:,:,i),10, 'measured');
end
noise1 = noiseStd*(randn(6,6,1024)+1j*randn(6,6,1024));
for z = 1:count
    
    for i = 1:NFFT
        %X_out_2(:,:,i,z) = awgn(X_out(:,:,i),0.1, 'measured');
        
        X_out_2(:,:,i,z) = X_out(:,:,i) + noise1(:,:,i) ;
    end
end

% 
for z = 1:count
    for j = 1:6
        for k = 1:6
            for l = 1:NFFT
                temp(l) = X_out_2(j,k,l,z);
            end
        
            for n = 1:NFFT
                t(n) = conj(X(n));
            end
        
            out = t.*temp;
        
            for m = 1:NFFT
                conv_res(j,k,m,z) = out(m);
            end
        end
    end
    
    for i = 1:NFFT
    %t(:,:) = eye(6);
    %t = conv_res(:,:,i)'*conv_res(:,:,i);
        conv_res(:,:,i,z) = conv_res(:,:,i,z)./(norm(conv_res(:,:,i,z)));
    end
end

for rup = 1:NFFT
    tempx = conv_res(:,:,rup,1);
    for z = 2:count
        tempx = tempx + conv_res(:,:,rup,z);
    end
    fin_conv_res(:,:,rup) = tempx./count;
end


save('convres', 'fin_conv_res');

%% using convulated propogation matrix to recover gaussian pulse
% 
% for p = 1:NFFT
%     rec_out_conv(:,:,p) = X_out_1(:,:,p)*inv(fin_conv_res(:,:,p));
% end
% 
% for o = 1:NFFT
%     rec_res_conv(o) = rec_out_conv(5,5,o);
% end
% rec_out_fin_conv = ((ifft(rec_res_conv)));
% 
% figure(1);
% hold on
% plot(t_axis,abs((rec_out_fin_conv))/norm(abs(rec_out_fin_conv)), 'k');
% %legend('input', 'pro_input', 'n_prop', 'n_prop_in', 'out', 'conv_out');
% hold off



