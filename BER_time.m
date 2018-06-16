clear all;
close all;

% generating gaussian pulse 
fs = 100e+9; %sampling frequency 100GHz
L = 1023; % Length of signal
sigma = 10e-10; % Pulse duration 1ns
variance=sigma^2;
NFFT = 1024;
T = 1/fs;
t = (-L/2:L/2)*T; % Time base
t0 = max(t)/2; % Used to centering the pulse
%Egauss1 = 1/(sqrt(2*pi*variance))*(exp(-2*log(2)*(t).^2/(sigma)^2));
Egauss = 1/(sqrt(2*pi*variance))*(exp(-2*log(2)*(t).^2/(sigma)^2));


% xlabel('Time (fs)');
% ylabel('Amplitude');
% grid on
%load('delay.mat');
% fft of gaussian pulse
freq = fs*linspace(-NFFT/2,NFFT/2,NFFT);
df = fs/NFFT;
dw = 2*pi*df;
dt = 1./fs;
t1 = 2047*dt;
t_axis = (-t1/2+dt:dt:t1/2);
load('Momega.mat');
load('convres.mat');
M_omega_t = ifft(M_omega,1024,3);
for i = 1:NFFT
    omega(i) = 2*pi*freq(i)/NFFT;
    delay(i) = exp(-1j*omega(i)*(5e-10));
    M_inv(:,:,i) = inv(M_omega(:,:,i));
end 
M_inv_t = ifft(M_inv,1024,3);
iter = 1;
total_bits = 6*iter;
num_bit = 6;
snr_c = [10];
for snr = 1:length(snr_c)
    err_bits = 0;
    for bit = 1:iter
        
        data=randint(1,num_bit);%random bit generation (1 or 0)
        s(bit,:)=2*data-1;%converssion of data for BPSK modulation
        %s = [1 1 1 -1 1 -1];
        for w = 1:6
            inp(w,:) = s(bit,w).*Egauss;
            figure(w);
            hold on;
            plot(t,inp(w,:)/norm(inp(w,:)),'b');
        end


        out1 = zeros(6,2047);
        for i = 1:6

            for j = 1:6
                for l = 1:NFFT
                    temp1(l) = M_omega_t(i,j,l);
                            
                end
                out1(i,:) = out1(i,:)+ conv(temp1,inp(j,:));
            end

        end
        
        sign_power_in = 0;
        for i = 1:NFFT
            sign_power_in = sign_power_in + out1(1,i)^2;
        end
        signp = sign_power_in;

        noiseEner = signp/(10^(snr_c(snr)/10));        % energy of noise to be added
        noiseVar = noiseEner/2;     % variance of noise to be added
        noiseStd = sqrt(noiseVar); 
        noise = noiseStd*(randn(6,2047)+1j*randn(6,2047));

        out1 = out1+noise;

        out2 = zeros(6,3070);

        for i = 1:6
            for j = 1:6
                for l = 1:NFFT
                    temp2(l) = M_inv_t(i,j,l);
                         
                end
                out2(i,:) = out2(i,:)+ conv(temp2,out1(j,:));
            end
            figure(i);
            hold on;
            
            plot(t_axis,(real(out1(i,:)))/norm(abs(out1(i,:))), 'r');
            legend('input','output');
            xlabel('Time(in ns)');
            ylabel('Magnitude');
        end

%



        for q =1:num_bit

                a = conv(Egauss,real(out2(q,:)));
                if real(a(2040)) < 0
                    b(bit,q) = -1;
                else
                    b(bit,q) = 1;
                end

                if b(bit,q) ~= s(bit,q)
                      err_bits = err_bits + 1;
                end
        end
end
   
    s;
    b;
    BER(snr) = err_bits/total_bits;
end
figure(8);
semilogy(snr_c,BER);
xlabel('SNR(in dB)');
ylabel('BER');