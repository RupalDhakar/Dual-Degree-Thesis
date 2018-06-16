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
sign_power_in = 0;
for i = 1:L
    sign_power_in = sign_power_in + Egauss(i)^2;
end
signp = sign_power_in;


iter = 1;
total_bits = 6*iter;

snr_c = [0];
%%
for snr_bit = 1:length(snr_c)
    err_bits = 0;
    for bit = 1:iter

        num_bit=6;%number of bit
        data=randint(1,num_bit);%random bit generation (1 or 0)
        s(bit,:)=2*data-1;%conversion of data for BPSK modulation
        noiseEner = signp/(10^(snr_c(snr_bit)/10));        % energy of noise to be added
        noiseVar = noiseEner/2;     % variance of noise to be added
        noiseStd = sqrt(noiseVar); 
        noise = noiseStd*(randn(6,NFFT)+1j*randn(6,NFFT)); 
        noise = noise';


        for i = 1:num_bit
            inp(i,:) = s(bit,i).*Egauss;
            temp_plot = inp(i,:);
%             figure(i);
%             plot(t,temp_plot/norm(temp_plot),'b');
%             hold on
        end

        % fft of gaussian pulse
        freq = fs*linspace(-NFFT/2,NFFT/2,NFFT);
        for i = 1:NFFT
            omega(i) = 2*pi*freq(i)/NFFT;
            delay(i) = exp(-1j*omega(i)*(5e-10));
        end 
        for i = 1:num_bit
            X(i,:) = fft(inp(i,:),NFFT);
            for fu = 1:NFFT/2
                X_wrp(i,NFFT/2+fu) = (X(i,fu));
            end
            for fu = NFFT/2+1:NFFT
                X_wrp(i,fu-NFFT/2) = (X(i,fu));
            end
            figure(i);
            plot(freq,abs(X_wrp(i,:))/norm(X_wrp(i,:)),'b');
            hold on
        end


        % passing gaussian through fiber
        %taking inverse fourier tranform of output

        load('Momega.mat');
        

        %for r = 1:length(SNR_c)
            for i = 1:NFFT
                %for z = 1:num_bit
                X_out(i,:) = X(:,i)'*M_omega(:,:,i);
                %end
                %X_out1(i,:) = awgn(X_out(i,:),-20);
                X_out1(i,:) = X_out(i,:)+noise(i,:);


            end

            df = fs/NFFT;
            dw = 2*pi*df;
            dt = 1./fs;
            t1 = NFFT*dt;
            t_axis = (-t1/2+dt:dt:t1/2);

            % using convulated propogation matrix to recover gaussian pulse
            load('convres.mat');

            for p = 1:NFFT
                %for q = 1:6
                rec_out_conv1(p,:) = X_out1(p,:)*inv(M_omega(:,:,p));
                rec_out_conv = rec_out_conv1';
                %end
            end

            for q = 1:num_bit 

                for o = 1:NFFT
                    rec_res_conv(o) = rec_out_conv(q,o);
                end
                %rec_res_conv = awgn(rec_res_conv,10);
                rec_out_fin_conv = (ifft(rec_res_conv));
                for rupal = 1:501
                    final(rupal) = rec_out_fin_conv(rupal);
                end
                a = conv(final, Egauss);
                if real(a(501)) < 0
                    b(bit,q) = -1 ;
                else
                    b(bit,q) = +1;
                end

                if b(bit,q) ~= s(bit,q)
                    err_bits = err_bits + 1;
                end


%                 figure(q);
%                 hold on
%                 plot(t_axis,((rec_out_fin_conv))/norm(abs(rec_out_fin_conv)), 'k');
%                 legend('input', 'pro_input', 'n_prop', 'n_prop_in', 'out', 'conv_out');
%                 hold off
            end
    end
    BER(snr_bit) = err_bits/total_bits
end

semilogy(snr_c,BER);
xlabel('SNR (in dB)');
ylabel('BER');
