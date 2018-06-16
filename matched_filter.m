clear all;
close all;
load('output_1.mat');
load('gauss.mat');
load('Momega.mat');
load('delay.mat');
NFFT = 1024;
%count = 0;

for j = 1:6
    for k = 1:6
        for l = 1:NFFT
            temp(l) = X_out(j,k,l);
        end
        
        for n = 1:NFFT
            t(i) = conj(X_wrp(i))*delay(i);
        end
        save('mat_t','t');
        out = t.*temp;
        for m = 1:NFFT
            conv_res(j,k,m) = out(m);
        end
        
    end
end

for i = 1:NFFT
    %t(:,:) = eye(6);
    %t = conv_res(:,:,i)'*conv_res(:,:,i);
    conv_res(:,:,i) = conv_res(:,:,1)./(norm(conv_res(:,:,1)));
end

save('convres', 'conv_res');
        