clear all;
close all;
len = 10000;
corr_len = 1;
Nseg = len/corr_len;
n = 6;
fs = 100e+8;
NFFT = 1024;
T = 1/fs;
freq = fs*linspace(-NFFT/2,NFFT/2,NFFT);

% omega vector
for i = 1:NFFT
    omega(i) = 2*pi*freq(i)/NFFT;
end 
     
%% matrix ofr Nseg segments , input and output coupling using haar measure

for i = 1:Nseg
        X = (randn(n) + i*randn(n))/sqrt(2);
        [Q1,R1] = qr(X);
        R1 = diag(diag(R1)./abs(diag(R1)));
        U = Q1*R1;
        Y = (randn(n) + i*randn(n))/sqrt(2);
        [Q2,R2] = qr(Y);
        R2 = diag(diag(R2)./abs(diag(R2)));
        V = Q2*R2;
        U_temp(:,:,i) = U;
        V_temp(:,:,i)= V;
end

save('U', 'U_temp');
save('V', 'V_temp');

%% propogation matrix
load('U.mat');
load('V.mat');
load('tvalues.mat');

for m = 1:NFFT
    M_omega_temp(:,:,m) = eye(n);
    
    for k = 1:n
        C(k) = exp(-1j*omega(m)*t_values(k)/Nseg);
    end
    D = diag(C);   
    
    for i = 1:Nseg
        temp = U_temp(:,:,i)*D*V_temp(:,:,i);
        M_omega_temp(:,:,m) = M_omega_temp(:,:,m)*temp;
    end
    
    M_omega(:,:,m) = M_omega_temp(:,:,m);
end
  
save('Momega', 'M_omega');