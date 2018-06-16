%two mode fibers
%six mode fiber
clear all
lambda = -6:12/1000:6;
len = size(lambda);
pdf_six_mode = zeros(1,len(2));

%f = (5*lambda.^2);
for i = 1:1:len(2)
    a = sqrt(30)*exp(-6/5*((lambda(i)^2)))/sqrt(pi);
    b = (13436928*(lambda(i))^10/244140625 - 4478976*(lambda(i))^8/9765625);
    c = (2581632*(lambda(i))^6/1953125 - 102816*(lambda(i))^4/78125);
    d = (7812*(lambda(i))^2/15625 + 644/15625);
    f = a*(b+c+d);
    pdf_six_mode(i) = f;
    %i = i+1;
    %append(pdf_six_mode, f);
end
plot(lambda, pdf_six_mode/norm(pdf_six_mode));
xlabel('Normalised group delay, \tau');
ylabel('p.d.f');
[pks, locs] = findpeaks(pdf_six_mode);
t_values = lambda(locs)
