clear all;clc;close all;
M=5000; % No of averagin iterations
N=20; % data record length
A=1;

for i=1:M,
    x=A+randn(1,N);    
    Aesti(i)=-.5+sqrt((sum(x.^2)/N)+.25); % M realization of a ML estimator
end
figure,plot(Aesti); 

Aesti_mean=sum(Aesti)/M;
Aesti_var=sum((Aesti-Aesti_mean).^2)/M;

disp(Aesti_mean);
disp(N*Aesti_var);