% Monte-carlo simulations
% Program to plot the PDF of the Slutsky's Theorem for the estimation
% of DC(A) in WGN

clear all;clc;close all;
M=5000; % No of averaging iterations
N=100; % data record length
A=0; 
vari=25;

for i=1:M,
    x=A+sqrt(vari)*randn(1,N);    % IID Zero mean rv's
    Aesti(i)=sum(x)/N; % M realization of a Sample Mean estimator(Random variable)
    sig2esti(i)=sum((x-Aesti(i)).^2)/N; % M realization of a Sample Variance estimator(Random variable)
    newrv(i)=Aesti(i)/(sqrt(sig2esti(i)/N));
end

% Finding the PDF
a=max(sig2esti);
b=min(sig2esti);
range=a-b;
Ncells=100;
Hepsi=(range/(2*Ncells));
values=(b+Hepsi):(2*Hepsi):a;

for j=1:Ncells,
    temp=sig2esti-values(j);
    counter =abs(temp) <= Hepsi;
    histo(j)=sum(counter)/M; % Histogram Estimator
    p(j)=histo(j)/(2*Hepsi); % PDF estimator
end

% Theoretical  Asymptotic PDF of the estimator for A=1

variance=2*((vari)^2)/N;
ptheory=inv(sqrt(2*pi*variance))*exp(-(values-vari).^2/(2*variance)); 

figure,hnd=plot(values,p,values,ptheory);
set(hnd,'linewidth',2);

hnd=title('PDF of Sample Variance Estimator');
set(hnd,'fontsize',15);
hnd=xlabel(' Random variable range ');
set(hnd,'fontsize',15);
hnd=ylabel('  Amplitude  ');
set(hnd,'fontsize',15); 
clear Hepsi values ptheory;

% Finding the PDF of the Slutsky's theorem
a=max(newrv);
b=min(newrv);
range=a-b;
Ncells=100;
Hepsi=(range/(2*Ncells));
values=(b+Hepsi):(2*Hepsi):a;

for j=1:Ncells,
    temp=newrv-values(j);
    counter =abs(temp) <= Hepsi;
    histo(j)=sum(counter)/M; % Histogram Estimator
    p(j)=histo(j)/(2*Hepsi); % PDF estimator
end

% Theoretical  Asymptotic PDF of the estimator for A=1

variance=1;
ptheory=inv(sqrt(2*pi*variance))*exp(-(values).^2/(2*variance)); 

figure,hnd=plot(values,p,values,ptheory);
set(hnd,'linewidth',2);

hnd=title('PDF of RV in Slutskys theorem');
set(hnd,'fontsize',15);
hnd=xlabel(' Random variable range ');
set(hnd,'fontsize',15);
hnd=ylabel('  Amplitude  ');
set(hnd,'fontsize',15); 