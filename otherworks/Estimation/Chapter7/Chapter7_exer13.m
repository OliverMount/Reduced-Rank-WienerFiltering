% Monte-carlo simulations
% Program to plot the PDF of the sample mean estimator for the estimation
% of DC(A) in WGN

clear all;clc;close all;
M=5000; % No of averaging iterations
N=50; % data record length
A=1; 
vari=.1;

for i=1:M,
    x=A+sqrt(vari)*randn(1,N);    
    Aesti(i)=sum(x)/N; % M realization of a Sample Mean estimator(Random variable)
end

% Finding the PDF
a=max(Aesti);
b=min(Aesti);
range=a-b;
Ncells=100;
Hepsi=(range/(2*Ncells));
values=(b+Hepsi):(2*Hepsi):a; % middle values of observations around epsilon

for j=1:Ncells,
    temp=Aesti-values(j);
    counter =abs(temp) <= Hepsi;
    histo(j)=sum(counter)/M; % Histogram Estimator
    p(j)=histo(j)/(2*Hepsi); % PDF estimator
end

% Theoretical  Asymptotic PDF of the estimator for A=1

variance=vari/N;
ptheory=inv(sqrt(2*pi*variance))*exp(-(values-A).^2/(2*variance)); 

figure,hnd=plot(values,p,values,ptheory);
set(hnd,'linewidth',2);

hnd=title('PDF of Sample Mean Estimator');
set(hnd,'fontsize',15);
hnd=xlabel(' Random variable range ');
set(hnd,'fontsize',15);
hnd=ylabel('  Amplitude  ');
set(hnd,'fontsize',15);  