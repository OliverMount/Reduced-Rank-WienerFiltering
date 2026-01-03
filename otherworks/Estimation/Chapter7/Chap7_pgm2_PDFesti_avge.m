% Program to plot the PDF of the sample mean estimator for the estimation
% of DC (A = 1) in WGN

clear all;clc;close all;
M=5000; % No of averagin iterations
N=20; % data record length
A=1; 
for loop=1:2,
for i=1:M,
    x=A+randn(1,N);    
    Aesti(i)=sum(x)/N; % M realization of a Sample Mean estimator(Random variable)
end

a=max(Aesti);
b=min(Aesti);
range=a-b;
Ncells=100;
Hepsi=(range/(2*Ncells)); % Half of Epsilon
values=(b+Hepsi):(2*Hepsi):a;

for j=1:Ncells,
    temp=Aesti-values(j);
    counter =abs(temp) <= Hepsi;
    histo(j)=sum(counter)/M;
    p(loop,j)=histo(j)/(2*Hepsi);
end
end

figure,plot(values,mean(p));
    
