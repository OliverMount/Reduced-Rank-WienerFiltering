% Example 3.5  CRLB for sinusoidal frequency estimation
% Refer page 36 of Steven Kay.

clear all;clc;close all;

n=0:9; % Data record Length
phi=0;
SNR=1;
temp=2*pi*n;
f=0.01:.001:0.49;

for i=1:length(f),
   CRLB(i)=1/sum((sin(2*pi*n*f(i)+phi).*temp).^2);
%    S(i,:)=sin(2*pi*f(i)*n + phi);
end

figure,plot(f,CRLB);
title('CRLB for Sinusoidal Frequency Estimation');