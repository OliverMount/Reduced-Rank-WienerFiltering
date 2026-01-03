% Example 3.5  CRLB for sinusoidal frequency estimation
% Refer page 36 of Steven Kay.

clear all;clc;close all;

n=0:9; % Data record Length
phi=0;
SNR=1:1000;
temp=2*pi*n;
f=0.25; % some preferred frequencies

for i=1:length(SNR),
   CRLB(i)=1/((SNR(i))*sum((sin(2*pi*n*f+phi).*temp).^2));
end

figure,plot(10*log10(SNR),CRLB);
title('CRLB for Sinusoidal Frequency Estimation');
xlabel('SNR(in dB)');
clear SNR;