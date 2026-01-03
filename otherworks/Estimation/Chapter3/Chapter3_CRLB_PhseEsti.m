% CRLB for PHASE estimation in AWGN

clear all;clc;close all;
fo=0:.01:.5;
N=50; % No of data samples
n=0:N-1; phi=0.2*pi; % Phi in radiance
noisevar=1;
A=1;

 for l=1:length(fo),
     x(l,:)=cos(4*pi*fo(l)*n+(2*phi));
     y(l)= (noisevar)/((A^2*(N/2)))*(1-(1/N)*sum(x(l,:))); % Original CRLB
     z(l)=(2*noisevar)/(N*(A^2)); % Approximate CRLB
 end

figure,plot(fo,y,fo,z);legend('Original','Approximation');