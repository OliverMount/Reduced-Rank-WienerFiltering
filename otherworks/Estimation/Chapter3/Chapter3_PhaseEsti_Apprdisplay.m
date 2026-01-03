% Approximation display in single Phase estimation 
% summation is zero except at 0 and .5 normalized frequency 

clear all;clc;close all;
fo=0:.01:.5;
N=50; % No of data samples
n=0:N-1; phi=0.02*pi; % Phi in radiance
noisevar=0.1;

 for l=1:length(fo),
     x(l,:)=cos(4*pi*fo(l)*n+(2*phi));
     y(l)= (1/N)*sum(x(l,:));
 end

figure,plot(fo,y);