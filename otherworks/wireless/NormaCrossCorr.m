% Normalized correlation function in rayleigh fading
clear all;clc;close all;
a=0:.1:5;
c=[0 1 1.6];

for i=1:length(c),
co(i,:)=(besselj(0,c(i)))^2./(1+(a.^2));
end
figure,plot(a,co');

b=0:.1:20;
co_a=(besselj(0,b).^2);
figure,plot(b,co_a);