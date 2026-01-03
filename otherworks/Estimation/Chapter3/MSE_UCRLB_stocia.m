% Stocia paper example for biased estimators Vs UCRLB

clear all; clc ;close all
sigma=1;
N=1:300;
a=2*(sigma^2)./N; %UCRLB
b=2*(sigma^2)./(N+2);  % MSE of biased estimator
plot(N,a,N,b);
% axis([0 300 -.0001 .04]);
legend('Variance 1','Variance 2');