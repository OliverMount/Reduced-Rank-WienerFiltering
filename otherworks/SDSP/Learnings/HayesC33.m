% program to compare the two auto correlation estimates

clear all;clc;close all;
len=1000;k=100;l=len/k;
noisevar=1;
w=sqrt(noisevar)*randn(1,len); % Part a

% Theoretical Values
r_th=[noisevar zeros(1,len-1)]; % Part b

r=TimeAC(w);
figure,stem(r(1:100),'r*');hold;
stem(r_th(1:100));
title('AutoCorrelation Comparison ');

for m=1:l,
    rx(m,:)=TimeAC(w(k*(m-1)+1:m*k)); % Averaging Estimates
end

rxx=mean(rx);
stem(rxx,'kp');

% Conclusion
% Averaging Autocorrelations have lesser variance than the other estimator