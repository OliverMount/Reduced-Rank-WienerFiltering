function x=MovA(b,noisevar,len)

% This program generates MA(q) Process for a given filter coefficients 
%  b ( numerator coeffs) and input Gaussian Noise variance. 
% specify the length of the MA sequence in len
% x is MA output and v is the Noise input

v=sqrt(noisevar)*randn(1,len);
% v1=[1 2 3];v=[v1 zeros(1,len-length(v1)-1)];
% v=[1 zeros(1,len-1)]; % To find impulse response
x=convol(v,b)';