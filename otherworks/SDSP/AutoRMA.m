function x=AutoRMA(a,b,noisevar,len)

% This program generates ARMA(p,q) Process for a given filter coefficients 
% a(denominator coeffs) and b ( numerator coeffs) and input Gaussian Noise variance. 
% specify the length of the ARMA sequence in len
% x is ARMA output and v is the Noise input

p=length(a)-1;q=length(b)-1;  % order of the ARMA --> ARMA(p,q) p-Poles q-Zeros
v=sqrt(noisevar)*randn(1,len);
% v=[1 zeros(1,len-1)]; % To find impulse response
%v=[1 1 2 -1 2 -2 1]; x=[1 1.2 2.5 -.11 2.745 -1.2495 1.52425]

Vmtx=convolmtx(v,length(b));
Xmtx(1,:)=zeros(1,p);
a1=-a(2:end);

for i=1:len,
    x(i)=Xmtx(i,:)*(a1)+Vmtx(i,:)*b; % convolution
    Xmtx(i+1,:)=[x(i) Xmtx(i,1:p-1)]; % Update Feedback 
end