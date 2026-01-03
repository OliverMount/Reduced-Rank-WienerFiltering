function x=AutoR_Kalman(a,noisevar,len)

% This program generates AR(P) Process for a given filter coefficients 
% a(denominator coeffs) and b ( numerator coeffs) and input Gaussian Noise variance. 
% specify the length of the ARMA sequence in len

p=length(a)-1;  % order of the AR --> AR(p) p-Poles
v=zeros(1,len);
% feedv=[1 1 2 -1 2 -2 1];
Xmtx(1,:)=[1 1];%randn(1,p); % Initial State N(0,1);
a1=-a(2:end);

for i=1:len,
    x(i)=Xmtx(i,:)*(a1')+v(i); % convolution
    Xmtx(i+1,:)=[x(i) Xmtx(i,1:p-1)]; % Update Feedback 
end