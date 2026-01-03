function r=TimeAC(x)

% Asymtotically unbiased estimator
% Divide by M-i for unbiased estimator but high variance

M=length(x);

for i=1:M,
     r(i)=(1/M)*(x(1:M-i+1)*x(i:M)'); 
end