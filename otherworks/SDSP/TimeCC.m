function r=TimeCC(x,y)

% Asymtotically unbiased estimator
% Divide by M-i for unbiased estimator 

M=length(x);

for i=1:M,
     r(i)=(1/M)*(y(1:M-i+1)*x(i:M)'); 
end
