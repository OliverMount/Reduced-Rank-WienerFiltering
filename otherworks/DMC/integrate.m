function I=integrate(x,ts)
% INTEGRATE(X,TS) X-Input signal and TS-sampling time
% This function calculates the Intergral of a Continuous time signal
%
% Numerical Integration using Composite Trapezoidal rule

y=x.^2;
M=length(x);
val=sum(y(2:M-1));

I=ts*(y(1)+y(M)+2*val)./2;    % Using Trapezoidal Formulae