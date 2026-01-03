function E=energy(x,ts)
% ENERGY(X,TS) X-Input signal and TS-sampling time
% This function calculates the energy of the discrete time signal
%
% Numerical Integration using Composite Trapezoidal rule
% Computation of energy of a discrete siganl x with sampling time in ts.

y=x.^2;
M=length(x);
val=sum(y(2:M-1));

E=ts*(y(1)+y(M)+2*val)./2;    % Using Trapezoidal Formulae