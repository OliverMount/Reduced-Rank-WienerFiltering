% Kalman filter to track AR(1) process

clear all;clc;close all;
len=100; % Length of the input sequence
a=.99;
var_u=.1; % Process Noise variance
var_w=.9; % Measurement noise variance
s=AutoR_Kalman([1 -a],var_u,len); % Process Model (AR(1))
for i=1:len,
    x(i)=s(i)+sqrt((var_w)^i+1)*randn;
end
% 
figure,plot(1:len,x,1:len,s);
legend('Observed Process','Original Process');

% Initialization

s_cap(1)=0 ; % Mean of the Initial state (Mean of the Parameter)
M(1)=1;

% Kalman filter steps
for j=2:len,

s_cap_prev(j)=(a)*s_cap(j-1);    % 1.Prediction
M_prev(j) =(a^2) *M(j-1) + var_u ; % 2. Prediction MSE
K(j)= M_prev(j)/ (M_prev(j) + (.99)^(j-2)); % 3.Kalman Gain
s_cap(j)=s_cap_prev(j)+ K(j)*(x(j-1)-s_cap_prev(j)); % 4. Correction
M(j)=(1-K(j))*M_prev(j); %5. corrected MSE

end

% Plotting of the MSE

figure,plot(1:len,M,'k*',1:len,M_prev,'r*');
legend('Corrected MSE','Prediction MSE');
title('Correction and Prediction MSE');

% Plotting of Original and Predicted process
figure,plot(1:len,s,'k-',1:len,x,'r-',1:len,s_cap,'b-');
legend('Original Process','Observed process','Estimated process');

figure,plot(1:len,s,'k-',1:len,s_cap,'b-');
legend('Original Process','Estimated process');