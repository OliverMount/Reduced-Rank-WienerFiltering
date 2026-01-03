% This program for Sequential Weighted Least squares (SWLS)

% source Generation
clear all;clc;close all;
n=0:299;
N=length(n);
theta_org=[1 .03]'; % Original parameters
W=eye(N); % Weighting Matrix


H=[ones(N,1) n'];
% H=[ones(N,1) n' n'.^2 n'.^3];

% Model for Source
x=H*theta_org + .25*randn(1,N)';  % variance of the Noise .0625

p=2; % No of parameters to estimate

% Initialization
sigma=inv(H(1,:)*W(1,1)*H(1,:)');    % Covariance Mtx in the case of BLUE
theta=(sigma*H(1,:)*W(1,1)*x(1))';   % Initial estimator based on one sample
% Jmin(1)=x(1)*(x(1)-H(1,:)*theta);    % First time error

for i=2:N,
    
    a=(1/W(i,i))+H(i,:)*sigma*H(i,:)'; % changed H'(i,:) to H(i,:) and Vice versa in the original algo
    K(:,i)=sigma*H(i,:)'/a;      % finding the Gain
    e(i)=x(i)-H(i,:)*theta;    % Prediction error
    theta=theta+K(:,i)*e(i);    % The estimator update
%     Jmin(i)=Jmin(i-1)+ e(i)/a;   % Minimum Error
    sigma=[eye(p)-K(:,i)*H(i,:)]*sigma;    % covariance update
    
end

figure,hnd=plot((e.^2));
hnd=title('Error of Recursive Least Squares');
set(hnd,'fontsize',15); 
hnd=xlabel('Time (t)');
set(hnd,'fontsize',15); 
