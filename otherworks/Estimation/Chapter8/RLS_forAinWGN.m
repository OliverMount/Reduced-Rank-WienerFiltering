% This program for Sequential Weighted Least squares (SWLS) for WGN ( Single
% Parameter Sequntial Least Squares )

% source Generation
clear all;clc;close all;
n=0:299;
N=length(n);
A=10;
theta_org=A; % Original parameters
W=eye(N); % Weighting Matrix

H=[ones(N,1)];

% Model for Source
x=H*theta_org + randn(1,N)';  % unit variance WG
p=1; % No of parameters to estimate

% Initialization
sigma(1)=inv(H(1,:)'*W(1,1)*H(1,:));    % Covariance Mtx in the case of BLUE
theta(1)=sigma(1)*H(1,:)'*W(1,1)*x(1);   % Initial estimator based on one sample
% Jmin(1)=x(1)*(x(1)-H(1,:)*theta);    % First time error

for i=2:N,
    
    a=(1/W(i,i))+H(i,:)*sigma(i-1)*H(i,:)'; % changed H'(i,:) to H(i,:) and Vice versa in the original algo
    K(i)=sigma(i-1)*H(i,:)'/a;      % finding the Gain
    e(i)=x(i)-H(i,:)'*theta(i-1);    % Prediction error
    theta(i)=theta(i-1)+K(i)*e(i);    % The estimator update
%     Jmin(i)=Jmin(i-1)+ e(i)/a;   % Minimum Error
    sigma(i)=[eye(p)-K(i)*H(i,:)']*sigma(i-1);  % covariance update
    
end

figure,plot(n,K);xlabel('Current sample, N');ylabel('Gain');axis([-2 length(n)-200 0 .6]);
figure,plot(n,sigma);xlabel('Current sample, N');ylabel('Variance');axis([-2 length(n)-200 0 1.1]);
figure,plot(n,theta,n,ones(1,length(n))*A);legend('Estimated','Original');xlabel('Current sample, N');ylabel('Estimate');axis([-2 length(n)-100 9.2 10.4]);