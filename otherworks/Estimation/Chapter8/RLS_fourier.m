% This program for Sequential Weighted Least squares (SWLS)

% source Generation
clear all;clc;close all;
n=0:299;
N=length(n);
W=eye(N); % Weighting Matrix
M=1;
fo=.1;

Theta_org=[1 2]'; % Original parameters
p=length(Theta_org); % No of parameters to estimate
H=[cos(2*pi*fo*n)' sin(2*pi*fo*n)']; % The system matrix
S=H*Theta_org;
% Model for Source (data observed)
x= S + randn(N,1);  % variance of the Noise is 1

% Initialization
sigma=inv(H(1:p,:)'*inv(W(1:p,1:p))*H(1:p ,:));    % Covariance Mtx in the case of BLUE (p by p matrix)
theta(:,1)=(sigma*H(1:p,:)*W(1,1)*x(1:p));   % Initial estimator based on one sample (p by 1 vector)
% Jmin(1)=x(1)*(x(1)-H(1,:)*theta);    % First time error

for i=1:N-p,
    
    a=(1/W(i,i))+H(p+i,:)*sigma*(H(p+i,:).'); % changed H'(i,:) to H(i,:) and Vice versa in the original algo
    K(:,p+i)=sigma*H(p+i,:).'/a;      % finding the Gain (p by 1 vector)
    e(p+i)=x(p+i)-H(p+i,:)*theta(:,i);    % Prediction error
    theta(:,i+1)=theta(:,i)+K(:,i)*e(p+i);    % The estimator update
%     Jmin(i)=Jmin(i-1)+ e(i)/a;   % Minimum Error
    sigma=[eye(p)-K(:,i)*H(p+i,:)]*sigma;    % covariance update 
    
end

% figure,plot(K');xlabel('Current sample, N');ylabel('Gain');legend('a','b');%axis([-2 length(n)-200 0 .6]);
% % figure,plot(sigma);xlabel('Current sample, N');ylabel('Variance');%axis([-2 length(n)-200 0 1.1]);
% figure,plot([[theta theta(:,end)]'  ones(length(n),p)*diag(Theta_org)]);xlabel('Current sample, N');ylabel('Estimate');%axis([-2 length(n)-100 9.2 10.4]);


theta(:,end)
Theta_org