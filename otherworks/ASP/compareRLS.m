function [w1 ,e1]=compareRLS(H,N,x)

W=eye(N);
% Initialization
sigma=inv(H(1,:)*W(1,1)*H(1,:)');    % Covariance Mtx in the case of BLUE
theta=(sigma*H(1,:)*W(1,1)*x(1))';   % Initial estimator based on one sample
Jmin(1)=x(1)*(x(1)-H(1,:)*theta);    % First time error
w1(1,:)=theta';
p=length(theta);
for i=2:N,
    
    a=(1/W(i,i))+H(i,:)*sigma*H(i,:)'; % changed H'(i,:) to H(i,:) and Vice versa in the original algo
    K(:,i)=sigma*H(i,:)'/a;      % finding the Gain
    e(i)=x(i)-H(i,:)*theta;    % Prediction error
    theta=theta+K(:,i)*e(i);    % The estimator update
%     disp(theta);
    Jmin(i)=Jmin(i-1)+ e(i)/a;   % Minimum Error
    sigma=[eye(p)-K(:,i)*H(i,:)]*sigma;    % covariance update
    w1(i,:)=theta';
end
e1=e;