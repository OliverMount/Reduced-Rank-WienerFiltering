function [w,e]=assign2LMSNew(u,d,delta)
%Adaptive filter part
% Initialization
N=2; % length of the adaptive filter
x=convolmtx(u,N);
w(1,:)=zeros(1,N); % zero Wight initialization
x=[zeros(1,N);x]; % delaying by one sample for the predictor
y(1)=x(1,:)*w(1,:)';
% o/p of adaptive filter at first time instant
e(1)=d(1)-y(1); % error signal
w(2,:)=w(1,:)+(delta*(x(1,:)*e(1))); % weight adapt at first time instant
z=1500; % no of iterations
% iterative algorithm
for n=2:z,
y(n)=x(n,:)*w(n,:)';
e(n)=d(n)-y(n);
w(n+1,:)=w(n,:)+(delta*(x(n,:)*e(n)));
end
clear n;