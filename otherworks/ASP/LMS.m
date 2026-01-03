function [w,e]=LMS(u,d,mu,N,z)
x=convolmtx(u,N);
w(1,:)=zeros(1,N); % zero Wight initialization
y(1)=x(1,:)*w(1,:)'; % o/p of adaptive filter at first time instant
e(1)=d(1)-y(1); % error signal
w(2,:)=w(1,:)+(mu*(x(1,:)*e(1))); % weight adapt at first time instant
% iterative algorithm
for n=2:z,
y(n)=x(n,:)*w(n,:)';
e(n)=d(n)-y(n);
w(n+1,:)=w(n,:)+(mu*(x(n,:)*e(n)));
end