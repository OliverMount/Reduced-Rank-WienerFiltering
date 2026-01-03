function [w,s]=LMSDDM(w0,yy,mu,N)

% LMS in Decision directed mode
% u --> input to LMS   d --> desired signal  mu --> step size  w0=Intial
% weigh from the training mode
% N --> Adaptive filter order 
% w --> Weight vector  e --> Estimation error  y --> output of adap filter

x=convolmtx(yy,N);
w(1,:)=w0; %  Weight initialization
y(1)=x(1,:)*w(1,:)'; % o/p of adaptive filter at first time instant
% s(1)=y(1) > 0; % Decision box
% s(1)=2*s(1)-1;

if y(1)>0,
    s(1)=1;
else s(1)=-1;
end
e(1)=s(1)-y(1); % error signal
w(2,:)=w(1,:)+(mu*(x(1,:)*e(1))); % weight adapt at first time instant
z=length(yy); % no of iterations
% iterative algorithm
for n=2:z,
y(n)=x(n,:)*w(n,:)';
% s(n)= y(n)>0;
% s(n)=2*s(n)-1;

if y(n)>0,
    s(n)=1;
else s(n)=-1;
end
    
e(n)=s(n)-y(n);
w(n+1,:)=w(n,:)+(mu*(x(n,:)*e(n)));
end