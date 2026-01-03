% MATLAB Program:
% Main Function :
clear all;clc;close all;
mu=[.025 .05 0.075 ]; % step sizes
for o=1:length(mu),
for ensemble=1:100,
u=randint(1,1000);
u=2*u-1;
% figure,stem(u);
% axis([1 20 -1.5 1.5]); % showing only the first 20 bits.
h=[.4 1 -.7 .2]; % Filter coefficients
mtx=convolmtx(u,length(h)); % input signal
d=mtx*h'; % Desired signal
% d=u;
[M N]=size(mtx);
% figure,plot(d);
%Adaptive filter part
% Initialization
N=4;
x=convolmtx(u,N);
w=zeros(1,N); % Adaptive filter coeffient
y(1)=x(1,:)*w'; % o/p of adaptive filter at first time instant
e(1)=d(1)-y(1); % error signal
w(1,:)=w+(mu(o)*(x(1,:)*e(1)));

% iterative algorithm
for n=2:1000,
y(n)=x(n,:)*w(n-1,:)';
e(ensemble,n)=d(n)-y(n);
w(n,:)=w(n-1,:)+(mu(o)*(x(n,:)*e(ensemble,n)));
end
clear n;
end

final(o,:)=(sum(e)/100).^2;
% figure,plot(w);
% legend('First tap','Second tap','Third tap','Fourth tap');
end
figure,semilogy(final');
legend('mu = .025','mu = .05','mu = .075');

