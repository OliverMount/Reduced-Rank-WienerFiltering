% This program for comparison of  RLS with LMS

% source Generation
clear all;clc;close all;
noisevar=.001;
mu=.075;
W=3.3;
% Channel Description
for n=1:3,
h(n)=.5*(1+cos(2*pi*(n-2)/W));
end

for ensemble=1:50,
L=1000; % Length of the desired signal
%Desired signal genration
d=randint(1,L);
d=2*d-1; % Mean zero Variance 1.
d_delayed=[zeros(1,7) d(1:L-7)];
%Channel Output.
mtx=convolmtx(d,length(h));
z=mtx*h';
v=sqrt(noisevar)*randint(1,length(z)); % Noise Variance .001
u=z+v';
u=u';
[w,e]=compareLMs(u,d_delayed,mu,L);
[w1,e1]=compareRLS(mtx,L,u);
err(ensemble,:)=e;
err_RLS(ensemble,:)=e1;
store_w(ensemble,:)=w(end,:);
store_W_RLS(ensemble,:)=w1(end,:);
end

figure,semilogy(mean(abs(err).^2));hold
hnd=title('Error convergence plot');
set(hnd,'fontsize',15);
hnd=xlabel('Number of iterations --> ');
set(hnd,'fontsize',15);
hnd=ylabel('Mean square error --> ');
set(hnd,'fontsize',15);


% semilogy(mean(err_RLS).^2);
semilogy(mean(abs(err_RLS)).^2,'r-');
legend('LMS','RLS');