% Assignment 3 - For Task 2: Effect of StepSize:

clear all;
clc;
close all;

for n=1:3,
h(n)=.5*(1+cos(2*pi*(n-2)/3.1));
end
mu=[.0025 .0075 .075];
for varyW=1:length(mu),
for ensemble=1:200,
L=1000; % Length of the desired signal
%Desired signal genration
d=randint(1,L);
d=2*d-1; % Mean zero Variance 1.
d_delayed=[zeros(1,7) d(1:L-7)];
%Channel Output.
mtx=convolmtx(d,length(h));
z=mtx*h';
v=sqrt(.001)*randint(1,length(z));
u=z+v';
u=u';
[w,e]=assign3LMS(u,d_delayed,mu(varyW));
err(ensemble,:)=e;
store_w(ensemble,:)=w(end,:);
end
Avg_weight(varyW,:)=mean(store_w);
final(varyW,:)=(mean(abs(err))).^2;
end
figure,semilogy(final');
legend('mu=.0025','mu=0.0075','mu=0.075');
hnd=title('Error convergence plot');set(hnd,'Fontsize',15);
hnd=xlabel('Number of iterations --> ');set(hnd,'Fontsize',15);
hnd=ylabel('Mean square error --> ');set(hnd,'Fontsize',15);

