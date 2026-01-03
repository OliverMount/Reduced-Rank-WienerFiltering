% For Newton and Newton-LMS:
% Main program:
% Adaptive Noise Cancellation
clear all;clc;close all;
N=16;
len=63*N; % Length of the input seqence
adapt_ord=2; % Adaptive Filter Order
nvar=0.01;
R=[0.5+nvar 0.5*cos(2*pi/N);0.5*cos(2*pi/N) .5+nvar];
P=[0 -sin(2*pi/N)]';dvar=2;
wei_aver=zeros(len+1,adapt_ord)';
len_ensemble=200;
for ensemble=1:len_ensemble,
for n=0:len,% 63 Cycles of Sinusoid with N=16 period
i=n+1;
x(i)=sin(2*pi*n/N);
d(i)=2*cos(2*pi*n/N);
end
v=sqrt(.01)*randn(1,len+1);
u=x+v; % Input to the adaptive filter
% Adaptive Part
adapt_ord=2;no_of_iter=len;mu=.05;
[w,e]=NEWTON_LMS(u,d,mu,adapt_ord,len,R);
% [w,e]=NEWTON(u,d,mu,adapt_ord,len,R,P);
err(:,ensemble)=e';
wei(:,ensemble)=w(:,end);
wei_aver=wei_aver+w;
end
final=mean(wei');
str=sprintf('The final filter weights are w1= %f w2= %f ', final(1,1) , final(1,2));
disp(str);
w_av=wei_aver'/len_ensemble;
% Error Surface plot
v=[0.49 2 6.3 15 25] ; % Constnt values of Surface
w0=-2:.1:8;
w1=-10:0.1:0;
for i=1:length(w0),
for j=1:length(w1),
J(i,j)=((.5+0.01)*(w0(i)^2+w1(j)^2))+(cos(2*pi/N)*w0(i)*w1(j))+(2*w1(j)*sin(2*pi/N));
+2;
end
end
clear i j;
for j=1:length(w_av(:,2)),
Jcap(j)=dvar-(2*P'*[w_av(j,1);w_av(j,2)])+([w_av(j,1) w_av(j,2)]*R*[w_av(j,1);
w_av(j,2)]);
end
figure,[c,h]=contour(w0,w1,J,v);
title('Surface Contours and Weight-value tracks');
clabel(c,h,v,'labelspacing',600);
xlabel('W0'); ylabel('W1');
hnd=plot3(w_av(:,1),w_av(:,2),Jcap,'k-');