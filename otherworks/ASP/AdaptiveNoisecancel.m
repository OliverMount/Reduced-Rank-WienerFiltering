% For LMS and NLMS:
% Main program:
% Adaptive Noise Cancellation
len_ensemble=100;
N=16;
len=63*N; % Length of the input seqence
adapt_ord=2; % Adaptive Filter Order
wei_aver=zeros(len+1,adapt_ord);
for ensemble=1:len_ensemble,
for n=0:len,% len Cycles of Sinusoid with N=16 period
i=n+1;
x(i)=sin(2*pi*n/N);
d(i)=2*cos(2*pi*n/N);
end
v=sqrt(.01)*randn(1,len+1);
u=x+v; % Input to the adaptive filter
% Adaptive Part
no_of_iter=len;mu=0.01;
[w,e]=NLMS(u,d,mu,adapt_ord,len);
[w,e]=LMS(u,d,mu,adapt_ord,len); % mu=.1
err(ensemble,:)=e;
wei(ensemble,:)=w(end,:);
wei_aver=wei_aver+w;
end
final=mean(wei);
str=sprintf('The final filter weights are w1= %f w2= %f ', final(1,1) , final(1,2));
disp(str);
pred_error=abs(sum(err)/ensemble).^2;
w_av=wei_aver/len_ensemble;
nvar=.01;
R=[0.5+nvar 0.5*cos(2*pi/N);0.5*cos(2*pi/N) .5+nvar];
P=[0 -sin(2*pi/N)]';dvar=2;
% Wopt And Jmin calculation
wopt=inv(R)*P;
Jmin=dvar-wopt'*R*wopt;
% Error Surface plot
v=[0.49 2 6.3 15 25] ; % Constnt values of Surface
w0=-2:.1:8;
w1=-10:0.1:0;
for i=1:length(w(:,1)),
for j=1:length(w(:,2)),
J(i,j)=Jmin+([wopt(1)-w_av(j,1) wopt(2)-w_av(j,2)]*R*[wopt(1)-w_av(j,1) ;
wopt(2)-w_av(j,2)]);
end
end
clear i j;
for j=1:length(w(:,2)),
Jcap(j)=dvar-(2*P'*[w_av(j,1);w_av(j,2)])+([w_av(j,1) w_av(j,2)]*R*[w_av(j,1);
w_av(j,2)]);
end
figure,[c,h]=contour(w0,w1,J,v);
title('Surface Contours and Weight-value tracks');hold on;
clabel(c,h,v,'labelspacing',600);
xlabel('W0'); ylabel('W1');
plot3(w_av(:,1),w_av(:,2),Jcap,'k-');