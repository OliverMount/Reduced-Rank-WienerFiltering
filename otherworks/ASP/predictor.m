%Specific program for second order predictor
clear all;
clc;
close all;
a1=0.1;a2=-0.8;
w1=-a1;w2=-a2;
for ensemble=1:200,
v=sqrt(0.27)*randn(1,1500);
% Generation of the AR2 process
u(1)=v(1);
u(2)=w1*u(1)+v(2);
for n=3:length(v),
u(n)=w1*u(n-1)+w2*u(n-2)+v(n);
end
% only for AR2 process
lamda1=(1-(a1/(1+a2)));
lamda2=(1+(a1/(1+a2)));
mu=.0099; % Selection of step size
d=u;
[w,e]=assign2LMSNew(u,d,mu);
err(ensemble,:)=e;
wei(ensemble,:)=w(end,:);
end
final=mean(wei);
str=sprintf('The final filter weights are w1= %f w2= %f ', final(1,1) , final(1,2));
disp(str);
pred_error=(sum(abs(err))/((ensemble))).^2;
figure,plot(10*log10(pred_error));
title('Prediction Error power plot');
xlabel('Number of iterations --> ');
ylabel('Mean square error in dB --> ');
final_weights=w;
figure,plot(final_weights);
legend('First tap','Second tap');
title('Wight vector convergence plot');
xlabel('Number of iterations --> ');
ylabel('Weight values --> ');