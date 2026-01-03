clear all;clc;close all;
y=0:.2:10;
avg=.2; % Mean of the original pdf 
a=abs(y-avg);
b=abs(-y-avg);

f=.5*(exp(-a) + exp(-b)); % PDF after tranformation

figure,plot(y,f);
figure,plot(y,.5*(exp(-a)));
figure,plot(y,.5*exp(-b));

% Plot of the Original Laplacian PDF

% x=-8:.2:8;
% a=-5;
% b=1; % variance is 2/(b^2) for Laplacian PDF
% 
% p=(b/2)*exp(-b*abs(x-a));
% 
% figure,plot(x,p);




