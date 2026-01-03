% Numerical Determination of MLE Newton-Raphson Approach
% Example problem in book (7.11) (PAge 180-181)
clear all;clc;close all;

N=50; % Data record length
n=0:N-1;
r=0.5; % Parameter to be estimated
noivar=0.01; % Noise Variance
x= (r).^n + sqrt(noivar)*randn(1,N);

% Plot of the cost function
rp=0:.05:1;
for i=1:length(rp),
J(i)=-sum((x-(rp(i)).^n).^2);
end
% figure,plot(rp,J);axis([0 1 -3 0]); % Function to be maximised 

% Newton Raphson method
noiter=30;
r_e(1)= 0.2 ;% 0.2; 1.2; % Initial estimate

for k=2:noiter,
    numer=sum((x-(r_e(k-1)).^n).*n.*(r_e(k-1).^(n-1)));
    denomi= sum(n.*(r_e(k-1).^(n-2)).*(((n-1).*x)-((2*n-1).*(r_e(k-1).^n))));
    r_e(k)=r_e(k-1) - (numer/denomi);
end
% figure,plot(1:noiter,r_e,1:noiter,ones(1,noiter)*r);

% Method oF Scoring (replacing the second derivative with its expected
% value will increase the stability of the iteration

% noiter=30;
r_es(1)= 1.2 ;% 0.2; 1.2;

for k=2:noiter,
    numer=sum((x-(r_es(k-1)).^n).*n.*(r_es(k-1).^(n-1)));
    denomi= sum(n.^(2).*(r_es(k-1).^(2*n-2)));
    r_es(k)=r_es(k-1) + (numer/denomi);
end
figure,plot(1:noiter,r_es,1:noiter,ones(1,noiter)*r,1:noiter,r_e);
legend('By Scoring Approach','Original Parameter','Newton-Raphson method');