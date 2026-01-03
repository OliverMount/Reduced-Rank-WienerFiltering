% Estimation of Known sinusodal amplitudes (Fporier series example in
% chapeter 4. We can find complete statistics of a and s 
clear all;clc;close all;
para=[1 .5 1 .5]';
M=length(para); %No of parameters
N=16;
n=0:5*N-1;
H=[cos(2*pi*n/N) ; cos(2*pi*n*3/N) ; sin(2*pi*n/N) ; sin(2*pi*n*3/N)]';
% s=cos(2*pi*n/8)+ .5*cos(2*pi*n*3/8) + sin(2*pi*n/8) + .5*sin(2*pi*n*3/8); %Here only four parameters
s=H*para;
noi=randn(1,length(s));
x=s+noi.';

a=inv(H'*H)*H'*x;

estimated_signal=H*a;
figure,plot(s);hold;
plot(x,'k');
plot(estimated_signal,'r')
legend('Original','Observed','Estimated');
hnd=title('Known sinusoids in Noise');
set(hnd,'fontsize',15); 

SNR=10*log10(var(s)/var(noi))