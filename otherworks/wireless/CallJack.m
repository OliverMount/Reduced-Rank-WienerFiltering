% simulation of Jakes Flat fading 
clear all;clc;close all;
N=34;% 34( 8 Oscillators) 66(16 Oscillators)
P=1000; % no of symbols we want
M=.5*((N/2)-1); % one less than no of No of oscillators
fm=100; % Maximum Doppler shift (Page 42 of stuber);
fs=10000;%sampling frequency

n=1:M;
f_n=fm*cos(2*pi*n/N)/fs;
alpha=0;beta=(pi*n/M);
ideal=besselj(0,2*pi*fm*(0:P-1)/fs);
for l=1:M,
    A(:,l)=cos(2*pi*f_n(l)*(0:P-1))';
    b(l)=2*cos(pi*l/M); c(l)=2*sin(pi*l/M);
end

r_I=A*b'+(sqrt(2)*cos(2*pi*(fm/fs)*(0:P-1)))';
r_Q=A*c' ;
r=(1/sqrt((2*M)+1))*(r_I+j*r_Q);

new_AC=TimeAC(r.');

kkk=fm*(0:P-1)/fs;
figure,plot(kkk,real(new_AC),kkk,ideal);
% figure,plot(1:length(r),real(new_AC),1:length(r),ideal);
legend('Simulated ','Ideal');
axis([0 14 -1 1.1]);

% fading envelope
% figure,plot(10*log10(abs(r)));
figure,semilogy((0:P-1),(abs(r)));
% figure,plot((0:P-1),10*log10(abs(r)));
hnd=title(' Fading Envelope ');
set(hnd,'fontsize',15);
axis([0 P-1 1e-2 5]);
