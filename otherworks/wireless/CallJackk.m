% simulation of Jakes Flat fading  for M waveforms
% Maximum we can generate M waveforms by this method
clear all;clc;close all;
N=34;% 34( 8 Oscillators) 66(16 Oscillators)
P=100; % no of symbols we want
M=.5*((N/2)-1); % one less than no of No of oscillators
fm=100; % Maximum Doppler frequency
fs=10*fm;%sampling frequency
n=1:M;
f_n=fm*cos(2*pi*n/N)/fs;
alpha=0;beta=(pi*n/(M)); gamm=2*pi*(0:M-1)'*n/M;
gamm_c=cos(gamm); gamm_s=sin(gamm);

ideal=besselj(0,2*pi*fm*(0:P-1)/fs);

A=cos(2*pi*(0:P-1)'*f_n); B=sin(2*pi*(0:P-1)'*f_n);
b=2*cos(pi*(n)/M); c=2*sin(pi*(n)/M);
for q=1:M, % M waveforms
r_I(:,q)=A*(gamm_c(q,:).*b)'+B*(gamm_s(q,:).*(-b))'+(sqrt(2)*cos(2*pi*(fm/fs)*(0:P-1)))';
r_Q(:,q)=A*(gamm_c(q,:).*c)'+B*(gamm_s(q,:).*(-c))';
end
r=(r_I+j*r_Q); % Low-pass equivalent time correlated multipath samples
z=abs(r);
figure,plot((0:P-1)./fs,10*log10(z(:,1))); % Fading envelope
% axis([0 P/fs -5 15]);% % axis([0 P/fs -5 15]);
% 
% PP=var(z)+(mean(z)^2);
temp1=TimeAC(r_I(:,1)');
temp2=TimeAC(r_Q(:,1)');
% temp1_n=temp1/(PP/2);
% temp2_n=temp2/(PP/2);
temp1_n=temp1/(mean(r_I(:,1).^2));
temp2_n=temp2/(mean(r_Q(:,1).^2));

kkk=fm*(0:P-1)/fs;
figure,plot(kkk,temp1_n,kkk,temp2_n,kkk,ideal);
legend('Simulated I-phase','Simulated Q-phase','Ideal');
axis([0 10 -1 1]);

figure,plot(kkk,temp1_n,kkk,ideal);
legend('Simulated I-phase','Ideal');
axis([0 10 -1 1]);

figure,plot(kkk,temp2_n,kkk,ideal);
legend('Simulated Q-phase','Ideal');
axis([0 10 -1 1]);%