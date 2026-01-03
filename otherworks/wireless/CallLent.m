% Simulation of Modified Jakes (Lent model) Flat fading  for k waveforms

clear all;clc;close all;
l=4; % specification of number of waveforms to be generated
M=2^l; N=4*M; 
W=WalGen(M);
P=100; % no of symbols we want

fm=100; % Maximum Doppler frequency
fs=10*fm;%sampling frequency
ideal=besselj(0,2*pi*fm*(0:P-1)/fs);
n=1:M;
alph=2*pi*(n-.5)/N;beta=(pi*n/M); 
f_n=fm*cos(alph)/fs;

% thet=rand(P,M); % Uniform random variables between (-pi,pi)
for te=1:l,
thet(:,te)=beta'+2*pi*(te-1)/(M+1);
end
b=sqrt(2/M)*cos(beta); c=sqrt(2/M)*sin(beta);
% A=cos(2*pi*(0:P-1)'*f_n );

for q=1:l, % l waveforms
r_I(:,q)=cos(2*pi*(0:P-1)'*f_n + ones(P,1)*thet(:,q)')*((b.*W(l,:))');
r_Q(:,q)=cos(2*pi*(0:P-1)'*f_n + ones(P,1)*thet(:,q)')*((c.*W(l,:))');
end
r=(r_I + j*r_Q); % Low-pass equivalent time correlated multipath samples
z=abs(r);
figure,plot((0:P-1)./fs,10*log10(z(:,3))); % Fading envelope
% axis([0 P/fs -5 15]);

% PP=var(z)+(mean(z)^2);
temp1=TimeAC(r_I(:,3)');
temp2=TimeAC(r_Q(:,3)');
% temp1_n=temp1/(PP/2);
% temp2_n=temp2/(PP/2);
temp1_n=temp1/(mean(r_I(:,3).^2));
temp2_n=temp2/(mean(r_Q(:,3).^2));

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