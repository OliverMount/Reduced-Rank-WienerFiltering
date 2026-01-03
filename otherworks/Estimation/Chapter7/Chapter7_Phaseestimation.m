% Simusoidal phase estimator (Single parameter estimation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This portion for calulalting no of samples (N) we want
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;clc;close all;
N=20:5:100;
A=1; % knowm Amplitude
fo=.08; % Normalized frequency
phi=.25*pi;
noisevar=.05;
nosamples=1000;

% plot(n,A*cos(2*pi*fo*n + phi),n,x); % Original and noisy observations
for l=1:length(N),
n=0:N(l)-1;
for k=1:nosamples,
x= A*cos(2*pi*fo*n + phi)+(sqrt(noisevar)*randn(1,N(l)));
Phi_e(l,k)=-(atan2((x*(sin(2*pi*fo*n'))),(x*(cos(2*pi*fo*n')))));

end
end

plot(N,mean(Phi_e'),N,ones(1,length(N))*phi); % variances can alos be tabulated
xlabel('Value of N');ylabel('Mean of the Phase estimator');legend('Aymptotic','Actual');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second part of the experiment (Fixed data recored and vary the SNR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

N=80; %
n=0:N-1;
A=1; % knowm Amplitude
fo=.08; % Normalized frequency
phi=.25*pi; % Unknown but constant phase
s=A*cos(2*pi*fo*n + phi); 

nosamples=1000;
SNR=-20:3:10; % SNR in dB
SNR_R=10.^(SNR/10);
noisevar=((A^2)./(2*SNR_R));
var_phi=(1./(N*SNR_R));

% plot(n,A*cos(2*pi*fo*n + phi),n,x); % Original and noisy observations
for l=1:length(noisevar),
for k=1:nosamples,
x= s+(sqrt(noisevar(l))*randn(1,N));
Phi_e(k,l)=-(atan2((x*(sin(2*pi*fo*n'))),(x*(cos(2*pi*fo*n')))));
end
end

figure,plot(SNR,mean(Phi_e),SNR,ones(1,length(SNR))*phi,'k-.');
title('Actual Vs Asymptotic mean for Phase estimator');
xlabel(' SNR (dB) ');ylabel('Mean');legend('Aymptotic','Actual')

% figure,semilogy(SNR,var(Phi_e),SNR,ones(1,length(SNR)).*(1./(N*SNR_R)),'k-.');
figure,plot(SNR,10*log10(var(Phi_e)),SNR,10*log10(ones(1,length(SNR)).*(1./(N*SNR_R))),'k-.');
title('Actual Vs Asymptotic variance for Phase estimator');
xlabel(' SNR (dB) ');ylabel(' Varaiance  (dB) ');legend('Aymptotic','Actual')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% third part of the experiment (Plot of Log- Liklihood for Different SNR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

N=80; %
n=0:N-1;
A=1; % knowm Amplitude
fo=.08; % Normalized frequency
phi=-pi:.01*pi:pi; % Unknown but constant phase
true_phi=.25*pi;

notimes=3; % No of realizations of the Log likihood we want
SNR=-15; % SNR in dB
SNR_R=10.^(SNR/10);
noisevar=((A^2)./(2*SNR_R));
% var_phi=(1./(N*SNR_R));

% plot(n,A*cos(2*pi*fo*n + phi),n,x); % Original and noisy observations
s= A*cos(2*pi*fo*n +true_phi);

for k=1:notimes,
x= s +(sqrt(noisevar)*randn(1,N));
for l=1:length(phi),
% LLF(l,k)= (-N/2)*log(2*pi*noisevar)-(1/(2*noisevar))*(sum((x-(A*cos(2*pi*fo*n +phi(l)))).^2));
LLFF(l,k)= -(1/(2*noisevar))*(sum((x-(A*cos(2*pi*fo*n +phi(l)))).^2));
end
end

figure,plot(phi,LLFF)%
title('At Low SNR ( -15 dB )');
xlabel(' Phase,  {\phi} ');ylabel('Log-liklihood function');
axis([-pi pi -100 0])