% OFDM LS channel estimation

clear all;clc;close all;
N=16; % No of sub-carriers
p_gap=5;
Pilot_pos=1:p_gap:N;
L=length(Pilot_pos); % For Comb type carriers
Z=4; % Signifies QPSK

% L=N; % for Block type pilots

F=sqrt(1/N)*dftmtx(N);

No=1;
% h=[1+.34j 1-.9j 1+.5j];
h=[1+j .9+j .25-j];
% h=h/norm(h);
H=CircMtx(h,N); % Channel matrix for AWGN
M=length(h);
% [EE,E]=eig(H);

% Q=F(1:M,1:N); % Q for Comb-type Pilots
QQ=F'; % Q for Block type pilot carriers
QQQ1=QQ(Pilot_pos,1:M);
QQQ=QQ(:,1:M);
QQQQ=QQ(Pilot_pos,Pilot_pos);
temp=(2*randint(1,L)-1).';
Pilot=(upsample(temp,p_gap));
Pilot=flipud(Pilot(1:N));


% Channel Estimation part
noi=randn(1,N)+j*randn(1,N);
noise=sqrt(No/2)*noi.';
x=F*Pilot; % taking IFFT before transmission
y=H*x+noise; % Passing thro Channel  + noise (received symbols)
% Y=F*y; % FFT after reception

% [LE,he]=channel_esti_block(Pilot,Y,N,M); %Freq domain channel estimator
% LE=F'*[he; zeros(N-M,1)];

 
% % Frequency domain channel estimation (using data Y)
% % h_capf=inv(Q*Q')*Q*(Pilot.*Y)/sqrt(N) % Freq domain CIR estimate
% b=diag(Pilot)*QQQ*sqrt(N);
% h_capf=inv(b'*b)*b'*Y; % Freq domain CIR estimate
% % h.';
% % H_cap=F'*[h_cap; zeros(N-M,1)]*sqrt(N)
% H_ori=diag(F*H*F');
% % H_capf=QQQ*h_capf*sqrt(N); % Freq domain TF estiamte ( method  1 )
% H_capf=(Pilot.*Y); % Freq domain TF estiamte ( method 2 )
% 
% 
% % H_p=[H_ori H_capf H_capf1];
% % figure,plot(abs(H_p));
% % legend('Original','Method 1','Method 2');
% 
% Time domain channel estimation (Using data y)
% y1=F'*diag(Pilot)*QQQ*h.'*sqrt(N)+noise;
% a=F'*diag(flipud(Pilot))*QQQ*sqrt(N);
% % a=F'*diag(Pilot)*QQQ*sqrt(N);
% % a=QQQQ*diag(flipud(temp))*QQQ*sqrt(N);
% % y_temp=flipud(y);
% % y_pilot=flipud(y_temp(Pilot_pos));
% % y_pilot_c=flipud(y(N:-p_gap:1));
% % h_capt=(inv(a'*a))*a'*y_pilot
% h_capt=(inv(a'*a))*a'*y;
% % P_F=inv(QQQ'*QQQ)*QQQ';
% % H_capt=(inv(P_F'*a'*a*P_F))*P_F'*a'*y/sqrt(N);
% H_capt=QQQ*h_capt*sqrt(N);
% y_p=abs([y y1]);
% figure,plot(y_p);
% legend('y','ymodel');
% 
% Y1=diag(Pilot)*QQQ*h.'*sqrt(N)+F*noise;
% figure,plot(abs([Y Y1]));
% legend('Y','Ymodel');
% h.';
% H_p=[H_ori H_capf H_capt];
% figure,plot(abs(H_p));
% legend('Original','H cap via Freq domain','H via Time domain ');
% figure,plot(angle(H_p));
% legend('Original','H cap via Freq domain','H via Time domain ');