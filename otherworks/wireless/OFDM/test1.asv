% OFDM LS channel estimation
% 16 (3.5.) 32  64-7  512-7  256 -3
clear all;clc;close all;
N=512; % No of sub-carriers 
p_gap=7;
Pilot_pos=1:p_gap:N;
L=length(Pilot_pos); % For Comb type carriers
Z=4; % Signifies QPSK

% L=N; % for Block type pilots

F=sqrt(1/N)*dftmtx(N);
No=1;
% h=[1+.34j 1-.9j 1+.5j];
h=[1+j -.8+j .25-j];
% h=h/norm(h);
H=CircMtx(h,N); % Channel matrix for AWGN
M=length(h);
% [EE,E]=eig(H);

% Q_temp1=F'*F;
% Q=F(1:M,Pilot_pos); % Q for Comb-type Pilots
% Q_t=Q';
QQ=F'; % Q for Block type pilot carriers
% QQQ1=QQ(fliplr(Pilot_pos),1:M);
QQQ1=QQ(Pilot_pos,1:M);
QQQ=QQ(:,1:M);
QQQQ=QQ(Pilot_pos,Pilot_pos);

% QQQQ=sqrt(1/L)*dftmtx(L)';
temp=(2*randint(1,L)-1).';
Pilot=(upsample(temp,p_gap));
Pilot=flipud(Pilot(1:N));
% Pilot1=(2*randint(1,L)-1).';

% Channel Estimation part
noi=randn(1,N)+j*randn(1,N);
noise=sqrt(No/2)*noi.';
x=F'*Pilot; % taking IFFT before transmission
y=H*x+noise; % Passing thro Channel  + noise (received symbols)
Y=F*y; % FFT after reception

% [LE,he]=channel_esti_block(Pilot,Y,N,M); %Freq domain channel estimator
% LE=F'*[he; zeros(N-M,1)];
 
% % Frequency domain channel estimation (using data Y)
% h_capf=inv(Q*Q')*Q*(Pilot.*Y)/sqrt(N) % Freq domain CIR estimate
% b=diag(Pilot)*QQQ*sqrt(N);
% h_capf=inv(b'*b)*b'*Y; % Freq domain CIR estimate
% h.';
% H_cap=F'*[h_cap; zeros(N-M,1)]*sqrt(N)
H_ori=diag(F*H*F');
% % H_capf=QQQ*h_capf*sqrt(N); % Freq domain TF estiamte ( method  1 )
H_capf=(Pilot.*Y); % Freq domain TF estiamte ( method 2 )
% 
% 
% % H_p=[H_ori H_capf H_capf1];
% % figure,plot(abs(H_p));
% % legend('Original','Method 1','Method 2');
% 
% Time domain channel estimation (Using data y)
y1=F'*diag(Pilot)*QQQ*h.'*sqrt(N)+noise;
a=F'*diag(Pilot)*QQQ*sqrt(N);
% a=F'*diag(Pilot)*QQQ*sqrt(N);
% a=QQQQ*diag(flipud(temp))*QQQ*sqrt(N);
% y_temp=flipud(y);
% y_pilot=flipud(y_temp(Pilot_pos));
% y_pilot_c=flipud(y(N:-p_gap:1));
% h_capt=(inv(a'*a))*a'*y_pilot
h_capt=(inv(a'*a))*a'*y;
% P_F=inv(QQQ'*QQQ)*QQQ';
% H_capt=(inv(P_F'*a'*a*P_F))*P_F'*a'*y/sqrt(N);
H_capt=QQQ*h_capt*sqrt(N);
% y_p=abs([y y1]);
% figure,plot(y_p);
% legend('y','ymodel');
% 
Y1=diag(Pilot)*QQQ*h.'*sqrt(N);%+F*noise;
% figure,plot(abs([Y Y1]));
% legend('Y','Ymodel');
% h.';
% H_p=[H_ori H_capf H_capt];
% figure,plot(abs(H_p));
% legend('Original','H cap via Freq domain','H via Time domain ');
% figure,plot(angle(H_p));
% legend('Original','H cap via Freq domain','H via Time domain ');
noi_temp=F*noise;
noi_temp1=noise;
% y2=QQQQ*diag(flipud(temp))*QQQ1*h.'*sqrt(N)+QQQQ*noi_temp(flipud(Pilot_pos));
% bb=diag(Pilot)*QQQ*sqrt(N);
% Y1=diag(Pilot)*QQQ1*h.'*sqrt(N)+F*noise;
Y2=diag(flipud(temp))*QQQ1*h.'*sqrt(N)+noi_temp(flipud(Pilot_pos));
% y2=QQQQ*diag(flipud(temp))*QQQ1*h.'*sqrt(N)+noi_temp1(flipud(Pilot_pos));
y2=flipud(y(N:-p_gap:1))
y_temp=flipud(y(N:-p_gap:1))
aa=QQQQ*diag(flipud(temp))*QQQ1*sqrt(N);
h_capt1=(inv(aa'*aa))*aa'*y_temp;
H_capt1=QQQ1*h_capt1*sqrt(N);

% Y3=diag(flipud(temp))*Q_t*h.'*sqrt(N)
Y_temp=flipud(Y(N:-p_gap:1));

% H_p=[H_ori H_capt]
% figure,plot(abs(H_p));
% figure,plot(angle(H_p));
% legend('Original','Block type');
% % legend('Original','Block type','Comb-type without interpolation');

H_temp1=flipud(H_capt(N:-p_gap:1));
H_temp=flipud(H_ori(N:-p_gap:1));

H_q=[H_temp H_temp1 H_capt1];
% 
% figure,plot(abs([H_q]));legend('Ori','ori1','esti');
% figure,plot(angle([H_q]));legend('Ori','ori1','esti');

 pickme_t=interp(H_capt1,p_gap);
 pickme=pickme_t(1:N);
 
H_p=[H_ori H_capt pickme]
% figure,plot(abs([H_p]));legend('Ori','ori1','esti');
% figure,plot(angle([H_p]));legend('Ori','ori1','esti');
nnn=flipud(H_capf);
f_H=nnn(Pilot_pos);
pickme_f=interp(f_H,p_gap);
 pickme_ff=pickme_f(1:N);
 
 H_pp=[H_ori flipud(pickme_ff)];
 figure,plot(abs([H_pp]));legend('Ori','esti');
 