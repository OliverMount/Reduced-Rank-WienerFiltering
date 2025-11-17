clear all;clc;close all;
% R=[10 4 ;4 10]; sigma_d=10; % Example from RR book
% R=[1.1 .5 ;.5 1.1]; sigma_d=0.9486; % From simon Haykin p=[.5272 -.4458]';
% R=[.5 .25 ;.25 .5];sigma_d=2; % From Widrow sterns
R=[2 1 ;1 2]; p=[7 8]';sigma_d=42; % Another critical from Widrow
[V,D]=eig(R);
% ori_p=[.5272 -.4458]';
% ori_p=[9 1]';
ori_p=[1 0]';
% Performance plotting
thet=(0:359)*(pi/180);
for ha=1:360,
RRR=[cos(thet(ha)) -sin(thet(ha));sin(thet(ha)) cos(thet(ha))];
P=RRR*ori_p;
disp(P);
err_full(ha)=(sigma_d-P'*inv(R)*P); % full rank error for all P

for m=1:2,
    temp(m,ha)=P'*V(:,m)*(1/D(m,m))*V(:,m)'*P;  
    er(m,ha)=(sigma_d-temp(m,ha));  % PC and CSM error for all P
end

T=[P'; -P(2) P(1)]/norm(P);
Rd=T(1,:)*R*T(1,:)';
err_RR(ha)=sigma_d-P'*inv(T'*Rd*T)*P;
end

fi=[er ;err_full;err_RR];
hnd=plot(fi');legend('Small','Large','Full-rank','MSNWF');grid;
set(hnd(1),'linewidth',1.2);set(hnd(2),'linewidth',1.2);set(hnd(3),'linewidth',1.2);
set(hnd(4),'linewidth',1.2,'color',[0 0 0]);
hnd=xlabel('Cross correlation vector angles in degrees');set(hnd,'fontsize',13,'color',[0 0 1]);
hnd=ylabel('MSE');set(hnd,'fontsize',13,'color',[0 0 1]);
% axis([ 0 365 9.8 10.07]) 
% 
% % err_full=10*log10(sigma_d- p'*inv(R)*p); in dB's
% err_full=(sigma_d- p'*inv(R)*p);
% for m=1:2,
%     temp(m)=p'*V(:,m)*(1/D(m,m))*V(:,m)'*p;
% %     er(m)=10*log10(sigma_d-temp(m)); % in dB's
% er(m)=(sigma_d-temp(m));
% end
% 
% % T=[9 1; -1 9]*(1/sqrt(82)); % T for example 1
% % T=[p';0.4458 0.5272]*(1/norm(p));  T for example 2
% T=[0 -.866025; .866025 0]/norm(p);
% 
% Rd=T(1,:)*R*T(1,:)';
% % disp(10*log10(sigma_d-p'*inv(T'*Rd*T)*p)); in dB's
% disp((sigma_d-p'*inv(T'*Rd*T)*p)); 

% Results for example 1
% err_full= 1.0952 ;
% er= 4.6667    6.4286;
% er_RR= 2.4619;

% Results for Simon haykin
% err_full =
% 
%     0.1576
% er =
% 
%     0.1597    0.9465
% er_RR=0.1632

% Results from Widrow sterns