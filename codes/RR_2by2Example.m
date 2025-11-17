clear all;clc;close all;
R=[10 4 ;4 10]; p=[9 1]';sigma_d=10; % Example from RR book
% R=[1.1 .5 ;.5 1.1]; p=[.5272 -.4458]';sigma_d=.9486; % From simon Haykin
% R=[.5 .20 ;.20 .5]; p=[0 -.866025]';sigma_d=2; % From Widrow sterns
% R=[2 1 ;1 2]; p=[7 8]';sigma_d=42; % Another critical from Widrow
[V,D]=eig(R);

% err_full=10*log10(sigma_d- p'*inv(R)*p); in dB's
err_full=(sigma_d- p'*inv(R)*p);
for m=1:2,
    temp(m)=p'*V(:,m)*(1/D(m,m))*V(:,m)'*p;
%     er(m)=10*log10(sigma_d-temp(m)); % in dB's
er(m)=(sigma_d-temp(m));
end

% T=[9 1; -1 9]*(1/sqrt(82)); % T for example 1
% T=[p';0.4458 0.5272]*(1/norm(p));  T for example 2
T=[0 -.866025; .866025 0]/norm(p);

Rd=T(1,:)*R*T(1,:)';
% disp(10*log10(sigma_d-p'*inv(T'*Rd*T)*p)); in dB's
disp((sigma_d-p'*inv(T'*Rd*T)*p)); 

ori_p=[1 0]';
% Performance plotting
thet=(0:359)*(pi/180);
for ha=1:360,
RRR=[cos(thet(ha)) -sin(thet(ha));sin(thet(ha)) cos(thet(ha))];
P(:,ha)=RRR*ori_p;

end

% Results for example 1
% err_full= 1.0952 ;
% er= 4.6667    6.4286;
% er_RR= 2.4619;

% Results for Simon haykin
% err_full =  0.1576
% er = 0.1597    0.9465
% er_RR=0.1632

% Results from Widrow sterns