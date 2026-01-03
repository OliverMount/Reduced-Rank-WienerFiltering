
function [H_capt,pickme,pickme1]=KTFt_esti_comb(y1,y2,Pilot,temp,Poss,M,D,p_gap)
% This function estimated TF in time domain
% function [H_capt,H_capt_RR]=KTFt_esti_comb(y,y2,Pilot,temp,Poss,M,D)
N=length(Pilot);
F=sqrt(1/N)*dftmtx(N);
Pilot_pos=Poss;
QQ=F'; QQQ=QQ(:,1:M); QQQ1=QQ(Pilot_pos,1:M);QQQQ=QQ(Pilot_pos,Pilot_pos);
a=F'*diag(Pilot)*QQQ*sqrt(N); % a in the data model

h_capt=(inv(a'*a))*a'*y1; % LS estimator for CIR (Full rank estimator)
% disp(cond(a'*a));
H_capt=QQQ*h_capt*sqrt(N); % TF estimator (Full rank)

% A=a'*a;b=a'*y;
aa=QQQQ*diag(temp)*QQQ1*sqrt(N);
h_capt1=(inv(aa'*aa))*aa'*y2;
H_capt1=QQQ1*h_capt1*sqrt(N);

pickme_t=interp(H_capt1,p_gap); % Lowpass interpolated for full rank
pickme=pickme_t(1:N); 

A=aa'*aa ; b=aa'*y2;
Kh_capt=LANCZOS_LS(D,M,A,b); % CIR estimate (Reduced rank)

for z=1:D,
H_capt_RR(:,z)=QQQ1*Kh_capt(:,z)*sqrt(N); % TF estimator (Reduced rank)
temp1=interp(H_capt_RR(:,z),p_gap); % Lowpass interpolated for reduced rank
pickme1(:,z)=temp1(1:N); 
end


