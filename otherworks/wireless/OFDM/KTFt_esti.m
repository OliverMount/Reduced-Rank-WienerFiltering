function [KH_capt,H_capt_RR]=KTFt_esti(y,Pilot,M,D)

% This function estimated TF in time domain
N=length(Pilot);
F=sqrt(1/N)*dftmtx(N);

Q=F'; % Q for Block type pilot carriers
QQ=Q(:,1:M);
a=F'*diag(Pilot)*QQ*sqrt(N); % a in the data model

h_capt=(inv(a'*a))*a'*y; % LS estimator for CIR (Full rank estimator)
KH_capt=QQ*h_capt*sqrt(N); % TF estimator (Full rank)

A=a'*a;b=a'*y;
Kh_capt=LANCZOS_LS(D,M,A,b); % CIR estimate (Reduced rank)
for z=1:D,
H_capt_RR(:,z)=QQ*Kh_capt(:,z)*sqrt(N); % TF estimator (Reduced rank)
end



