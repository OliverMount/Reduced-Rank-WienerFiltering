function H_capt=TFt_esti(y,Pilot,M)

% This function estimated TF in time domain
N=length(y);
F=sqrt(1/N)*dftmtx(N);

Q=F'; % Q for Block type pilot carriers
QQ=Q(:,1:M);
a=F'*diag(Pilot)*QQ*sqrt(N); % a in the data model
h_capt=(inv(a'*a))*a'*y; % LS estimator for CIR
H_capt=QQ*h_capt*sqrt(N); % TF estimator
