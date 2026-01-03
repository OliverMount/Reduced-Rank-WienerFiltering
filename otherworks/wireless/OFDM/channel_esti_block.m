function [LE,he]=channel_esti_block(Pilot,Y,N,M)

L=N; % for Block type pilots
F=(1/sqrt(N))*dftmtx(N);
Q=F(1:M,1:N); % Q for Block type pilot carriers

 % Channel Estimation part
he=(Q*(Y./Pilot))/sqrt(N);
LE=F*(F'*[he; zeros(N-M,1)]);

