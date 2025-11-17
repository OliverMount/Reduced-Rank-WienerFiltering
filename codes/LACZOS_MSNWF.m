% LANCZOS based MSNWF

%________________________________________________________________________
% Conventional Wiener solution
%________________________________________________________________________
clc;clear all;close all;
w0=0.05*pi;n=0:199;;phi=0;
d=sin(w0*n+phi); % Desired signal
k=length(d);
% figure,plot(d);

noi=randn(1,k);

v1=AutoRnoise([ 1 -.8]',noi,k); % Noise from Primary sensor
v2=AutoRnoise([1 .6]',noi,k); % Noise from secondary sensor
x=d+v1;
% figure,plot(n,x,n,d);

rv1=TimeAC(v1);
rv2=TimeAC(v2);
% figure,plot(rv2);

rxv2=TimeCC(x,v2);
% figure,plot(rxv2);

p=12; % Order of the noise estimator;
D=6; % Reduced rank

Rv2=toeplitz(rv2(1:p)); % Formation of the autocorrelation matrix
w_opt=inv(Rv2)*rxv2(1:p)'; % finding the optimal solution

conv_rv2=convolmtx(v2,12);
estimated_noise=convol(v2,w_opt')';
% figure,plot(n,v1,n,estimated_noise(1:k));
% legend('Original','Estimated');
% figure,plot(v1);
recovered_signal=x-estimated_noise(1:k);
% figure,plot(n,d,'r',n,recovered_signal);
% legend('Original','Estimated');
MSE=rv1(1)-w_opt'*Rv2*w_opt;

% for del=0:p-1,
%     err(del+1)=sum((x-estimated_noise(1+del:k+del)).^2);
% end

%__________________________________________________________________________
% Laczos based MSNWF
%__________________________________________________________________________

Cfirst=zeros(D);
Clast=zeros(D);
% Initilaization
T(:,1)=zeros(p,1);
no=norm(rxv2(1:p));
T(:,2)=(rxv2(1:p)')/no;
r(1,2)=0;r(2,2)=T(:,2)'*Rv2*T(:,2);b(2)=r(2,2);
Cfirst(2,2)=(1/r(2,2));Clast(2,2)=(1/r(2,2));
MSE_RR(1)=rv1(1)-(no).^2*Cfirst(2,2);

for i=3:D+1,
    u=Rv2*T(:,i-1)-r(i-1,i-1)*T(:,i-1)-r(i-2,i-1)*T(:,i-2);
    r(i-1,i)=norm(u);
    T(:,i)=u/r(i-1,i);
    r(i,i)=T(:,i)'*Rv2*T(:,i);
    b(i)=r(i,i)-r(i-1,i).^2*(1/b(i-1));
    Cfirst(2:i,i)=[Cfirst(2:i-1,i-1) ; 0] + (1/b(i))*Clast(2,i-1)*[r(i-1,i).^2*Clast(2:i-1,i-1) ; -r(i-1,i)];
    Clast(2:i,i)=(1/b(i))*[-r(i-1,i)*Clast(2:i-1,i-1);1] ;
    MSE_RR(i-1)=rv1(1)-(no^2)*Cfirst(2,i);
end
ZZ=T(:,2:D+1);
end

disp(MSE);disp(MSE_RR(6));
% Obtained the Subspace(V) ; Now Do the RR Process

C_first=Cfirst(2:end,2:end);
RR_wopt=C_first*no;

for m=1:D,
    convolmtx_d=conv_rv2*ZZ(:,1:m);
    estimated_RR(m,:)=(convolmtx_d*RR_wopt(1:m,m))';
    recovered_signal_RR(m,:)=x-estimated_RR(m,1:k);
end

figure,plot(1:k,recovered_signal,'r',1:k,recovered_signal_RR(6,:),'k*');
title('Full rank And MSNWF ( r=6 ) Comparison');
str1=num2str(MSE); str2=num2str(MSE_RR);
legend(str1,str2);

LANCZOS(7,12,Rv2,rxv2(1:p).',rv1(1));

