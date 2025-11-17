function V=LANCZOS(D,N,A,b,sig_d)

% LANCZOS a method for solving Ax=b in Kyrlov subspace
Cfirst=zeros(D);
Clast=zeros(D); p=N;
% Initialization
temp=norm(b);
T(:,1)=zeros(p,1);
T(:,2)=(b)/norm(b);
r(1,2)=0;r(2,2)=T(:,2)'*A*T(:,2);b(2)=r(2,2);
Cfirst(2,2)=(1/r(2,2));Clast(2,2)=(1/r(2,2));
% MSE_RR(1)=sig_d-(norm(b)).^2*Cfirst(2,2);

for i=3:D+1,
    u=A*T(:,i-1)-r(i-1,i-1)*T(:,i-1)-r(i-2,i-1)*T(:,i-2);
    r(i-1,i)=norm(u);
    T(:,i)=u/r(i-1,i);
    r(i,i)=T(:,i)'*A*T(:,i);
    b(i)=r(i,i)-r(i-1,i).^2*(1/b(i-1));
    Cfirst(2:i,i)=[Cfirst(2:i-1,i-1) ; 0] + (1/b(i))*Clast(2,i-1)*[r(i-1,i).^2*Clast(2:i-1,i-1) ; -r(i-1,i)];
    Clast(2:i,i)=(1/b(i))*[-r(i-1,i)*Clast(2:i-1,i-1);1] ;
%     MSE_RR(i-1)=sig_d-(norm(b)^2)*Cfirst(2,i);
end
V=T(:,2:D+1);

C_first1=Cfirst(2:end,2:end)

RR_wopt1=V*(C_first1*temp);


