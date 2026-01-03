% Adpative noise reduction
%________________________________________________________________________
% Conventional Wiener solution
%________________________________________________________________________
clc;clear all;close all;
w0=0.05*pi;n=0:199;;phi=0;
d=sin(w0*n+phi); % Desired signal
k=length(d);

% Noise source generator
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

Rv2=toeplitz(rv2(1:p)); % Formation of the autocorrelation matrix
w_opt=inv(Rv2)*rxv2(1:p)'; % finding the optimal solution

conv_rv2=convolmtx(v2,p);
estimated_noise=convol(v2,w_opt')';
% figure,plot(n,v1,n,estimated_noise(1:k));
% legend('Original','Estimated');
% figure,plot(v1);
recovered_signal=x-estimated_noise(1:k);
% figure,plot(n,d,'r',n,recovered_signal);
% legend('Original','Estimated');
MSE=rv1(1)-w_opt'*Rv2*w_opt;

for del=0:p-1,
    err(del+1)=sum((x-estimated_noise(1+del:k+del)).^2);
end
% disp(err);
% _________________________________________________________________________
% PC method of Rank reduction
%__________________________________________________________________________

temp=rxv2(1:p)';
I(1,1)=0;  % To store the Eigen Values
[V D]=eig(Rv2);
disp(diag(D));

for i=p:-1:1,
    T=V(:,p:-1:i); % eigen vectors
    lamda(p-i+1)=D(i,i); % Eigen values
    I(p-i+1,p-i+1)=lamda(p-i+1);
    Trans_auto=T'*Rv2*T;
    Trans_cross=T'*temp;
    eig_wopt{1,i}=inv(Trans_auto)*Trans_cross; % Optimum solution
    conv_d=conv_rv2*T;
    estimated_signal(p-i+1,:)=(conv_d*eig_wopt{1,i})';
    desired_signal_eigen(p-i+1,:)=x-estimated_signal(p-i+1,1:k);
    MSE_eig(p-i+1)=rv1(1)-temp'*T*inv(I)*T'*temp;
end
% figure,plot(1:200,desired_signal_eigen(6,:)',1:k,recovered_signal);
% figure,plot(1:200,recovered_signal);
disp(T);

%__________________________________________________________________________
% Cross Spectral (CS) Metric method
%__________________________________________________________________________

for o=1:p,

D_temp=zeros(o,o);
M=nchoosek(p,o);
combi=nchoosek(1:p,o);

for l=1:M,
    for q=1:o,
        V_temp(:,q)=V(:,combi(l,q));
        D_temp(q,q)=D(combi(l,q),combi(l,q));
    end
        CS(l)=temp'*V_temp*inv(D_temp)*V_temp'*temp; % CS Metric
    end
    clear l;
    [val pos]=max(CS);
    eig_pos=combi(pos,:);

    for m=1:o,
        eig_vectors_CS(:,m)=V(:,eig_pos(1,m)); % Choosen Vectors which maximizes CS Metric
        lamda_CS(m,m)=D(eig_pos(1,m),eig_pos(1,m));
    end
% disp(D);
% disp(lamda_CS);
% disp(MSE);disp(MSE_eig);

% Do RR now

    T=eig_vectors_CS;
    Trans_auto=T'*Rv2*T;
    Trans_cross=T'*temp;
    eig_wopt=inv(Trans_auto)*Trans_cross; % Optimum solution
    conv_d=conv_rv2*T;
    estimated_signal_CS=(conv_d*eig_wopt)';
    desired_signal_CS(o,:)=x-estimated_signal_CS(1:k);
    
    MSE_CS(o)=rv1(1)-temp'*T*inv(lamda_CS)*T'*temp;
    clear V_temp D_temp q combi CS eig_pos eig_vectors_CS lamda_CS T  Trans_cross eig_wopt conv_d;
end

% figure,plot(1:k,recovered_signal,1:k,desired_signal_eigen(6,:),1:k,desired_signal_CS(6,:)); Trans_auto
% legend('Full rank Wiener','Eigen Decomposition (6)','CS (6)'); 
% 
% % Original ,ED ,CS
% figure,plot(1:k,d,1:k,desired_signal_eigen(6,:),1:k,desired_signal_CS(6,:));
% legend('Original','Eigen Decomposition (6)','CS (6)'); 
% Original ,FRW ,CS
% figure,plot(1:k,d,1:k,recovered_signal,1:k,desired_signal_CS(6,:));
% legend('Original','Full rank Wiener','CS (6)'); 
fprintf('Full rank MSE %f \n',MSE);
fprintf('Eigen MSE \n');disp(MSE_eig);
fprintf('CS Metric MSE \n');disp(MSE_CS);
 
%__________________________________________________________________________
% Laczos based MSNWF
%__________________________________________________________________________
clear D V T;
D=12;
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
V=T(:,2:D+1);
end

fprintf('MSNWF MSE \n');disp(MSE_RR);
% Obtained the Subspace(V) ; Now Do the RR Process

C_first=Cfirst(2:end,2:end);
RR_wopt=C_first*no;

for m=1:D,
    convolmtx_d=conv_rv2*V(:,1:m);
    estimated_RR(m,:)=(convolmtx_d*RR_wopt(1:m,m))';
    recovered_signal_RR(m,:)=x-estimated_RR(m,1:k);
end

disp(T);
figure,plot(1:k,recovered_signal,'r',1:k,recovered_signal_RR(6,:),'k*',1:k,desired_signal_CS(6,:),'b*');
% axis([-1 205 -2 2]);
hnd=title('Full rank Wiener CS and MSNWF ( r = 6 ) Comparison');
set(hnd,'fontsize',15);
hnd=xlabel('Sample Index');set(hnd,'fontsize',15);
hnd=ylabel('Amplitude');set(hnd,'fontsize',15);

str1=num2str(MSE); str2=num2str(MSE_RR(6));str3=num2str(MSE_CS(6));
legend(str1,str2,str3);

%__________________________________________________________________________
% MSE as a function of RANK
%__________________________________________________________________________

MSE_fullrank=MSE*ones(1,p);
figure,hnd=plot(1:p,MSE_fullrank,1:p,MSE_eig,1:p,MSE_CS,1:p,MSE_RR,'kp');
set(hnd(1),'linewidth',2);
axis([ 0 13 -.5 3]) 
hnd=title('MSE as a Function of RANK');
set(hnd,'fontsize',15);
hnd=xlabel('Rank');set(hnd,'fontsize',15);
hnd=ylabel('MSE');set(hnd,'fontsize',15);
legend('Full Rank Wiener','PC','CS','MSNWF');