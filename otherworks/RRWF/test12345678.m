% This program for MMSE RR Equlaization 
clc;clear all;close all;
tic;
noiter=25; % number of averaing iterations
len=1000; % length of the transmiited sequence

del=7;  % Decision delay
M=11;  
D=8; % Dimension of the Subspace
p=M; % Full rank Dimension

Eb=1; % Energy per bit in V
EbNodB=0:15;
noisevar=1./10.^(EbNodB/10);

% h=[1 0.4]; % Assumed ISI Channel model
W=2.9;
for n=1:3,
h(n)=.5*(1+cos(2*pi*(n-2)/W)); % RC channel
end
L=length(h)-1; % ISI length
% figure,freqz(h);figure,zplane(h)

h_corr=corr(h,h)';
h_ori=[h_corr(L+1:end) zeros(1,M)]; % To eliminate the negative time lag 

for loop=1:noiter,   
    
    seq=randint(1,len);
    seq=2*seq-1; % Transmitted sequence

    r_dx_mtx=convolmtx(h,M); % Cross correlation vector (DOES NOT DEPEND ON NOISE)
    r_dx=r_dx_mtx(del+1,:);   % Change here the decision variable

        for i=1:length(noisevar),
            
            N=noisevar(i)*eye(M);  % Noise covariance matrix
            X=var(seq)*toeplitz(h_ori(1:M));  % Unit variance transmitted correlation
            Y=X+N; % Received Correlation matrix
            MMSE_w=inv(Y)*r_dx';  % Optimum MMSE Full rank solution
            
            y=convol(seq,h)';      % Channel output 
            rr=y+sqrt(noisevar(i)/2)*randn(1,length(y)); % RXed signal
            conv_Y=convolmtx(rr,M);
            s=convol(rr,MMSE_w')';  % MMSE filter 
            rx = s > 0;  % Decision Box 
            rx1=(2*rx-1);
            rx_all=rx1(del+1:del+len); % L= ISI length
            number=(seq==rx_all); %Counting the errors. (Change if del is an array) 
            pe(loop,i)=(len-sum(number))/len;
            MSE(loop,i)=h_ori(1)-MMSE_w'*X*MMSE_w;
%             SNR(loop,i)=10*log10((conj(MMSE_w')*X*MMSE_w)/(conj(MMSE_w')*N*MMSE_w)); % For equalized
%             SNRU(loop,i)=10*log10((var(seq)*h_ori(1)/(noisevar(i)))); % For Unequalized

            %__________________________________________________________________________
            % Laczos based MSNWF
            %__________________________________________________________________________

            Cfirst=zeros(D);
            Clast=zeros(D);
            % Initilaization
            T(:,1)=zeros(p,1);
            no=norm(r_dx);
            T(:,2)=(r_dx')/no;
            r(1,2)=0;r(2,2)=T(:,2)'*Y*T(:,2);b(2)=r(2,2);
            Cfirst(2,2)=(1/r(2,2));Clast(2,2)=(1/r(2,2));
            MSE_RR(1)=h_ori(1)-(no).^2*Cfirst(2,2);

            for z=3:D+1,
            u=Y*T(:,z-1)-r(z-1,z-1)*T(:,z-1)-r(z-2,z-1)*T(:,z-2);
            r(z-1,z)=norm(u);
            T(:,z)=u/r(z-1,z);
            r(z,z)=T(:,z)'*Y*T(:,z);
            b(z)=r(z,z)-r(z-1,z).^2*(1/b(z-1));
            Cfirst(2:z,z)=[Cfirst(2:z-1,z-1) ; 0] + (1/b(z))*Clast(2,z-1)*[r(z-1,z).^2*Clast(2:z-1,z-1) ; -r(z-1,z)];
            Clast(2:z,z)=(1/b(z))*[-r(z-1,z)*Clast(2:z-1,z-1);1] ;
            MSE_RR(z-1)=h_ori(1)-(no^2)*Cfirst(2,z);
            end
            V=T(:,2:D+1);
            % Obtained the Subspace(V) ; Now Do the RR Process

            C_first=Cfirst(2:end,2:end);
            RR_wopt=C_first*no;

            for m=1:D,
            convolmtx_d=conv_Y*V(:,1:m);
            estimated_RR(m,:)=(convolmtx_d*RR_wopt(1:m,m))';
            end
            
            rx_RR = estimated_RR(D,:) > 0;  % Decision Box 
            rx1=(2*rx_RR-1);
            rx_all=rx1(del+1:del+len); % L= ISI length
           
            number=(seq==rx_all); %Counting the errors. (Change if del is an array) 
            pe_RR(loop,i)=(len-sum(number))/len;
            clear V T;
            % _________________________________________________________________________
            % PC method of Rank reduction
            %__________________________________________________________________________
            
            temp=r_dx';
            I(1,1)=0;  % To store the Eigen Values
            [VV DD]=eig(Y);

                for Z=1:p,
                    TT=VV(:,1:Z); % eigen vectors
                    lamda(Z)=DD(Z,Z); % Eigen values
                    I(Z,Z)=lamda(Z);
                    Trans_auto=TT'*Y*TT;
                    Trans_cross=TT'*temp;
                    eig_wopt{1,Z}=inv(Trans_auto)*Trans_cross; % Optimum solution
                    conv_d=conv_Y*TT;
                    estimated_signal(Z,:)=(conv_d*eig_wopt{1,Z})';
%                     MSE_eig(Z)=h_ori(1)-temp'*TT*inv(I)*TT'*temp;
                end
            rx_eig = estimated_signal(D,:) > 0;  % Decision Box 
            rx1=(2*rx_eig-1);
            rx_all=rx1(del+1:del+len); % L= ISI length
           
            number=(seq==rx_all); %Counting the errors. (Change if del is an array) 
            pe_eig(loop,i)=(len-sum(number))/len; 
                              
                        
            %__________________________________________________________________________
             % Cross Spectral (CS) Metric method
            %__________________________________________________________________________

            for o=1:p,

            D_temp=zeros(o,o);
            NC=nchoosek(p,o);
            combi=nchoosek(1:p,o);

                for l=1:NC,
                    for q=1:o,
                    V_temp(:,q)=VV(:,combi(l,q));
                    D_temp(q,q)=DD(combi(l,q),combi(l,q));
                    end
                    CS(l)=temp'*V_temp*inv(D_temp)*V_temp'*temp; % CS Metric
                end
            clear l;
            [val pos]=max(CS);
            eig_pos=combi(pos,:);

                for m=1:o,
                eig_vectors_CS(:,m)=VV(:,eig_pos(1,m)); % Choosen Vectors which maximizes CS Metric
                lamda_CS(m,m)=DD(eig_pos(1,m),eig_pos(1,m));
                end
% disp(D);
% disp(lamda_CS);
% disp(MSE);disp(MSE_eig);

% Do RR now

                T=eig_vectors_CS;
                Trans_auto=T'*Y*T;
                Trans_cross=T'*temp;
                eig_wopt=inv(Trans_auto)*Trans_cross; % Optimum solution
                conv_d=conv_Y*T;
                estimated_signal_CS(o,:)=(conv_d*eig_wopt)';
                  
                 MSE_CS(o)=h_ori(1)-temp'*T*inv(lamda_CS)*T'*temp;
                clear V_temp D_temp q combi CS eig_pos eig_vectors_CS lamda_CS T Trans_auto Trans_cross eig_wopt conv_d;
                
            end
            
            rx_CS = estimated_signal_CS(D,:) > 0;  % Decision Box 
            rx1=(2*rx_CS-1);
            rx_all=rx1(del+1:del+len); % L= ISI length
           
            number=(seq==rx_all); %Counting the errors. (Change if del is an array) 
            pe_CS(loop,i)=(len-sum(number))/len; 
            
%             clear V T VV TT DD;
%                        
        end
    fprintf('Finished %d Iteration for %d bits  in %f seconds  \n',loop,len,toc);
end

No=noisevar;
EBNO=10*log10(Eb./No);
ensemble_avg=mean(pe);
ensemble_avg_RR=mean(pe_RR);
ensemble_avg_eig=mean(pe_eig);
ensemble_avg_CS=mean(pe_CS);
% avgSNR=mean(SNR);
% avgSNRU=mean(SNRU);

% ImprovedSNR=10*log10(

figure,hnd=semilogy(EBNO,.5*erfc(sqrt((Eb./No))),'b',EBNO,ensemble_avg,'rs-',EBNO,ensemble_avg_RR,'kp-',EBNO,ensemble_avg_eig,'k-.',EBNO,ensemble_avg_CS,'b:');grid;
set(hnd(1),'linewidth',2);
set(hnd(2),'linewidth',2);
set(hnd(3),'linewidth',2);
set(hnd(4),'linewidth',2);
set(hnd(5),'linewidth',2);
str1=num2str(D);
str=strcat('MSNWF D = ',str1);
legend('NO ISI','MMSE Equalized',str,'PC','CS');
axis([0 11 10e-8 0.5]);

% hnd=title('Equlaizer Performance Chart');
% set(hnd,'fontsize',15);
hnd=xlabel(' Eb / No   in   dB ');
set(hnd,'fontsize',13,'color',[0 0 1]);
hnd=ylabel('  Probability  of  Error (Pe)  ');
set(hnd,'fontsize',13,'color',[0 0 1]);

% SNR improvement chart

% figure,hnd=semilogy(EBNO,avgSNRU,'b',EBNO,avgSNR,'ro-');grid;
% set(hnd(1),'linewidth',2);
% set(hnd(2),'linewidth',2);
% legend('Un Equalized','MMSE Equalized');% This program for MMSE RR Equlaization 
