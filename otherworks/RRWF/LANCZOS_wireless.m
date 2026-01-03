% This program for MMSE equalization of Known channel (with delayed
% decision ) in a BPSK Transmission ( no unequlai comparison)

% RR_QPSK
clear all;clc;tic;
noiter=25; Z=4; len=1000;
startr=6;endr=4;
del=0;  % Decision delay (From del=5 started getting results) 
M=25;  % Number of Taps of the equalizer
% D=5; % Dimension of the Subspace
p=M; % Full rank Dimension
no_frames=100;
bits_frame=200; % for one frame channel is constant
Pow_delay=[-3 -5 -8 -10  -12 -11 -8 -10 -12 -9 -9.23]; % Power delay profile values in dB
% Pow_delay=[-10.3621 -9.3930 -6.3639 -8.9620  -9.3930  -11.3077  -13.3724 -11.3077 -12.9243 -14.9485 -17.4473 ];
Eb=1; EbNodB=0:22; No=1./10.^(EbNodB/10);  Es=log2(Z)*Eb;
codebook=sqrt(2)*exp(j*2*pi*(0:Z-1)/Z + j*pi/4);
len_no=length(No);
len_r=startr-endr+1;
% cell arrays for counting all the rank pe
pe_RR=cell(1,len_r);
for tee=1:len_r,
    pe_RR{1,tee}=zeros(noiter,len_no);
end

h_temp=CallJack_function(no_frames,11);
h=(diag(10.^(Pow_delay./10))*((h_temp(5,:).')/norm(h_temp(5,:)))).';
% h=(diag(10.^(Pow_delay./10))*((h_temp(1,:).'))).';
% full rank 25, !1-23 !2-24 3-24 !4-0 !5-0 !6-0 7-0 8-15 9-  10  11 12-1
% 13- 14-
L=length(h)-1; % ISI length
h_corr=corr(h,h).';
h_ori=[h_corr(L+1:end) zeros(1,M)]; % To eliminate the negative time lag 

for loop=1:noiter,
      
    seq=randint(1,len,[0 Z-1]); % Mapping of Bits            
    TRANSseq=sqrt(Es)*exp((j*2*pi*seq/Z)+(j*pi/4));
    
    r_dx_mtx=convolmtx(h,M); % Cross correlation vector (DOES NOT DEPEND ON NOISE)
    r_dx=r_dx_mtx(del+1,:);   % Change here the decision variable

    for i=1:length(No),
               
        N=No(i)*eye(M);  % Noise covariance matrix
        X=var(TRANSseq)*toeplitz(h_ori(1:M));  % Unit variance transmitted correlation
        Y=X+N; % Received Correlation matrix
        MMSE_w=inv(Y)*r_dx.';  % Optimum MMSE Full rank solution
        
        y=convol(TRANSseq,h).';      % Channel output 
        leny=length(y);
        noi=randn(1,leny)+j*randn(1,leny);
        noise=sqrt(No(i)/2)*noi;
        rr=y+noise; % RXed signal
        conv_Y=convolmtx(rr,M); % needed if RR is used
        s=convol(rr,MMSE_w.').';  % MMSE extimator
%         RECEIVEDseq=TRANSseq+noise;
        
    for k=1:length(s),
        [val pos]=min(abs(s(k)-codebook));
        decoded(k)=codebook(pos);
    end
    
    number=(decoded(del+1:del+len)==TRANSseq); %Counting the errors. 
    pe(loop,i)=(len-sum(number))/len;
    clear decoded;
    for D=startr:-1:endr,
    
%           %__________________________________________________________________________
            % Lanczos based MSNWF
            %__________________________________________________________________________

            Cfirst=zeros(D);
            Clast=zeros(D);
            % Initilaization
            T(:,1)=zeros(p,1);
            no=norm(r_dx);
            T(:,2)=(r_dx.')/no;
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
            % Signal detection for different ranks
            for m=1:D,
            convolmtx_d=conv_Y*conj(V(:,1:m));
            estimated_RR(m,:)=(convolmtx_d*conj(RR_wopt(1:m,m))).';
            clear convolmtx_d;
            end
                    
            for k=1:length(estimated_RR(D,:)),
                 [val pos]=min(abs(estimated_RR(D,k)-codebook));
                 decoded(D,k)=codebook(pos);
            end
           
            number=(decoded(D,del+1:del+len)==TRANSseq); %Counting the errors. (Change if del is an array) 
            pe_RR{1,D}(loop,i)=(len-sum(number))/len;
            clear V T; 
        end

end
%     fprintf('Finished %d Iteration for %d bits  in %f seconds  for %d ranks\n',loop,len,toc,len_r);
fprintf('Finished %d Iteration for %d symbols  in %f seconds\n',loop,len,toc);
end

EBNO=10*log10(Eb./No);
ensemble_avg=mean(pe/2);
for D=startr:-1:endr,
ensemble_avg_RR(D,:)=mean(pe_RR{1,D}/2);
end

% figure,hnd=semilogy(EBNO,erfc(sqrt((Eb./No))).*(1-.25*erfc(sqrt((Eb./No)))),'k',EBNO,ensemble_avg,'b-');grid;
% set(hnd(1),'linewidth',2);
% set(hnd(2),'linewidth',2);
% legend('Theoretical No-ISI','Practical');
% axis([0 20 10e-7 .5]);
% 
% hnd=title('Coherent QPSK Performance Chart');
% set(hnd,'fontsize',15);
% hnd=xlabel(' Eb / No   in   dB ');
% set(hnd,'fontsize',15);
% hnd=ylabel('  Probability  of  Symbol  Error  ');
% set(hnd,'fontsize',15);

%2*Q(sqrt(2*Eb./No)).*(1-.5*Q(sqrt(2*Eb./No)))

% figure,scatter(real(rr),imag(rr));
% figure,hnd=semilogy(EBNO,(erfc(sqrt((Eb./No))).*(1-.25*erfc(sqrt((Eb./No)))))/2,'k-',EBNO,ensemble_avg,'b-');
% figure,hnd=semilogy(EBNO,erfc(sqrt((Eb./No))).*(1-.25*erfc(sqrt((Eb./No)))),'k-',EBNO,ensemble_avg,'bp-',EBNO,ensemble_avg_RR(startr,:),'ko-',EBNO,ensemble_avg_RR(startr-1,:),'k*-',EBNO,ensemble_avg_RR(startr-2,:),'k>-');grid;
% figure,hnd=semilogy(EBNO,erfc(sqrt((Eb./No))).*(1-.25*erfc(sqrt((Eb./No)))),'k-',EBNO,ensemble_avg,'bp-',EBNO,ensemble_avg_RR(startr,:),'ko-',EBNO,ensemble_avg_RR(startr-1,:),'k*-',EBNO,ensemble_avg_RR(startr-2,:),'k>-',EBNO,ensemble_avg_RR(startr-3,:),'b>-');grid;

% USe this rank reduced upto four ranks
figure,hnd=semilogy(EBNO,ensemble_avg,'bp-',EBNO,ensemble_avg_RR(startr,:),'ko-',EBNO,ensemble_avg_RR(startr-1,:),'k*-',EBNO,ensemble_avg_RR(startr-2,:),'k>-',EBNO,ensemble_avg_RR(startr-3,:));grid;
axis([0 22 10e-5 1]);
legend('Full Rank','MSNWF D=6','MSNWF D=5','MSNWF D=4'); % use for 4 ranks
hnd=xlabel(' SNR  in   dB ');
set(hnd,'fontsize',12,'color',[0 0 0]);
hnd=ylabel('  Probability  of  Bit Error (Pe)  ');
set(hnd,'fontsize',12,'color',[0 0 0]);
% figure,hnd=semilogy(EBNO,(erfc(sqrt((Eb./No))).*(1-.25*erfc(sqrt((Eb./No)))))/2,'k-',EBNO,ensemble_avg,'bp-',EBNO,ensemble_avg_RR(startr,:),'ko-');%,EBNO,ensemble_avg_RR(startr-1,:),'k*-',EBNO,ensemble_avg_RR(startr-2,:),'k>-',EBNO,ensemble_avg_RR(startr-3,:));grid;
% set(hnd(1),'linewidth',1.5);
% set(hnd(2),'linewidth',1.5);

% set(hnd(3),'linewidth',1.5);
% set(hnd(4),'linewidth',1.5);
% set(hnd(5),'linewidth',1.5);
% set(hnd(6),'linewidth',1.5);
%set(hnd(7),'linewidth',1.5);
% axis([0 20 10e-5 1]);
% 
% legend('NO ISI','Full Rank','MSNWF D=4','MSNWF D=3','MSNWF D=2','MSNWF D=2'); % use for 4 ranks
% legend('NO ISI','Full Rank','MSNWF D=6','MSNWF D=5','MSNWF D=4'); % use for 4 ranks
% % axis([0 20 10e-5 1]);
% % 
% hnd=xlabel(' SNR  in   dB ');
% set(hnd,'fontsize',12,'color',[0 0 0]);
% hnd=ylabel('  Probability  of  Bit Error (Pe)  ');
% set(hnd,'fontsize',12,'color',[0 0 0]);

% program to choose the order of the filter
% desired_point=25;
% for newvar=1:desired_point, 
%     Rv2=toeplitz(Y(1,1:newvar)); 
%     w_opt=Levinson(Y(1,1:newvar),r_dx(1:newvar)); % finding the optimal solution
%     MSE_new(newvar)=var(seq)-w_opt'*Rv2*w_opt;
% end
% 
% figure,plot(abs(MSE_new));
% *************************************************************************
%Results
