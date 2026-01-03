% RR_QPSK
clear all;clc;tic;
noiter=5; Z=4; len=1000;

del=15;  % Decision delay
M=31;  
% D=5; % Dimension of the Subspace
p=M; % Full rank Dimension

 Eb=1; EbNodB=0:1:15; No=1./10.^(EbNodB/10);  Es=log2(Z)*Eb;
 codebook=sqrt(2)*exp(j*2*pi*(0:Z-1)/Z + j*pi/4);
 
% channels 
% W=2.9;
%     for n=1:3,
%         h(n)=.5*(1+cos(2*pi*(n-2)/W)); % RC channel
%     end

h=[.04 -.05 .07 -.21 -.5 .72 .36 0 .21 .03 .07]; % Proakis Channel A
% (Normal)
% h=[.407 .815 .407]; % Proakis Channel B ( Zero Close to unit circle)
% h=[.227 .460 .668 .460 .227]; % Proakis Channel C (Spectral Nulls) 
L=length(h)-1; % ISI length
h_corr=corr(h,h)';
h_ori=[h_corr(L+1:end) zeros(1,M)]; % To eliminate the negative time lag 

for D=6:-1:3,
for loop=1:noiter,
      
    seq=randint(1,len,[0 Z-1]); % Mapping of Bits            
    TRANSseq=sqrt(Es)*exp((j*2*pi*seq/Z)+(j*pi/4));
    
    r_dx_mtx=convolmtx(h,M); % Cross correlation vector (DOES NOT DEPEND ON NOISE)
    r_dx=r_dx_mtx(del+1,:);   % Change here the decision variable

    for i=1:length(No),
               
        N=No(i)*eye(M);  % Noise covariance matrix
        X=var(TRANSseq)*toeplitz(h_ori(1:M));  % Unit variance transmitted correlation
        Y=X+N; % Received Correlation matrix
        MMSE_w=inv(Y)*r_dx';  % Optimum MMSE Full rank solution
        
        y=convol(TRANSseq,h)';      % Channel output 
        leny=length(y);
        noi=randn(1,leny)+j*randn(1,leny);
        noise=sqrt(No(i)/4*log2(Z))*noi;
        rr=y+noise; % RXed signal
        conv_Y=convolmtx(rr,M);
        s=convol(rr,MMSE_w')';  % MMSE filter 
%         RECEIVEDseq=TRANSseq+noise;
        
    for k=1:length(s),
        [val pos]=min(abs(s(k)-codebook));
        decoded(k)=codebook(pos);
    end
    
    number=(decoded(del+1:del+len)==TRANSseq); %Counting the errors. 
    pe(loop,i)=(len-sum(number))/len;
    
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
            % Signal detection for different ranks
            for m=1:D,
            convolmtx_d=conv_Y*V(:,1:m);
            estimated_RR(m,:)=(convolmtx_d*RR_wopt(1:m,m))';
            end
                    
            for k=1:length(estimated_RR(D,:)),
                 [val pos]=min(abs(estimated_RR(D,k)-codebook));
                 decoded(D,k)=codebook(pos);
            end
           
            number=(decoded(D,del+1:del+len)==TRANSseq); %Counting the errors. (Change if del is an array) 
            pe_RR(loop,i)=(len-sum(number))/len;
            clear V T; 
          
        end
    fprintf('Finished %d Iteration for %d bits  in %f seconds  \n',loop,len,toc);
end
ensemble_avg_RR(D,:)=mean(pe_RR);
end

EBNO=10*log10(Eb./No);
ensemble_avg=mean(pe);

% figure,hnd=semilogy(EBNO,erfc(sqrt((Eb./No))).*(1-.25*erfc(sqrt((Eb./No)))),'k',EBNO,(ensemble_avg),'bo-');grid;
% set(hnd(1),'linewidth',2);
% set(hnd(2),'linewidth',2);
% legend('Theoretical','Practical');
% axis([0 11 10e-7 .5]);
% 
% hnd=title('Coherent QPSK Performance Chart');
% set(hnd,'fontsize',15);
% hnd=xlabel(' Eb / No   in   dB ');
% set(hnd,'fontsize',15);
% hnd=ylabel('  Probability  of  Error (Pe)  ');
% set(hnd,'fontsize',15);

%2*Q(sqrt(2*Eb./No)).*(1-.5*Q(sqrt(2*Eb./No)))

% figure,scatter(real(rr),imag(rr));

% figure,hnd=semilogy(EBNO,erfc(sqrt((Eb./No))).*(1-.25*erfc(sqrt((Eb./No)))),'k-',EBNO,ensemble_avg,'bp-',EBNO,ensemble_avg_RR(4,:),'ko-',EBNO,ensemble_avg_RR(3,:),'k*-',EBNO,ensemble_avg_RR(2,:),'k>-');%,EBNO,ensemble_avg_RR(2,:),'b>-');grid;
% grid;
% set(hnd(1),'linewidth',1.5);
% set(hnd(2),'linewidth',1.5);
% set(hnd(3),'linewidth',1.5);
% % set(hnd(4),'linewidth',1.5);
% % set(hnd(5),'linewidth',1.5);
% % set(hnd(6),'linewidth',1.5);
% %set(hnd(7),'linewidth',1.5);
% 
% 
% legend('NO ISI','Full Rank','MSNWF D=4','MSNWF D=3','MSNWF D=2');%,'MSNWF D=2');
% axis([0 12 10e-9 0.5]);
% 
% hnd=xlabel(' Eb/No  in   dB ');
% set(hnd,'fontsize',13,'color',[0 0 1]);
% hnd=ylabel('  Probability  of  Symbol  Error  ');
% set(hnd,'fontsize',13,'color',[0 0 1]);
% 
% Channel 1
% ensemble_avg =[ 0.1819    0.1415    0.1040    0.0715    0.0451    0.0261    0.0131    0.0057    0.0020    0.0006 0.0001    0.0000    0.0000         0         0         0];
%     
% ensemble_avg_RR =[0         0         0         0         0         0         0         0         0         0  0         0         0         0         0         0;...
%     0.1876    0.1496    0.1146    0.0845    0.0599    0.0404    0.0259    0.0161    0.0095    0.0054  0.0027    0.0013    0.0006    0.0002    0.0001    0.0000 ;...
%     0.1820    0.1422    0.1050    0.0726    0.0476    0.0277    0.0150    0.0070    0.0029    0.0010 0.0003    0.0001    0.0000         0         0         0;...
%     0.1820    0.1413    0.1041    0.0719    0.0459    0.0263    0.0134    0.0059    0.0021    0.0006  0.0001    0.0000    0.0000         0         0         0;...
%     0.1812    0.1416    0.1037    0.0717    0.0451    0.0260    0.0132    0.0056    0.0021    0.0005 0.0001    0.0000         0         0         0         0];
    
    
     
    
     
% % for proakis 31 tap equalizer
% ensemble_avg = [0.1892    0.1496    0.1129    0.0810    0.0543    0.0340    0.0193    0.0096  0.0043    0.0016    0.0005    0.0001    0.0000    0.0000         0         0];
% ensemble_avg_RR =[ 0         0         0         0         0         0         0         0    0         0         0         0         0         0         0         0;...
%          0         0         0         0         0         0         0         0   0         0         0         0         0         0         0         0 ;...
%     0.1895    0.1502    0.1140    0.0829    0.0569    0.0372    0.0227    0.0130  0.0071    0.0036    0.0017    0.0008    0.0003    0.0002    0.0001    0.0000; ...
%     0.1897    0.1496    0.1131    0.0813    0.0550    0.0347    0.0203    0.0109  0.0053    0.0023    0.0009    0.0003    0.0001    0.0000    0.0000    0.0000;...
%     0.1901    0.1499    0.1130    0.0810    0.0545    0.0340    0.0192    0.0100 0.0045    0.0017    0.0006    0.0001    0.0000    0.0000    0.0000         0;...
%     0.1902    0.1495    0.1132    0.0812    0.0542    0.0338    0.0193    0.0098  0.0043    0.0016    0.0005    0.0001    0.0000    0.0000    0.0000         0];
%   
%    
% figure,hnd=semilogy(EBNO,erfc(sqrt((Eb./No))).*(1-.25*erfc(sqrt((Eb./No)))),'k-',EBNO,ensemble_avg,'kp-',EBNO,ensemble_avg_RR(5,:),'k*-',EBNO,ensemble_avg_RR(4,:),'ko-',EBNO,ensemble_avg_RR(3,:),'k^-');%,EBNO,ensemble_avg_RR(3,:),'k>-');%,EBNO,ensemble_avg_RR(2,:),'b>-');grid;
% grid;
% set(hnd(1),'linewidth',1.5);
% set(hnd(2),'linewidth',1.5);
% % set(hnd(3),'linewidth',1.5);
% % set(hnd(4),'linewidth',1.5);
% % set(hnd(5),'linewidth',1.5);
% % set(hnd(6),'linewidth',1.5);
% %set(hnd(7),'linewidth',1.5);
% 
% 
% legend('NO ISI','Full Rank','MSNWF D=5','MSNWF D=4','MSNWF D=3');%,'MSNWF D=2');%,'MSNWF D=2');
% axis([0 14 10e-9 0.5]);
% 
% hnd=xlabel(' Eb/No  in   dB ');
% set(hnd,'fontsize',13,'color',[0 0 1]);
% hnd=ylabel('  Probability  of  Error (Pe)  ');
% set(hnd,'fontsize',13,'color',[0 0 1]);
