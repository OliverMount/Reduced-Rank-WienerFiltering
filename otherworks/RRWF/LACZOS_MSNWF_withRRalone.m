            
% This program for MMSE RR Equlaization 
clc;clear all;close all;
tic;
noiter=75; % number of averaing iterations
len=100000; % length of the transmiited sequence

del=15;  % Decision delay
M=31;  
% D=5; % Dimension of the Subspace
p=M; % Full rank Dimension

Eb=1; % Energy per bit in V
EbNodB=0:12;
noisevar=1./10.^(EbNodB/10);

% h=[1 0.5]; % Assumed ISI Channel model % Ali sayeed channel
h=[.04 -.05 .07 -.21 -.5 .72 .36 0 .21 .03 .07]; % Proakis Channel A
% (Normal)
% h=[.407 .815 .407]; % Proakis Channel B ( Zero Close to unit circle)
% h=[.227 .460 .668 .460 .227]; % Proakis Channel C (Spectral Nulls) 
L=length(h)-1; % ISI length
% figure,freqz(h);figure,zplane(h)
% figure,stem(h);

h_corr=corr(h,h)';
h_ori=[h_corr(L+1:end) zeros(1,M)]; % To eliminate the negative time lag 

for D=7:-1:2,
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

%              % unequalized stream
%             rx1_u=rr >0;
%             rx1_u=(2*rx1_u-1);
%             rx_all_u=rx1_u(1:L+len); % L= ISI length
%             number_u=(seq==rx_all_u(1:len)); %Counting the errors. (Change if del is an array) 
%             pe_u(loop,i)=(len-sum(number_u))/len;
%                        

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
            
            rx_RR = estimated_RR(D,:) > 0;  % Decision Box 
            rx1=(2*rx_RR-1);
            rx_all=rx1(del+1:del+len); % L= ISI length
           
            number=(seq==rx_all); %Counting the errors. (Change if del is an array) 
            pe_RR(loop,i)=(len-sum(number))/len;
            clear V T;
            
             end
    fprintf('Finished %d Iteration for %d bits  in %f seconds  \n',loop,len,toc);
end

ensemble_avg_RR(D,:)=mean(pe_RR);
end

No=noisevar;
EBNO=10*log10(Eb./No);
ensemble_avg=mean(pe);
% ensemble_avg_u=mean(pe_u);

figure,hnd=semilogy(EBNO,.5*erfc(sqrt((Eb./No))),'b',EBNO,ensemble_avg,'r',EBNO,ensemble_avg_RR(6,:),'g-',EBNO,ensemble_avg_RR(5,:),'b',EBNO,ensemble_avg_RR(4,:),'ro',EBNO,ensemble_avg_RR(3,:),'k>-',EBNO,ensemble_avg_RR(2,:),'b>-');grid;
set(hnd(1),'linewidth',1.5);
set(hnd(2),'linewidth',2);
set(hnd(3),'linewidth',1.5);
set(hnd(4),'linewidth',1.5);
set(hnd(5),'linewidth',1.5);
set(hnd(6),'linewidth',1.5);
set(hnd(7),'linewidth',1.5);


legend('NO ISI','Full Rank','MSNWF D=6','MSNWF D=5','MSNWF D=4','MSNWF D=3','MSNWF D=2');
axis([0 14 10e-9 0.5]);

hnd=xlabel(' Eb/No  in   dB ');
set(hnd,'fontsize',13,'color',[0 0 1]);
hnd=ylabel('  Probability  of  Error (Pe)  ');
set(hnd,'fontsize',13,'color',[0 0 1]);
% ensemble_avg =
% 
%   Columns 1 through 12 
% 
%     0.0986    0.0762    0.0557    0.0386    0.0243    0.0141    0.0071    0.0032    0.0011    0.0003    0.0001    0.0000
% 
%   Columns 13 through 16 
% 
%     0.0000         0         0         0
%     
% ensemble_avg_RR =
% 
%   Columns 1 through 12 
% 
%          0         0         0         0         0         0         0         0         0         0         0         0
%     0.1009    0.0796    0.0604    0.0444    0.0309    0.0209    0.0134    0.0083    0.0048    0.0026    0.0014    0.0007
%     0.0983    0.0763    0.0564    0.0391    0.0254    0.0151    0.0079    0.0037    0.0015    0.0005    0.0001    0.0000
%     0.0980    0.0762    0.0559    0.0387    0.0245    0.0143    0.0072    0.0031    0.0011    0.0003    0.0001    0.0000
% 
%   Columns 13 through 16 
% 
%          0         0         0         0
%     0.0003    0.0001    0.0000    0.0000
%     0.0000    0.0000         0         0
%     0.0000         0         0         0


% ensemble_avg =
% 
%   Columns 1 through 12 
% 
%     0.0984    0.0760    0.0559    0.0384    0.0245    0.0141    0.0072    0.0031    0.0011    0.0003    0.0001    0.0000
% 
%   Columns 13 through 16 
% 
%          0         0         0         0
%          
%   ensemble_avg_RR =
% 
%   Columns 1 through 12 
% 
%          0         0         0         0         0         0         0         0         0         0         0         0
%     0.1006    0.0795    0.0605    0.0441    0.0312    0.0209    0.0134    0.0083    0.0049    0.0027    0.0014    0.0007
%     0.0987    0.0766    0.0565    0.0389    0.0253    0.0150    0.0079    0.0038    0.0015    0.0005    0.0002    0.0000
%     0.0984    0.0762    0.0560    0.0389    0.0245    0.0145    0.0072    0.0033    0.0012    0.0003    0.0001    0.0000
% 
%   Columns 13 through 16 
% 
%          0         0         0         0
%     0.0003    0.0001    0.0000    0.0000
%     0.0000    0.0000         0         0
%          0         0         0         0



% ensemble_avg_RR =
% 
%   Columns 1 through 12 
% 
%          0         0         0         0         0         0         0         0         0         0         0         0
%     0.1006    0.0797    0.0606    0.0444    0.0311    0.0209    0.0133    0.0082    0.0048    0.0027    0.0014    0.0007
%     0.0986    0.0763    0.0563    0.0392    0.0253    0.0148    0.0080    0.0037    0.0015    0.0005    0.0002    0.0000
%     0.0986    0.0762    0.0562    0.0387    0.0246    0.0141    0.0072    0.0032    0.0012    0.0003    0.0001    0.0000
% 
%   Columns 13 through 16 
% 
%          0         0         0         0
%     0.0003    0.0001    0.0000    0.0000
%     0.0000    0.0000         0         0
%     0.0000         0         0         0
%     
%     ensemble_avg =
% 
%   Columns 1 through 12 
% 
%     0.0983    0.0763    0.0560    0.0386    0.0245    0.0141    0.0072    0.0031    0.0011    0.0003    0.0001    0.0000
% 
%   Columns 13 through 16 
% 
%     0.0000    0.0000         0         0

% 
% For Proakis Channel
% 
% ensemble_avg_RR1 =[  0         0         0         0         0         0         0         0         0   0         0         0         0 ;
%     0.1025    0.0814    0.0625    0.0464    0.0330    0.0229    0.0154    0.0100    0.0064   0.0040    0.0025    0.0015    0.0010;
%     0.1015    0.0797    0.0601    0.0433    0.0299    0.0193    0.0118    0.0068    0.0036  0.0019    0.0009    0.0004    0.0002;
%     0.1016    0.0795    0.0600    0.0430    0.0292    0.0186    0.0109    0.0059    0.0029 0.0012    0.0005    0.0002    0.0001;
%     0.1016    0.0794    0.0599    0.0429    0.0291    0.0183    0.0106    0.0055    0.0025    0.0010    0.0003    0.0001    0.0000;
%     0.1014    0.0795    0.0599    0.0429    0.0290    0.0182    0.0106    0.0054    0.0025 0.0010    0.0003    0.0001    0.0000;
%     0.1013    0.0797    0.0598    0.0429    0.0290    0.0182    0.0105    0.0054    0.0025    0.0009    0.0003    0.0001    0.0000];
% ensemble_avg1=[ 0.1015    0.0798    0.0600    0.0431    0.0290    0.0183    0.0106    0.0055    0.0025    0.0009    0.0003    0.0001    0.0000 ];
