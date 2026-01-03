% This program for MMSE equalization of Known channel (with delayed
% decision ) in a BPSK Transmission ( no unequlai comparison)

% RR_QPSK
clear all;clc;tic;
noiter=2; Z=4; len=1000;
startr=5;endr=2;
del=35;  % Decision delay (From del=5 started getting results) (Note that decision dely could be longer than equlaizer length)
M=31;  % Number of Taps of the equalizer
% D=5; % Dimension of the Subspace
p=M; % Full rank Dimension

 Eb=1; EbNodB=0:14; No=1./10.^(EbNodB/10);  Es=log2(Z)*Eb;
 codebook=sqrt(2)*exp(j*2*pi*(0:Z-1)/Z + j*pi/4);
 
% channels 
% W=2.9;
%     for n=1:3,
%         h(n)=.5*(1+cos(2*pi*(n-2)/W)); % RC channel
%     end

% 11-tap Data quality telephone channel
h=[0+j*0 0.0485+j*0.0194 0.0573+j*0.0253 0.0786+j*0.0282 0.0874+j*0.0447 0.9222+j*0.3031 0.1427+j*0.0349 0.0835+j*0.0157 0.0621+j*0.0078 0.0359+j*0.0049 0.0214+j*0.0019];
% h=[.04 -.05 .07 -.21 -.5 .72 .36 0 .21 .03 .07]; % Proakis Channel A
% (Normal)
% h=[.407 .815 .407]; % Proakis Channel B ( Zero Close to unit circle)
% h=[.227 .460 .668 .460 .227]; % Proakis Channel C (Spectral Nulls) 
L=length(h)-1; % ISI length
h_corr=corr(h,h).';
h_ori=[h_corr(L+1:end) zeros(1,M)]; % To eliminate the negative time lag 

for D=startr:-1:endr,
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
        noise=sqrt(No(i)/4*log2(Z))*noi;
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
            pe_RR(loop,i)=(len-sum(number))/len;
            clear V T; 
          
        end
    fprintf('Finished %d Iteration for %d bits  in %f seconds  \n',loop,len,toc);
end
ensemble_avg_RR(D,:)=mean(pe_RR);
end

EBNO=10*log10(Eb./No);
ensemble_avg=mean(pe);

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
% hnd=ylabel('  Probability  of  Error (Pe)  ');
% set(hnd,'fontsize',15);

%2*Q(sqrt(2*Eb./No)).*(1-.5*Q(sqrt(2*Eb./No)))

% figure,scatter(real(rr),imag(rr));

figure,hnd=semilogy(EBNO,erfc(sqrt((Eb./No))).*(1-.25*erfc(sqrt((Eb./No)))),'k-',EBNO,ensemble_avg,'bp-',EBNO,ensemble_avg_RR(startr-1,:),'ko-',EBNO,ensemble_avg_RR(startr-2,:),'k*-',EBNO,ensemble_avg_RR(startr-3,:),'k>-');grid;
% figure,hnd=semilogy(EBNO,erfc(sqrt((Eb./No))).*(1-.25*erfc(sqrt((Eb./No)))),'k-',EBNO,ensemble_avg,'bp-',EBNO,ensemble_avg_RR(startr,:),'ko-',EBNO,ensemble_avg_RR(startr-1,:),'k*-',EBNO,ensemble_avg_RR(startr-2,:),'k>-',EBNO,ensemble_avg_RR(startr-3,:),'b>-');grid;

set(hnd(1),'linewidth',1.5);
set(hnd(2),'linewidth',1.5);
set(hnd(3),'linewidth',1.5);
% set(hnd(4),'linewidth',1.5);
% set(hnd(5),'linewidth',1.5);
% set(hnd(6),'linewidth',1.5);
%set(hnd(7),'linewidth',1.5);
% 
% 
% legend('NO ISI','Full Rank','MSNWF D=4','MSNWF D=3','MSNWF D=2','MSNWF D=2'); % use for 4 ranks
legend('NO ISI','Full Rank','MSNWF D=5','MSNWF D=4','MSNWF D=3'); % use for 4 ranks
axis([0 13 10e-7 1]);
% 
% hnd=xlabel(' Eb/No  in   dB ');
% set(hnd,'fontsize',13,'color',[0 0 1]);
% hnd=ylabel('  Probability  of  Error (Pe)  ');
% set(hnd,'fontsize',13,'color',[0 0 1]);

% Results for D=5:-1:2 (N= 21 tap equalizer delay 25) ( 4 3 2 graph is good
% although we have result for rank 5 also)
% ensemble_avg_RR =
% 
%   Columns 1 through 13 
% 
%          0         0         0         0         0         0         0         0         0         0         0         0         0
%     0.2102    0.1634    0.1222    0.0865    0.0582    0.0371    0.0222    0.0124    0.0068    0.0035    0.0017    0.0008    0.0004
%     0.2071    0.1594    0.1163    0.0799    0.0520    0.0314    0.0175    0.0090    0.0043    0.0018    0.0008    0.0002    0.0001
%     0.2070    0.1595    0.1157    0.0795    0.0512    0.0308    0.0170    0.0085    0.0041    0.0016    0.0006    0.0002    0.0000
%     0.2072    0.1581    0.1159    0.0800    0.0511    0.0303    0.0170    0.0085    0.0038    0.0015    0.0006    0.0002    0.0001
% 
%   Columns 14 through 15 
% 
%          0         0
%     0.0002    0.0001
%     0.0000         0
%     0.0000    0.0000
%     0.0000    0.0000