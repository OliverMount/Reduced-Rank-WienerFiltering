% OFDM simulation (QPSK) modulation with MMSE symbol estimation,symbol by symbol detection (ML)
% obtained QPSK performance

clear all;clc;close all;tic;
N=32; % No of sub-carriers
O=100; % no of OFDM symbols
len=N*O; % total number of QPSK symbols 
noiter=5;
Z=4; % Signifies QPSK
Eb=1; EbNodB=0:1:15; No=1./10.^(EbNodB/10);  Es=log2(Z)*Eb;
codebook=sqrt(Es)*exp(j*2*pi*(0:Z-1)/Z + j*pi/4); % QPSK codebook
F=(1/sqrt(N))*dftmtx(N);

 H=CircMtx([1],N); % Channel matrix for AWGN
 [EE,E]=eig(H);
 
 for loop=1:noiter,
  
 seq=randint(1,len,[0 Z-1]);  % generation of symbols      
 X_s=sqrt(Es)*exp((j*2*pi*seq/Z)+(j*pi/4)); % Mapping of symbols
 X=SP(X_s,N); % serial to parallel conversion
 
 for i=1:length(No),
    
     for q=1:O, % upto  no of OFDM symbols
        noi=randn(1,N)+j*randn(1,N);
        noise=sqrt(No(i)/2)*noi.';
        x=F'*X(:,q); % taking IFFT before transmisssion
        y(:,q)=H*x+noise; % Passing thro Channel  + noise (received symbols)
        Y_p(:,q)=F*y(:,q); % FFT after reception
        Y_decoded(:,q)=E'*inv(E*E'+(No(i)/var(X_s))*eye(N))*Y_p(:,q); % MMSE estimator
    end
    Y=PS(Y_decoded,N); % Parallel to serial conversion
    %ML decoding
    for k=1:length(Y),
        [val pos]=min(abs(Y(k)-codebook));
        decoded(k)=codebook(pos);
    end    

     number=(decoded==X_s); %Counting the errors. 
     pe(loop,i)=(len-sum(number))/len;
end
fprintf('Finished %d Iteration for %d symbols in %f seconds  \n',loop,len,toc);
end

EBNO=10*log10(Eb./No);
ensemble_avg=mean(pe);

figure,hnd=semilogy(EBNO,erfc(sqrt((Eb./No))).*(1-.25*erfc(sqrt((Eb./No)))),'k-',EBNO,ensemble_avg);grid;
set(hnd(1),'linewidth',1.5);
set(hnd(2),'linewidth',1.5);
% set(hnd(3),'linewidth',1.5);
% set(hnd(4),'linewidth',1.5);
% set(hnd(5),'linewidth',1.5);
% set(hnd(6),'linewidth',1.5);
%set(hnd(7),'linewidth',1.5);
legend('Theoretical','Practical');
axis([0 12 10e-8 0.5]);
hnd=xlabel(' Eb/No  in   dB ');
set(hnd,'fontsize',13,'color',[0 0 1]);
hnd=ylabel('  Probability  of  Error (Pe)  ');
set(hnd,'fontsize',13,'color',[0 0 1]);

