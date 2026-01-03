% This program for LS equalization of Known channel (with0ut delayed
% decision ) in a BPSK Transmission ( with unequlai comparison)

clc;clear all;close all;
tic;
noiter=10; % number of averaing iterations
len=100000; % length of the transmiited sequence
LT=64; % Length of the training sequence for channel estimation

del=0;  % Decision delay
M=3;    %MMSE filter length (Why 3?)

Eb=1; % Energy per bit in V
EbNodB=0:1:15;
noisevar=1./10.^(EbNodB/10);

h=[1 0.4]; % Assumed ISI Channel model
% h=[.5 1.2 1.5 -1];
% h=[1 1.2 1.5 -1];
L=length(h)-1; % ISI length
% figure,freqz(h);figure,zplane(h)

h_corr=corr(h,h)';
h_ori=[h_corr(L+1:end) zeros(1,M)]; % To eliminate the negative time lag 

temp(1)=1;
for i=2:L+1,
    temp(i)=xor(1,temp(i-1));
end
temp=2*temp-1;

for loop=1:noiter,   
    
    seq=randint(1,len);
    seq=2*seq-1; % Transmitted sequence

    r_dx_mtx=convolmtx(h,M); % Cross correlation vector (DOES NOT DEPEND ON NOISE)
    r_dx=r_dx_mtx(del+1,:);   % Change here the decision variable

        for i=1:length(noisevar),
            
            N=noisevar(i)*eye(M);  % Noise covariance matrix
            X=var(seq)*toeplitz(h_ori(1:M));  % Unit variance transmitted correlation
            Y=X+N; % Received Correlation matrix

            MMSE_w=inv(Y)*r_dx';  % Optimum MMSE solution
            
            % LS Channel estimation part
            train_seq=randn(1,LT);
            H=LSDataform(train_seq,L+1,3);
            channel_train_out=H*h'+sqrt(noisevar(i)/2)*randn(1,length(H))';
            LS_esti=inv(H'*H)*H'*channel_train_out;   
            LS_QMF=temp.*LS_esti';

            y=convol(seq,h)';      % Channel output 
            r=y+sqrt(noisevar(i)/2)*randn(1,length(y)); % RXed signal
            
            s=convol(r,MMSE_w')';  % MMSE filter 
            LS_s=convol(r,LS_QMF);
            
            rx = s > 0;  % Decision Box
            LS_rx=LS_s > 0;
            
            rx1=(2*rx-1);
            LS_rx1=(2*LS_rx)-1;
            
            rx_all=rx1(1:len); % L= ISI length
            LS_all=LS_rx1(1:len);
            
            number=(seq==rx_all); %Counting the errors. (Change if del is an array) 
            LS_number=(seq==LS_all');
            
            pe(loop,i)=(len-sum(number))/len;
            LS_pe(loop,i)=(len-sum(LS_number))/len;
%             SNR(loop,i)=10*log10((conj(MMSE_w')*X*MMSE_w)/(conj(MMSE_w')*N*MMSE_w)); % For equalized
%             SNRU(loop,i)=10*log10((var(seq)*h_ori(1)/(noisevar(i)))); % For Unequalized

        end
    fprintf('Finished %d Iteration for %d bits  in %f seconds  \n',loop,len,toc);
end

No=noisevar;
EBNO=10*log10(Eb./No);
ensemble_avg=mean(pe);
ensemble_avg_LS=mean(LS_pe);
% avgSNR=mean(SNR);
% avgSNRU=mean(SNRU);

% ImprovedSNR=10*log10(

figure,hnd=semilogy(EBNO,.5*erfc(sqrt((Eb./No))),'b',EBNO,(ensemble_avg),'ro-',EBNO,(ensemble_avg_LS),'ko-');grid;
set(hnd(1),'linewidth',2);
set(hnd(2),'linewidth',2);
set(hnd(3),'linewidth',2);
legend('NO ISI','MMSE Equalized','LS Equalized');
axis([0 11 10e-9 0.5]);

hnd=title('Equlaizer Performance Chart');
set(hnd,'fontsize',15);
hnd=xlabel(' Eb / No   in   dB ');
set(hnd,'fontsize',15);
hnd=ylabel('  Probability  of  Error (Pe)  ');
set(hnd,'fontsize',15);

% SNR improvement chart

% figure,hnd=semilogy(EBNO,avgSNRU,'b',EBNO,avgSNR,'ro-');grid;
% set(hnd(1),'linewidth',2);
% set(hnd(2),'linewidth',2);
% legend('Un Equalized','MMSE Equalized');
