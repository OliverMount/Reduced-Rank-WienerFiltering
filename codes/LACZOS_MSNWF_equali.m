% This program for MMSE equalization of Known channel (with delayed
% decision ) in a BPSK Transmission ( no unequlai comparison)

clc;clear all;close all;
tic;
noiter=100; % number of averaing iterations
len=10000; % length of the transmiited sequence


del=7;  % Decision delay
M=11;    %MMSE filter length (Why 3?)

Eb=1; % Energy per bit in V
EbNodB=0:1:15;
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

            MMSE_w=inv(Y)*r_dx';  % Optimum MMSE solution
            
            y=convol(seq,h)';      % Channel output 
            r=y+sqrt(noisevar(i)/2)*randn(1,length(y)); % RXed signal
            s=convol(r,MMSE_w')';  % MMSE filter 
            rx = s > 0;  % Decision Box 
            rx1=(2*rx-1);
            rx_all=rx1(del+1:del+len); % L= ISI length
           
            number=(seq==rx_all); %Counting the errors. (Change if del is an array) 
            pe(loop,i)=(len-sum(number))/len;
%             SNR(loop,i)=10*log10((conj(MMSE_w')*X*MMSE_w)/(conj(MMSE_w')*N*MMSE_w)); % For equalized
%             SNRU(loop,i)=10*log10((var(seq)*h_ori(1)/(noisevar(i)))); % For Unequalized

            MSE(loop,i)=h_ori(1)-MMSE_w'*X*MMSE_w;

        end
    fprintf('Finished %d Iteration for %d bits  in %f seconds  \n',loop,len,toc);
end

No=noisevar;
EBNO=10*log10(Eb./No);
ensemble_avg=mean(pe);
% avgSNR=mean(SNR);
% avgSNRU=mean(SNRU);

% ImprovedSNR=10*log10(

figure,hnd=semilogy(EBNO,.5*erfc(sqrt((Eb./No))),'b',EBNO,(ensemble_avg),'ro-');grid;
set(hnd(1),'linewidth',2);
set(hnd(2),'linewidth',2);
legend('NO ISI','MMSE Equalized');
axis([0 11 10e-9 0.5]);

hnd=title('Equalizer Performance Chart');
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


Avg_MSE=mean(MSE)