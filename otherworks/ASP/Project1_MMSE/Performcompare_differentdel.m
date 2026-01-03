% This program for MMSE equalization of Known channel (with delayed
% decision ) in a BPSK Transmission ( with unequlai comparison for different del);

clc;clear all;close all;
tic;
noiter=50; % number of averaing iterations
len=10000; % length of the transmiited sequence

del=[0 1 2 3 ];  % Decision delay
M=3;    %MMSE filter length (Why 3?)
Eb=1; % Energy per bit in V
EbNodB=0:1:15;

% h=[1 0.5]; % Assumed ISI Channel model
h=[-.5 1.2 1.5 -1]; 
% h=[1 1.2 1.5 -1];
L=length(h)-1; % ISI length
% figure,freqz(h);figure,zplane(h)

h_corr=corr(h,h)';
h_ori=[h_corr(L+1:end) zeros(1,M)]; % To eliminate the negative time lag 
r_dx_mtx=convolmtx(h,M); % Cross correlation vector (DOES NOT DEPEND ON NOISE VARIANCE)

for m=1:length(del),
    r_dx=r_dx_mtx(m,:);
for loop=1:noiter,   
    
    seq=randint(1,len);
    seq=2*seq-1; % Transmitted sequence

    noisevar=1./10.^(EbNodB/10);   
      
      for i=1:length(noisevar),
            
            N=noisevar(i)*eye(M);  % Noise covariance matrix
            X=toeplitz(h_ori(1:M));  % Unit variance transmitted correlation
            Y=X+N; % Received Correlation matrix

            MMSE_w=inv(Y)*r_dx'; % Optimum MMSE solution

            y=convol(seq,h)';      % Channel output 
            r=y+sqrt(noisevar(i)/2)*randn(1,length(y)); % RXed signal
            s=convol(r,MMSE_w')';  % MMSE filter 
            rx = s > 0;  % Decision Box (equalized sequence)
            rx1=(2*rx-1);
                  
            rx_all=rx1(m:len+m-1); % L= ISI length
            number=(seq==rx_all); %Counting the errors.
                          
            pe(loop,i)=(len-sum(number))/len;
          
        end
        
    fprintf('Finished %d Iteration for %d bits  in %f seconds  \n',loop,len,toc);
    
end

fprintf('Finished decoding for %d delayed decision sequence  \n',m-1);
ensemble_avg(m,:)=mean(pe);
end

No=noisevar;
EBNO=10*log10(Eb./No);

% figure,hnd=semilogy(EBNO,.5*erfc(sqrt((Eb./No))),'k',EBNO,(ensemble_avg),'ko-',EBNO,ensembleunequali,'ro-');grid;
figure,hnd=semilogy(EBNO,ensemble_avg(1,:),'ro-',EBNO,ensemble_avg(2,:),'go-',EBNO,ensemble_avg(3,:),'bo-',EBNO,ensemble_avg(4,:),'r*-.');grid;
set(hnd(1),'linewidth',2);
set(hnd(2),'linewidth',2);
set(hnd(3),'linewidth',2);
set(hnd(4),'linewidth',2);
legend('\Delta = 0','\Delta = 1','\Delta = 2','\Delta = 3');
axis([0 17 10e-6 0.5]);

hnd=title('MMSE equalizer with Delayed decision');
set(hnd,'fontsize',15);
hnd=xlabel(' Eb / No   in   dB ');
set(hnd,'fontsize',15);
hnd=ylabel('  Probability  of  Error (Pe)  ');
set(hnd,'fontsize',15);

% figure,hnd=scatter(seq,zeros(1,length(seq)),60,[1 0 0]);
% set(hnd,'linewidth',1.8);
% axis([-2 2 -.1 .1]);grid;hold on;
% scatter(y,zeros(1,length(y)),20,[0 0 1]);
% axis([-10 10 -.1 .1]); 