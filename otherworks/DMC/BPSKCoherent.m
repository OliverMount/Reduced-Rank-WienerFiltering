% Performance evaluation of BPSK
clear all;clc;tic;
noiter=500; % no of iterations
seq_length=1000; % no of bits to be transmitted

for loop=1:noiter,    
    
    seq=randint(1,seq_length); % sequence of 0's and 1's
    seq=2*seq-1;            % Antipodal mapper (-1 and 1)
    
    Eb=1; % Energy per bit in V
    SNR=0:1:12; % range of SNR we want to perfrom 
    No=Eb./10.^(SNR/10);   % Noise variance values SNR=10*log10(Eb/No)
    
    % Generation of PAM signal
    
    for i=1:length(No),
        TRANSsignal(i,:)=seq;
        % Noise generation
        noi=randn(1,length(TRANSsignal)); %zero mean unit variance Gaussina noise
        noise=sqrt(No(i)/2)*noi; % Converting the variance to (no/2) from 1
        RECEIVEDsignal(i,:)=TRANSsignal(i,:)+noise; % Signal + Noise
        z(i,:)=RECEIVEDsignal(i,:) > 0;     % Decision Threshold = 0 ;
        z1(i,:)=(2*z(i,:))-1;
        number=(seq==z1(i,:)); %Counting the errors. 
        pe(loop,i)=(seq_length-sum(number))/seq_length;
        
    end
    fprintf('Finished %d Iteration for %d bits  in %f seconds  \n',loop,seq_length,toc);
end

ensemble_avg=mean(pe); % average Probability of Error

figure,hnd=semilogy(SNR,0.5*erfc(sqrt((Eb./No))),'k',SNR,(ensemble_avg),'ko-');grid;
set(hnd(1),'linewidth',2);
set(hnd(2),'linewidth',2);
legend('Theoretical','Practical');
axis([0 12 10e-8 0.5]);

hnd=title('Coherent BPSK Performance Chart');
set(hnd,'fontsize',15);
hnd=xlabel(' SNR   in   dB ');
set(hnd,'fontsize',15);
hnd=ylabel('  Probability  of  Bit Error (Pe)  ');
set(hnd,'fontsize',15);
