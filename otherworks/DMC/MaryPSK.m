% Performance evaluation of QPSK alone here
clear all;clc;tic;
Z=4; % Signifies the contellation size
noiter=100;  
len=1000;

SNR=0:1:16;  Es=1; No=Es./10.^(SNR/10);
codebook=(1/sqrt(2))*[-1-j -1+j 1-j 1+j];
% avgE=mean(abs(codebook).^2); % Average Energy of the constellation 

for loop=1:noiter,
      
    seq=randint(1,len,[1 Z]); % Mapping of Bits            
    TRANSseq=codebook(seq);

    % Generation of MPSK signal
    for i=1:length(No),
        % Noise generation
        noi=sqrt(1/2)*(randn(1,len)+j*randn(1,len));
        noise=sqrt(No(i))*noi;
        RECEIVEDseq=TRANSseq+noise;
        
    for k=1:len,
        [val pos]=min(abs(RECEIVEDseq(k)-codebook));
        decoded(k)=codebook(pos);
    end
    
    number=(decoded==TRANSseq); %Counting the errors. 
    pe(loop,i)=(len-sum(number))/len;
end

    fprintf('Finished %d Iteration for %d bits  in %f seconds  \n',loop,len,toc);
end

ensemble_avg=mean(pe/log2(Z));

figure,hnd=semilogy(SNR,ensemble_avg,'k-');grid;
set(hnd(1),'linewidth',2);
% set(hnd(2),'linewidth',2);
% legend('Theoretical','Practical');
axis([min(SNR) max(SNR) 10e-7 1]);

hnd=title('Coherent M-PSK Performance Chart');
set(hnd,'fontsize',15);
hnd=xlabel('SNR   in   dB ');
set(hnd,'fontsize',15);
hnd=ylabel('  Probability  of  Symbol  Error  ');
set(hnd,'fontsize',15);