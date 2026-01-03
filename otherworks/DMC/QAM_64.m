% Performance evaluation of M-QAM
clear all;clc;tic;
noiter=100; M=4; len=1000;

SNR=0:2:20;  Es=1; No=Es./10.^(SNR/10); 
codebook=QAM_codebook(M);

% codebook=QAM_codebook(M);
% avgE=mean(abs(codebook).^2); % Average Energy of the constellation (Usually one)

for loop=1:noiter,
      
    seq=randint(1,len,[1 M]); % Mapping of Bits            
    TRANSseq=codebook(seq);
    % Generation of M-QAM signal
    for i=1:length(No),
        % Noise generation
        noi=sqrt(1/2)*(randn(1,len)+j*randn(1,len));
        noise=sqrt(No(i)/log2(M))*noi;
        RECEIVEDseq=TRANSseq+noise;
        
    for k=1:len,
        [val pos]=min(abs(RECEIVEDseq(k)-codebook));
        decoded(k)=codebook(pos);
    end
    
    number=(decoded==TRANSseq); %Counting the errors. 
    pe(loop,i)=(len-sum(number))/len;
end

    fprintf('Finished %d Iteration for %d symbols  in %f seconds  \n',loop,len,toc);
end

ensemble_avg=mean(pe);

figure,hnd=semilogy(SNR,ensemble_avg,'k-');grid;
set(hnd(1),'linewidth',2);
% set(hnd(2),'linewidth',2);
% legend('Theoretical','Practical');
axis([min(SNR) 20 10e-6 1]);

hnd=title('Coherent 64-QAM Performance Chart');
set(hnd,'fontsize',15);
hnd=xlabel('SNR   in   dB ');
set(hnd,'fontsize',15);
hnd=ylabel('  Probability  of  Symbol  Error  ');
set(hnd,'fontsize',15);

% hnd=scatter(real(codebook),imag(codebook)); set(hnd,'linewidth',2);