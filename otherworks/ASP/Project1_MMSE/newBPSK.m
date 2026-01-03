% Performance evaluation of BPSK
clear all;
clc;
tic;
noiter=50; seq_length=10000;

Eb=1; % Energy per bit in V
EbNodB=0:1:15;
No=1./10.^(EbNodB/10);   
    
for loop=1:noiter,    
    
    seq=randint(1,seq_length);
    seq=2*seq-1;            % Antipodal mapper
    
    % Generation of PAM signal
    TRANSsignal=sqrt(Eb)*seq;
    for i=1:length(No),
        
        % Noise generation
        noi=randn(1,length(TRANSsignal));
        noise=sqrt(No(i)/2)*noi;
        RECEIVEDsignal=TRANSsignal+noise;
        z=RECEIVEDsignal > 0;     % Decision Threshold = 0 ;
        z1=(2*z)-1;
        number=(seq==z1); %Counting the errors. 
        pe(loop,i)=(seq_length-sum(number))/seq_length;
        
    end
    fprintf('Finished %d Iteration for %d bits  in %f seconds  \n',loop,seq_length,toc);
end

EBNO=10*log10(Eb./No);
ensemble_avg=mean(pe);

figure,hnd=semilogy(EBNO,.5*erfc(sqrt((Eb./No))),'k',EBNO,(ensemble_avg),'ko-');grid;
set(hnd(1),'linewidth',2);
set(hnd(2),'linewidth',2);
legend('Theoretical','Practical');
axis([0 11 10e-9 0.5]);

hnd=title('Coherent BPSK Performance Chart');
set(hnd,'fontsize',15);
hnd=xlabel(' Eb / No   in   dB ');
set(hnd,'fontsize',15);
hnd=ylabel('  Probability  of  Error (Pe)  ');
set(hnd,'fontsize',15);
