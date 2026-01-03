% This program for LMS equalization of Known channel (without delayed
% decision ) in a BPSK Transmission ( no unequalized comparison)

clc;clear all;close all;
tic;
noiter=50; % number of averaging iterations
len=1000; % length of the transmitted sequence

del=0;  % Decision delay
M=5;    %MMSE filter length (Why 3?)
o=1000; % Length of the training sequence
% stepsize=.2; % For LMS algorithm

Eb=1; % Energy per bit in V
EbNodB=0:1:15;
noisevar=1./10.^(EbNodB/10);   

h=[1 0.5]; % Assumed ISI Channel model
% h=[.5 1.2 1.5 -1];
% h=[1 1.2 1.5 -1];
L=length(h)-1; % ISI length
% figure,freqz(h);figure,zplane(h)

for loop=1:noiter,   
    
    seq=randint(1,len);
    seq=2*seq-1; % Transmitted sequence
    
        for i=1:length(noisevar),
                     
            % Chennel Estimation  part (Training Period)
            
            train_sym=randint(1,o);  % Training sysmbols
            filtered=convol(train_sym,h)';      % Channel output 
            noisyfiltered=filtered+sqrt(noisevar(i)/2)*randn(1,length(filtered));
%             stepsize=1/(sum(noisyfiltered.^2));
            wopt=LMSequali(noisyfiltered(1:o),train_sym,.1025,M);
            
            y=convol(seq,h)';      % Channel output 
            r=y+sqrt(noisevar(i)/2)*randn(1,length(y)); % RXed signal
%             s=convol(r,wopt(end,:));  % LMS filter 
            rx_u= r > 0; % unequalized sequence
            rx_u1=2*rx_u-1;
            
            % Decision Directed mode.
            [w_opt,s]=LMSDDM(wopt(end,:),r,.01,M);    
            rx_uall=rx_u1(1:len); % unequalized 
            rx_all=s(1:len); 
          
            number=(seq==rx_all); %Counting the errors. 
            number1=(seq==rx_uall);
            
            pe(loop,i)=(len-sum(number))/len;
            pe1(loop,i)=(len-sum(number1))/len;

        end
    fprintf('Finished %d Iteration for %d bits  in %f seconds  \n',loop,len,toc);
end

No=noisevar;
EBNO=10*log10(Eb./No);
ensemble_avg=mean(pe);
ensembleunequali=mean(pe1);

figure,hnd=semilogy(EBNO,.5*erfc(sqrt((Eb./No))),'k',EBNO,(ensemble_avg),'ko-',EBNO,(ensembleunequali),'k*-');grid;
set(hnd(1),'linewidth',2);
set(hnd(2),'linewidth',2);
set(hnd(3),'linewidth',2);
legend('NO ISI','LMS Equalized','Unequalized');
axis([0 11 10e-7 0.5]);

hnd=title('Equlaizer Performance Chart');
set(hnd,'fontsize',15);
hnd=xlabel(' Eb / No   in   dB ');
set(hnd,'fontsize',15);
hnd=ylabel('  Probability  of  Error (Pe)  ');
set(hnd,'fontsize',15);

figure,plot(wopt);
hnd=title('Weight Vector Convergence (Training Phase)');
set(hnd,'fontsize',15);
hnd=xlabel(' Number of Iterations ');
set(hnd,'fontsize',15);
hnd=ylabel('  Weight values  ');
set(hnd,'fontsize',15);

figure,plot(w_opt);
hnd=title('Weight Vector Stability (Equalization Phase)');
set(hnd,'fontsize',15);
hnd=xlabel(' Number of Iterations ');
set(hnd,'fontsize',15);
hnd=ylabel('  Weight values ');
set(hnd,'fontsize',15);