% This program for LS equalization of Known channel (with0ut delayed
% decision ) in a BPSK Transmission ( with unequlai comparison)

clc;clear all;close all;
tic;
noiter=50; % number of averaing iterations
len=10000; % length of the transmiited sequence
LT=64; % Length of the training sequence for channel estimation

del=0;  % Decision delay


Eb=1; % Energy per bit in V
EbNodB=0:1:15;
noisevar=1./10.^(EbNodB/10);

h=[1 0.4]; % Assumed ISI Channel model
% h=[.5 1.2 1.5 -1];
% h=[1 1.2 1.5 -1];
L=length(h)-1; % ISI length
% figure,freqz(h);figure,zplane(h)

% temp(1)=1;
% for i=2:L+1,
%     temp(i)=xor(1,temp(i-1));
% end
% temp=2*temp-1;

for loop=1:noiter,   
    
    seq=randint(1,len);
    seq=2*seq-1; % Transmitted sequence

       for i=1:length(noisevar),
            
            % LS Channel estimation part
            train_seq=randn(1,LT);
            H=LSDataform(train_seq,L+1,3);
            channel_train_out=H*h'+sqrt(noisevar(i))*randn(1,length(H))';
            LS_esti=inv(H'*H)*H'*channel_train_out;   
%             LS_QMF=temp.*LS_esti';

            y=convol(seq,h)';      % Channel output 
            r=y+sqrt(noisevar(i)/2)*randn(1,length(y)); % RXed signal
            
%             LS_s=convol(r,LS_QMF);
            LS_s=filter(1,LS_esti',r);
%             LS_s=filter(1,h,r);
            
            rx_u= r > 0; % unequalized sequence
            LS_rx=LS_s > 0;
            
            rx1=(2*rx_u-1);
            LS_rx1=(2*LS_rx)-1;
            
            rx_all=rx1(1:len); % L= ISI length
            LS_all=LS_rx1(1:len);
            
            number=(seq==rx_all); %Counting the errors. (Change if del is an array) 
            LS_number=(seq==LS_all);
            
            pe(loop,i)=(len-sum(number))/len;
            LS_pe(loop,i)=(len-sum(LS_number))/len;

            
end
    fprintf('Finished %d Iteration for %d bits  in %f seconds  \n',loop,len,toc);
end

No=noisevar;
EBNO=10*log10(Eb./No);
ensemble_avg=mean(pe);
ensemble_avg_LS=mean(LS_pe);

figure,hnd=semilogy(EBNO,.5*erfc(sqrt((Eb./No))),'b',EBNO,(ensemble_avg),'ro-',EBNO,(ensemble_avg_LS),'ko-');grid;
set(hnd(1),'linewidth',2);
set(hnd(2),'linewidth',2);
set(hnd(3),'linewidth',2);
legend('NO ISI','Unequalized','LS Equalized');
axis([0 13 10e-6 0.5]);

hnd=title('Equalizer Performance Chart');
set(hnd,'fontsize',15);
hnd=xlabel(' Eb / No   in   dB ');
set(hnd,'fontsize',15);
hnd=ylabel('  Probability  of  Error (Pe)  ');
set(hnd,'fontsize',15);