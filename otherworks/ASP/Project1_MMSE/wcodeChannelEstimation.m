% This program for LS equalization of Known channel with Wavelet Codes (with0ut delayed
% decision ) in a BPSK Transmission ( with unequlai comparison)

clc;clear all;close all;
tic;
noiter=100; % number of averaing iterations
len=10000; % length of the transmiited sequence
% LT=64; % Length of the training sequence for channel estimation

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


train_seq1= [ -0.1093   -0.3027    0.0464    0.2502   -0.0193   -0.1619    0.0195    0.1141   -0.0019   -0.0616   0.0118    0.0495   -0.0121   -0.0470    0.0050    0.0337    0.0056   -0.0073    0.0112    0.0214 ...
  -0.0185   -0.0450    0.0184    0.0578   -0.0232   -0.0719    0.0205    0.0751   -0.0205   -0.0758   0.0223    0.0794   -0.0403   -0.1126   -0.0125    0.0414    0.0406    0.0544   -0.0497   -0.1222 ...
   0.0712    0.1934   -0.0737   -0.2366    0.0829    0.2731   -0.0971   -0.3172    0.1440    0.4226   -0.0452   -0.3130   -0.0028    0.1525    0.0117   -0.0497   -0.0461   -0.0599    0.0414    0.1119 ...
  -0.0504   -0.1503    0.0698    0.2041 ];

train_seq=[train_seq1 2*train_seq1(1:32)];


for loop=1:noiter,   
    
    seq=randint(1,len);
    seq=2*seq-1; % Transmitted sequence

       for i=1:length(noisevar),
            
            % LS Channel estimation part
%             train_seq=randn(1,LT);
            H=LSDataform(train_seq,L+1,3);
            channel_train_out=H*h'+sqrt(noisevar(i))*randn(1,length(H))';
            LS_esti=inv(H'*H)*H'*channel_train_out;   
%             LS_QMF=temp.*LS_esti';

            y=convol(seq,h)';      % Channel output 
            r=y+sqrt(noisevar(i)/2)*randn(1,length(y)); % RXed signal
            
%             LS_s=convol(r,LS_QMF);
            LS_s=filter(1,LS_esti',r);
            
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
axis([0 13 10e-11 0.5]);

hnd=title('Equlaizer Performance Chart');
set(hnd,'fontsize',15);
hnd=xlabel(' Eb / No   in   dB ');
set(hnd,'fontsize',15);
hnd=ylabel('  Probability  of  Error (Pe)  ');
set(hnd,'fontsize',15);