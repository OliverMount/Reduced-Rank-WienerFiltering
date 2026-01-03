% OFDM with (QPSK) modulation and RR LS channel estimation

clear all;clc;close all;tic;
N=64; % No of sub-carriers
O=10; % no of OFDM symbols
len=N*O; % total number of QPSK symbols 

% Pilot positions in Comb type carriers
p_gap=3; %16 (3,5) 32-31 64(3,7)  512-7  256 -3
Pilot_pos=1:p_gap:N;
L=length(Pilot_pos); % For Comb type carriers

noiter=3;
Z=4; % Signifies QPSK
D=4; % Desired reduced rank

Eb=1; EbNodB=0:5:18; No=1./10.^(EbNodB/10);  Es=log2(Z)*Eb;
codebook=exp(j*2*pi*(0:Z-1)/Z + j*pi/4); % QPSK codebook

% Static Telephone Channel
% h=[0+j*0 0.0485+j*0.0194 0.0573+j*0.0253 0.0786+j*0.0282 0.0874+j*0.0447 0.9222+j*0.3031 0.1427+j*0.0349 0.0835+j*0.0157 0.0621+j*0.0078 0.0359+j*0.0049 0.0214+j*0.0019];
h=[1+j -.8+j .25-j];
% h=[1+j];
M=length(h);
H=CircMtx(h,N); % Channel matrix for AWGN

F=(1/sqrt(N))*dftmtx(N);
QQ=F'; % Q for Block type pilot carriers
QQQ1=QQ(Pilot_pos,1:M); % Q for Comb-type carriers
% QQQ=QQ(:,1:M); % For Block type carriers
QQQQ=QQ(Pilot_pos,Pilot_pos);

for loop=1:noiter,
  
 seq=randint(1,len,[0 Z-1]);  % generation of symbols      
 X_s=exp((j*2*pi*seq/Z)+(j*pi/4)); % Mapping of symbols
 X=SP(X_s,N); % serial to parallel conversion
 
 for i=1:length(No),
   for q=1:O, % upto  no of OFDM symbols
    temp=(2*randint(1,L)-1).';
    Pilot=(upsample(temp,p_gap));
    Pilot=flipud(Pilot(1:N));
    % Time domain Channel Estimation part
    noi=randn(1,N)+j*randn(1,N);
    noise=sqrt(No(i)/2)*noi.';
    x=F'*Pilot; % taking IFFT before transmission
    y1=H*x+noise; % Passing thro Channel  + noise (received symbols)
    y2=flipud(y1(N:-p_gap:1)); % Received values at the  
%     [H_capt,H_capt_RR]=KTFt_esti_comb(y,y2,Pilot,flipud(temp),Pilot_pos,M,D); E_capt=diag(H_capt);% Time doamin TF estimator Block type
     [H_capt,H_capt_rr,H_captrr]=KTFt_esti_comb(y1,y2,Pilot,flipud(temp),Pilot_pos,M,D,p_gap); ;% Time doamin TF estimator Block type
%     Y=F*y;
%     H_capf=TFf_esti(Y,Pilot);  E_capf=diag(H_capf);% Freq doamin TF estimator
    E_capt=diag(H_capt);
    E_capt_RR1=diag(H_capt_rr);% Time doamin TF estimator
    E_captRR=diag(H_captrr(:,3));% Time doamin TF estimator
    E=F*H*F'; % For Full channel knowledge
         
%         noi=randn(1,N)+j*randn(1,N);
%         noise=sqrt(No(i)/(2*log2(Z)))*noi.';
        x=F'*X(:,q); % taking IFFT before transmisssion
        y(:,q)=H*x+noise; % Passing thro Channel  + noise (received symbols)
        Y_p(:,q)=F*y(:,q); % FFT after reception
        Y_decoded(:,q)=E'*inv(E*E'+(No(i)/var(X_s))*eye(N))*Y_p(:,q);
        Y_decoded_capt(:,q)=E_capt'*inv(E_capt*E_capt'+(No(i)/var(X_s))*eye(N))*Y_p(:,q);
%         Y_decoded_capf(:,q)=E_capf'*inv(E_capf*E_capf'+(No(i)/var(X_s))*eye(N))*Y_p(:,q);
        Y_decoded_capt_RR1(:,q)=E_capt_RR1'*inv(E_capt_RR1*E_capt_RR1'+(No(i)/var(X_s))*eye(N))*Y_p(:,q); % Full rank solution
        Y_decoded_captRR(:,q)=E_captRR'*inv(E_captRR*E_captRR'+(No(i)/var(X_s))*eye(N))*Y_p(:,q); % RR solution
    end
    Y=PS(Y_decoded,N); % Parallel to serial conversion (with known channel) 
    Y_capt=PS(Y_decoded_capt,N); % With estiamted channel
%     Y_capf=PS(Y_decoded_capf,N); % With estimated channel
    Y_capt_RR1=PS(Y_decoded_capt_RR1,N); 
    Y_captRR=PS(Y_decoded_captRR,N); 

    for k=1:length(Y),
        [val pos]=min(abs(Y(k)-codebook));
        decoded(k)=codebook(pos);
        [val_1 pos_1]=min(abs(Y_capt(k)-codebook));
         decoded_capt(k)=codebook(pos_1);
%         [val_2 pos_2]=min(abs(Y_capf(k)-codebook));
%         decoded_capf(k)=codebook(pos_2);
        [val pos_rr1]=min(abs(Y_capt_RR1(k)-codebook));
        decoded_capt_rr1(k)=codebook(pos_rr1);
        
        [val posrr]=min(abs(Y_captRR(k)-codebook));
        decoded_captrr(k)=codebook(posrr);
        
    end    

     number=(decoded==X_s); %Counting the errors. 
     pe(loop,i)=(len-sum(number))/len;
     
     number1=(decoded_capt==X_s);
     pe_capt(loop,i)=(len-sum(number1))/len;
%      
%      number2=(decoded_capf==X_s);
%      pe_capf(loop,i)=(len-sum(number2))/len;

    number_rr1=(decoded_capt_rr1==X_s);  % Full rank 
     pe_capt_rr1(loop,i)=(len-sum(number_rr1))/len;
     
     numberrr=(decoded_captrr==X_s); % For Reduced rank
     pe_captrr(loop,i)=(len-sum(numberrr))/len;
end
fprintf('Finished %d Iteration for %d symbols in %f seconds  \n',loop,len,toc);
end

EBNO=10*log10(Eb./No);
ensemble_avg=mean(pe)
ensemble_avg_capt=mean(pe_capt)
% ensemble_avg_capf=mean(pe_capf);
ensemble_avg_capt_rr1=mean(pe_capt_rr1)
ensemble_avg_captrr=mean(pe_captrr)

% figure,hnd=semilogy(EBNO,erfc(sqrt((Eb./No))).*(1-.25*erfc(sqrt((Eb./No)))),'k-',EBNO,ensemble_avg,EBNO,ensemble_avg_capt,EBNO,ensemble_avg_capf);grid;

% figure,hnd=semilogy(EBNO,erfc(sqrt((Eb./No))).*(1-.25*erfc(sqrt((Eb./No)))),'k-',EBNO,ensemble_avg,EBNO,ensemble_avg_capt,EBNO,ensemble_avg_capt_rr1);grid;
% figure,hnd=semilogy(EBNO,ensemble_avg,'*-',EBNO,ensemble_avg_capt,'p-',EBNO,ensemble_avg_capt_rr1,'o-');grid;
figure,hnd=semilogy(EBNO,ensemble_avg,'*-',EBNO,ensemble_avg_capt,'p-',EBNO,ensemble_avg_capt_rr1,'o-',EBNO,ensemble_avg_captrr,'s-');grid;
set(hnd(1),'linewidth',1.5);
set(hnd(2),'linewidth',1.5);
set(hnd(3),'linewidth',1.5);
set(hnd(4),'linewidth',1.5);
% set(hnd(5),'linewidth',1.5);
% set(hnd(6),'linewidth',1.5);
%set(hnd(7),'linewidth',1.5);
legend('Perfect Channel','Time-domain Estimator','Comb-type Lowpass interpolated','Low rank');
axis([0 16 10e-5 1]);
hnd=xlabel(' Eb/No  in   dB ');
set(hnd,'fontsize',13,'color',[0 0 1]);
hnd=ylabel('  Probability  of  Error (Pe)  ');
set(hnd,'fontsize',13,'color',[0 0 1]);

% plot(abs([E,diag(F*H*F')])); legend('Estimated','Original');


% figure,stem(abs([diag(F*H*F') H_capt]));legend('Original','Estimated');