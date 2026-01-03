% OFDM with (QPSK) modulation and RR LS channel estimation

clear all;clc;close all;tic;
N=512; % No of sub-carriers
O=10; % no of OFDM symbols
len=N*O; % total number of QPSK symbols 

% Pilot positions in Comb type carriers
p_gap=7; %16 (3,5) 32-31 64(3,7)  256 -3 512-7  
Pilot_pos=1:p_gap:N;
L=length(Pilot_pos); % For Comb type carriers

noiter=20;
Z=4; % Signifies QPSK
D=4; % Desired reduced rank

Eb=1; EbNodB=10:3:40; No=1./10.^(EbNodB/10);  Es=log2(Z)*Eb;
codebook=exp(j*2*pi*(0:Z-1)/Z + j*pi/4); % QPSK codebook

% Static Telephone Channel
% h=[0+j*0 0.0485+j*0.0194 0.0573+j*0.0253 0.0786+j*0.0282 0.0874+j*0.0447 0.9222+j*0.3031 0.1427+j*0.0349 0.0835+j*0.0157 0.0621+j*0.0078 0.0359+j*0.0049 0.0214+j*0.0019];
% h=[1+j -.8+j .25-j];
% Time varyimg mobile channel
Pow_delay=[-3 -5 -8 -10  -12 -11 -8 -10 -12 -9 -9.23]; % Power delay profile values in dB
h_temp=CallJack_function(O,11); % this is 'capital O' not zero
for xx=1:O,
h(:,xx)=(diag(10.^(Pow_delay./10))*((h_temp(xx,:).')/norm(h_temp(xx,:))));
end

% M=length(h);
M=11; % serious parameter..
% H=CircMtx(h,N); 

F=(1/sqrt(N))*dftmtx(N);
% QQ=F'; % Q for Block type pilot carriers
% QQQ1=QQ(Pilot_pos,1:M); % Q for Comb-type carriers
% % QQQ=QQ(:,1:M); % For Block type carriers
% QQQQ=QQ(Pilot_pos,Pilot_pos);

for loop=1:noiter,
  
 seq=randint(1,len,[0 Z-1]);  % generation of symbols      
 X_s=exp((j*2*pi*seq/Z)+(j*pi/4)); % Mapping of symbols
 X=SP(X_s,N); % serial to parallel conversion
 
 for i=1:length(No),
   for q=1:O, % upto  no of OFDM symbols
     H=CircMtx(h(:,q).',N); % channel matrix for Frequency selective Fading
    temp=(2*randint(1,L)-1).';
    Pilot=(upsample(temp,p_gap));
    Pilot=flipud(Pilot(1:N));
    % Time domain Channel Estimation part
    noi=randn(1,N)+j*randn(1,N);
    noise=sqrt(No(i)/2)*noi.';
    x=F'*Pilot; % taking IFFT before transmission
    y1=H*x+noise; % Passing thro Channel  + noise (received symbols)
%     Y1=F*y1;
    y2=flipud(y1(N:-p_gap:1)); % Received values at the  
%     [H_capt,H_capt_RR]=KTFt_esti_comb(y,y2,Pilot,flipud(temp),Pilot_pos,M,D); E_capt=diag(H_capt);% Time doamin TF estimator Block type
     [H_capt,H_capt_rr,H_captrr]=KTFt_esti_comb(y1,y2,Pilot,flipud(temp),Pilot_pos,M,D,p_gap); ;% Time doamin TF estimator Block type
%     Y=F*y;
%     H_capf=TFf_esti(Y1,Pilot);  ;% Freq doamin TF estimator
    % Estimated TF ready for decoding
%     E_capt=diag(H_capt);
%     E_capf=diag(H_capf);
%     E_capt_RR1=diag(H_capt_rr);% Time doamin TF estimator
    E_captRR=diag(H_captrr(:,1));% RR  TF estimator
%     E_captRR3=diag(H_captrr(:,3));% RR TF estimator
%     E=F*H*F'; % For Full channel knowledge
         
%         noi=randn(1,N)+j*randn(1,N);
%         noise=sqrt(No(i)/(2*log2(Z)))*noi.';
        x=F'*X(:,q); % taking IFFT before transmisssion
        y(:,q)=H*x+noise; % Passing thro Channel  + noise (received symbols)
        Y_p(:,q)=F*y(:,q); % FFT after reception
        
        
%         Y_decoded(:,q)=E'*inv(E*E'+(No(i)/var(X_s))*eye(N))*Y_p(:,q);
%         Y_decoded_capt(:,q)=E_capt'*inv(E_capt*E_capt'+(No(i)/var(X_s))*eye(N))*Y_p(:,q);
%         Y_decoded_capf(:,q)=E_capf'*inv(E_capf*E_capf'+(No(i)/var(X_s))*eye(N))*Y_p(:,q);
%         Y_decoded_capt_RR1(:,q)=E_capt_RR1'*inv(E_capt_RR1*E_capt_RR1'+(No(i)/var(X_s))*eye(N))*Y_p(:,q); % Full rank solution
        Y_decoded_captRR(:,q)=E_captRR'*inv(E_captRR*E_captRR'+(No(i)/var(X_s))*eye(N))*Y_p(:,q); % RR solution
%         Y_decoded_captRR3(:,q)=E_captRR3'*inv(E_captRR3*E_captRR3'+(No(i)/var(X_s))*eye(N))*Y_p(:,q); % RR solution
        
        
    end
%     Y=PS(Y_decoded,N); % Parallel to serial conversion (with known channel) 
%     Y_capt=PS(Y_decoded_capt,N); % With estiamted channel
%     Y_capf=PS(Y_decoded_capf,N); % With estimated channel
%     Y_capt_RR1=PS(Y_decoded_capt_RR1,N); 
    Y_captRR=PS(Y_decoded_captRR,N); 
%      Y_captRR3=PS(Y_decoded_captRR3,N); 

    for k=1:length(Y_captRR),
        
%         [val pos]=min(abs(Y(k)-codebook));
%         decoded(k)=codebook(pos);
        
%         [val_1 pos_1]=min(abs(Y_capt(k)-codebook));
%          decoded_capt(k)=codebook(pos_1);
         
%         [val_2 pos_2]=min(abs(Y_capf(k)-codebook)); % For frequency domain estimator
%         decoded_capf(k)=codebook(pos_2);

%         [val pos_rr1]=min(abs(Y_capt_RR1(k)-codebook));
%         decoded_capt_rr1(k)=codebook(pos_rr1);
        
        [val posrr]=min(abs(Y_captRR(k)-codebook));
        decoded_captrr(k)=codebook(posrr);
        
%         [val posrr3]=min(abs(Y_captRR3(k)-codebook));
%         decoded_captrr3(k)=codebook(posrr3);
        
    end    

%      number=(decoded==X_s); %Counting the errors. 
%      pe(loop,i)=(len-sum(number))/len;
     
%      number1=(decoded_capt==X_s);
%      pe_capt(loop,i)=(len-sum(number1))/len;
     
%      number2=(decoded_capf==X_s);
%      pe_capf(loop,i)=(len-sum(number2))/len;

%     number_rr1=(decoded_capt_rr1==X_s);  % Full rank 
%      pe_capt_rr1(loop,i)=(len-sum(number_rr1))/len;
     
     numberrr=(decoded_captrr==X_s); % For Reduced rank
     pe_captrr(loop,i)=(len-sum(numberrr))/len;
     
%      numberrr3=(decoded_captrr3==X_s); % For Reduced rank
%      pe_captrr3(loop,i)=(len-sum(numberrr3))/len;
end
fprintf('Finished %d Iteration for %d symbols in %f seconds  \n',loop,len,toc);
end

EBNO=10*log10(Eb./No);
% ensemble_avg=mean(pe)
% ensemble_avg_capt=mean(pe_capt)
% ensemble_avg_capf=mean(pe_capf)
% ensemble_avg_capt_rr1=mean(pe_capt_rr1) % Fullrank
ensemble_avg_captrr=mean(pe_captrr) % RAnk 4
ensemble_avg_captrr3=mean(per_captrr3) % RAnk 3

% figure,hnd=semilogy(EBNO,erfc(sqrt((Eb./No))).*(1-.25*erfc(sqrt((Eb./No)))),'k-',EBNO,ensemble_avg,EBNO,ensemble_avg_capt,EBNO,ensemble_avg_capf);grid;

% figure,hnd=semilogy(EBNO,erfc(sqrt((Eb./No))).*(1-.25*erfc(sqrt((Eb./No)))),'k-',EBNO,ensemble_avg,EBNO,ensemble_avg_capt,EBNO,ensemble_avg_capt_rr1);grid;
% figure,hnd=semilogy(EBNO,ensemble_avg,'*-',EBNO,ensemble_avg_capt,'p-',EBNO,ensemble_avg_capt_rr1,'o-');grid;
% figure,hnd=semilogy(EBNO,ensemble_avg,'*-',EBNO,ensemble_avg_capt,'p-',EBNO,ensemble_avg_capt_rr1,'k0--',EBNO,ensemble_avg_captrr,'ks-');grid;


figure,hnd=semilogy(EBNO,ensemble_avg,'ko-',EBNO,ensemble_avg_capt_rr1,'k-',EBNO,ensemble_avg_captrr1,'k+-',EBNO,ensemble_avg_captrr2,'k*-',EBNO,ensemble_avg_captrr3,'kp-');grid;
legend('Perfect Channel','Full-rank','Rank 1','Rank 2','Rank 3');
axis([10 40 10e-5 1]);
hnd=xlabel(' SNR  in  dB ');
set(hnd,'fontsize',13,'color',[0 0 1]);
hnd=ylabel('  Probability  of  Error (Pe)  ');
set(hnd,'fontsize',13,'color',[0 0 1]);
% set(hnd(1),'linewidth',1.5);
% set(hnd(2),'linewidth',1.5);
% set(hnd(3),'linewidth',1.5);
% set(hnd(4),'linewidth',1.5);
% set(hnd(5),'linewidth',1.5);
% set(hnd(6),'linewidth',1.5);
%set(hnd(7),'linewidth',1.5);
% legend('Perfect Channel','Time-domain Estimator','Comb-type Lowpass interpolated','Rank 4');
% figure,hnd=semilogy(EBNO,ensemble_avg,'ko-',EBNO,ensemble_avg_capt_rr1,'k-',EBNO,ensemble_avg_captrr,'k+-',EBNO,ensemble_avg_captrr3,'k*-');grid;
% legend('Perfect Channel','Full-rank','Rank 2','Rank 3');
% % set(hnd(1),'linewidth',1.5);
% % set(hnd(2),'linewidth',1.5);
% % set(hnd(3),'linewidth',1.5);
% % set(hnd(4),'linewidth',1.5);
% % legend('Perfect Channel','Comb-type Time domain Full-rank','Rank 4','Rank 3','Freq domain Full rank');

% plot(abs([E,diag(F*H*F')])); legend('Estimated','Original');


% figure,stem(abs([diag(F*H*F') H_capt]));legend('Original','Estimated');
% 
ensemble_avg =[   0.2707    0.1683    0.0999    0.0598    0.0332    0.0152    0.0068    0.0026    0.0010    0.0005    0.0002];
ensemble_avg_capt_rr1 = [  0.4216    0.2834    0.1832    0.1051    0.0635    0.0342    0.0174    0.0059    0.0026    0.0012    0.0004];
ensemble_avg_captrr2 = [ 0.4203    0.2827    0.1827    0.1062    0.0643    0.0347    0.0181    0.0072    0.0037    0.0019    0.0010];
 ensemble_avg_captrr3 =[  0.4215    0.2832    0.1831    0.1052    0.0636  0.0342    0.0175    0.0060    0.0028    0.0013    0.0003];
ensemble_avg_captrr1 =[    0.4229    0.2962    0.1955    0.1261    0.0859    0.0628    0.0489    0.0396    0.0387    0.0357    0.0343];