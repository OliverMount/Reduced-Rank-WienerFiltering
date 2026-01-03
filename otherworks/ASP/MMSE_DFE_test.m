% This program for MMSE equalization of Known channel (with delayed
% decision ) in a BPSK Transmission ( no unequlai comparison)

clc;clear all;close all;
tic;
noiter=25; % number of averaing iterations
len=1000; % length of the transmiited sequence

del=1;  % Decision delay
M=3;    %MMSE filter length (Why 3?)
Q=1; % order of the Feedback filter
Eb=1; % Energy per bit in V
EbNodB=0:1:12;
noisevar=1./10.^(EbNodB/10);

h=[1 0.5]; % Assumed ISI Channel model % Ali sayeed channel
% h=[.04 -.05 .07 -.21 -.5 .72 .36 0 .21 .03 .07]; % Proakis Channel A
% (Normal)
% h=[.407 .815 .407]; % Proakis Channel B ( Zero Close to unit circle)
% h=[.227 .460 .668 .460 .227]; % Proakis Channel C (Spectral Nulls) 
% h=[.5 1.2 1.5 -1];
% h=[1 1.2 1.5 -1];
L=length(h)-1; % ISI length
% figure,freqz(h);figure,zplane(h)

h_corr=corr(h,h)';
h_ori=[h_corr(L+1:end) zeros(1,M)]; % To eliminate the negative time lag 

r_dx_mtx=convolmtx(h,M); % Cross correlation vector 
R_DX=r_dx_mtx(del+1:del+Q,:)';
d_d(1,:)=zeros(1,Q); % Previously detected symbol at initial time 

for loop=1:noiter,   
    
    seq=randint(1,len);
    seq=2*seq-1; % Transmitted sequence
%     conv_x=convolmtx(seq,M);
%     conv_xt(1,:)=conv_x(1,:); % transformed observation at first time instant
    r_dx=r_dx_mtx(del+1,:);   % Change here the decision variable

        for i=1:length(noisevar),
            
            N=noisevar(i)*eye(M);  % Noise covariance matrix
            X=var(seq)*toeplitz(h_ori(1:M));  % Unit variance transmitted correlation
            Y=X+N; % Received Correlation matrix
            MMSE_w=inv(Y)*r_dx';  % Optimum MMSE solution
            
            y=convol(seq,h)';      % Channel output 
            r=y+sqrt(noisevar(i)/2)*randn(1,length(y)); % RXed signal
            s=convol(r,MMSE_w')';  % MMSE filter 
            rx = s > 0;  % Decision Box 
            rx1=(2*rx-1);
            rx_all=rx1(1+del:len+del); % L= ISI length
           
            number=(seq==rx_all); %Counting the errors. (Change if del is an array) 
            pe(loop,i)=(len-sum(number))/len;
            
            conv_x=convolmtx(r,M);
            conv_xt(1,:)=conv_x(1,:);
%           Finding DFE estimate
            R_t=Y-(1/var(seq))* R_DX*R_DX';
            g= inv(R_t)*r_dx'; % Feed foreward filter 
            % f=(1/var(seq))* R_DX'*g;  % Feedback filter
            
            % DFE Detection
            for z=2:len+M,
                d_cap(z-1)=conv_xt(z-1,:)*g; %  Estimation
                d_t(z-1)=2*(d_cap(z-1)>0)-1; % Decision
                d_d(z,:)=[d_t(z-1) d_d(z-1,1:end-1)];   % Detection
                conv_xt(z,:)=conv_x(z,:)-(1/var(seq))*d_d(z,:)*R_DX'; % Transformed observation
            end      
            
            rx_DFE=d_t(1+del:len+del); % L= ISI length
            number_DFE=(seq==rx_DFE); %Counting the errors. (Change if del is an array) 
            pe_DFE(loop,i)=(len-sum(number_DFE))/len;

        end
    fprintf('Finished %d Iteration for %d bits  in %f seconds  \n',loop,len,toc);
end

No=noisevar;
EBNO=10*log10(Eb./No);
ensemble_avg=mean(pe);
ensemble_avg_DFE=mean(pe_DFE);

figure,hnd=semilogy(EBNO,.5*erfc(sqrt((Eb./No))),'b',EBNO,(ensemble_avg),'ro-',EBNO,ensemble_avg_DFE);grid;
% set(hnd(1),'linewidth',2);
% set(hnd(2),'linewidth',2);
% legend('NO ISI','MMSE Equalized');
axis([0 12 10e-9 0.5]);

hnd=title('Equlaizer Performance Chart');
set(hnd,'fontsize',15);
hnd=xlabel(' Eb / No   in   dB ');
set(hnd,'fontsize',15);
hnd=ylabel('  Probability  of  Error (Pe)  ');
set(hnd,'fontsize',15);