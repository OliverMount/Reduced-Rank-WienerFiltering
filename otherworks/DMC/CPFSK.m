clear all;clc;close all; 
tic; 
Eb=.1; % Energy per bit in V 
% No=.1:.1:2; % Noise Spectral Density in Watts/ Hz; 
i=1:2:31; 
No=Eb./i; 

h=0.715;  % Modulation Index 
% h=0.5; 
Ns=2;   % Detector Observation interval 
te=2^(Ns-1); 

% M=2; T=1e-3; 
% fd=h/(2*T); % Peak frequency deviation 
% ts=0.1*T;fs=1/ts; 
% fc=1/(10*ts); 
% E=Eb*log2(M); 
% Tb=T/log2(M); 

M=2;T=1e-3; 
fc=5/T; 
fs=10*fc;ts=1/fs; 
E=Eb*log2(M); 
Tb=T/log2(M); 

noiter=10; 

% Frequency g and Phase q pulse genration 
t=0:ts:T;  % increment is the sampling time here 
for i=1:length(t), 
   g(i)=1/(2*T); 
end 
for j=1:length(t), 
   q(j)=integral(g(1:j),ts); 
end 
q(1)=0; 

N=length(q); 

zzz=MF(M,Ns); 
clear i j; 
%################################################################## 
for loop=1:noiter, 

seq=randint(1,2500); 
% seq=[0 0 0 0 0]; 
% seq=[1 0 1 0 1]; 
% seq=[1 1 0 0 1 1 0]; 
% seq=[1 -1 -1  1  1  1 -1 -1  1 -1]; 
% seq=[-1  1 -1 -1 -1 -1 -1  1 -1 -1]; 
% seq=[zeros(1,25) ones(1,25)]; 
seq1=2*seq-1; 
k=length(seq1); 
% Generation of transmission signal 
phase(1)=0;    % Initial Phase at the transmitter 

j=1;m=1; 
for i=1:k, 
   pulse=2*pi*h*seq1(i)*q; 
   phi(1,j:i*N)=phase(i)+pulse; 
   phase(i+1)=phi(1,end); 
   j=j+N; 
end 
% figure,plot(phi); 
% clear i j; 
% j=1; 
% for i=1:k-1, 
%     new_phi(1,j:i*(N-1))=phi(1,(i*N)+2:((i+1)*N)); 
%     j=j+(N-1); 
% end 
% final_phi=[seq1(1)*pulse new_phi]; 

for l=1:length(No), 
n=0:(length(phi)-1); 
val=2*pi*fc*n*ts; 

s(l,:)=sqrt(2*E/T)*cos(val+phi); % output of the transmitter 
% Channel Part 
noisevar=No(l)*(fs/2); 
noi=randn(1,length(s(l,:))); 
noise=sqrt(noisevar)*noi; 
r(l,:)=s(l,:)+noise; 
% figure,plot(s); 
% figure,plot(r); 

% clear val n; 
clear n; 
% n=0:(Ns*N)-1; 
% val=2*pi*fc*n*ts; 
% disp(size(val));disp(size(n)); 

% clear i val ; 
clear j; 
% RECEIVER part 
j=1;pp=1; 
% Matched filtering 
for i=1:(k-Ns+1), 

   temp=r(l,j:(Ns+i-1)*N); 

       for p=1:te, 
          
ss(p,:)=sqrt(2*E/T)*cos(val(pp:((Ns+i-1)*N))+prototype(zzz(p,:),q,h,phase(i))); 
          
ss(p+te,:)=sqrt(2*E/T)*cos(val(pp:((Ns+i-1)*N))+prototype(zzz(p+te,:),q,h,phase(i))); 
       end 

       for m=1:(2*te), 
           x(i,m)=inteexp(temp,ss(m,:),ts,noisevar); 
       end 

       j=j+N;pp=pp+N; 
     clear p; 
end 

clear i j m val; 

for i=1:(k-Ns+1), 
   lamda1(i)=sum(x(i,1:te)); 
   lamda2(i)=sum(x(i,te+1:end)); 
end 

le=lamda1-lamda2; 
d=(le<0); 
% d=lamda1<lamda2; 

e=2*d-1; 
val=(seq1(1:k-(Ns-1))==e); 

perror(loop,l)=(k-(Ns-1)-sum(val))/(k-(Ns-1)); 

end 
% clear phase phi; 
clear phi; 

fprintf('Finished %d Iteration for %d bits  in %f seconds  \n',loop,k,toc); 
end 
% 
finalsum=mean(perror); 

EBNO=10*log10(Eb./No); 
% figure,plot(EBNO,pe,'bh'); 
figure,hnd=semilogy(EBNO,finalsum,'k');grid; 
set(hnd,'linewidth',2); 
axis([0 14 10e-7 2]); 

hnd=title('CPFSK Performance Chart'); 
set(hnd,'fontsize',15); 
hnd=xlabel(' Eb / No   in   dB '); 
set(hnd,'fontsize',15); 
hnd=ylabel('  Probability  of  Error (Pe)  '); 
set(hnd,'fontsize',15); 
