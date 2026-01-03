function [r]=fade_gen(N,fd,fs)%,fade_real,fade_imag); % 
%fd=1;%input("enter fade rate fd : ");

%Tsym=100e-6;%try 10e-6 (10 micro seconds)

%fs=1/Tsym;%can be twice or higher times symbol rate 
%fs >= 2*fm

%N1=128;%samples within the Doppler PSD range N1/2 on either side

%N=1e5;

% => (N1/2)*fs/N (approx) < fd
% => N (approx)>(N1*fs)/(2*fd)
% => N=ceil((N1*fs)/(2*fd))

Nsamp=1e4;
Nsamp(Nsamp>N)=N;

N1=floor((2*N*fd)/fs);
niter=ceil(N/Nsamp);

Nend=N;
if (N1<20 & fd<=100)
    N=(fs)/2;
    N1=ceil((2*2*N*fd)/fs);
    N1=20;
    N=N1*fs/(2*fd);
elseif (N1<20 & fd>100)
    % N=(fs)/6.6;N1=floor((2*N*fd)/fs); 
    N1=100;
    N=N1*fs/(2*fd);
end;
N1=100;
N=N1*fs/(2*fd);

% niter=ceil(N/Nsamp);

% Nend=N;

% N=ceil(N1*fs/(2*fd))

% N=1000;

f=(0:N1/2)*(fs/N);

psd=zeros(1,N1/2+1);

for m=1:N1/2
    psd(1,m)= 1./(sqrt(1-(f(1,m)/fd)^2));%jakes model
end

norm=(sum(psd)+sum(psd(2:N1/2+1)));%0.5*pi*2.5*
psd=psd./norm;
%plot(abs(psd));
%keyboard;
d=[];

a1=randn(1,((N1/2)+1));
a1(1,1)=0;
a2=randn(1,N1/2);

a3=randn(1,N1/2+1);
a3(1,1)=0;
a4=randn(1,N1/2);

b1=(a1+i*[0 a2])./2;
b2=(a3+i*[0 a4])./2;

b1=b1.*sqrt(psd);
b2=b2.*sqrt(psd);



%s1=fopen(fade_real,'w');
%s2=fopen(fade_imag,'w');


for m=1:niter-1
    
    d1=[ b1+i*b2 conj(b1(1,N1/2+1:-1:2))+i*conj(b2(1,N1/2+1:-1:2))]*[exp(i*2*pi*((0:N1/2)'*((m-1)*Nsamp:m*Nsamp-1))/N);exp(i*2*pi*(N-(N1/2):(N-1))'*((m-1)*Nsamp:m*Nsamp-1)/N)];
    %fprintf(s1,'%6.6f\n',real(d1));
    %fprintf(s2,'%6.6f\n',imag(d1));
    
    d=[d d1];
end;

m=niter-1;
d1=[ b1+i*b2  conj(b1(1,N1/2+1:-1:2))+i*conj(b2(1,N1/2+1:-1:2))]*[exp(i*2*pi*((0:N1/2)'*((m)*Nsamp:Nend-1))/N);exp(i*2*pi*(N-(N1/2):(N-1))'*((m)*Nsamp:Nend-1)/N)];

%fprintf(s1,'%6.6f\n',real(d1));
%fprintf(s2,'%6.6f\n',imag(d1));

r=[d d1];
clear d;
clear d1;

%fclose(s1);
%fclose(s2);

%cor=xcorr(abs(d),0)/(N*sum(psd)*fd);
fm=fd;
%rvar=d*d'/N;
%rvar
ideal=besselj(0,2*pi*fm*(0:N-1)/fs);
temp1=TimeAC(real(r));
temp2=TimeAC(imag(r));
% temp1_n=temp1/(PP/2);
% temp2_n=temp2/(PP/2);
temp1_n=temp1/(mean(real(r).^2));
temp2_n=temp2/(mean(imag(r).^2));

kkk=fm*(0:N-1)/fs;
% figure,plot(kkk,temp1_n,kkk,temp2_n,kkk,ideal);
% legend('Simulated I-phase','Simulated Q-phase','Ideal');
% axis([0 10 -1 1]);

figure,plot(kkk,temp1_n,kkk,ideal);
legend('Simulated I-phase','Ideal');
axis([0 10 -1 1]);

figure,plot(kkk,temp2_n,kkk,ideal);
legend('Simulated Q-phase','Ideal');
axis([0 10 -1 1]);%