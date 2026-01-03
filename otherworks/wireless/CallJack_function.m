% simulation of Jakes Flat fading  for M waveforms
% Maximum we can generate M waveforms by this method

function r=CallJack_function(P,M,fm,fs)

N=2*((2*M)+1);
n=1:M;
f_n=fm*cos(2*pi*n/N)/fs;
alpha=0;beta=(pi*n/(M)); gamm=2*pi*(0:M-1)'*n/M;
gamm_c=cos(gamm); gamm_s=sin(gamm);

ideal=besselj(0,2*pi*fm*(0:P-1)/fs);
% disp(ideal);

A=cos(2*pi*(0:P-1)'*f_n); B=sin(2*pi*(0:P-1)'*f_n);
b=2*cos(pi*(n)/M); c=2*sin(pi*(n)/M);
for q=1:M, % M waveforms
r_I(:,q)=A*(gamm_c(q,:).*b)'+B*(gamm_s(q,:).*(-b))'+(sqrt(2)*cos(2*pi*(fm/fs)*(0:P-1)))';
r_Q(:,q)=A*(gamm_c(q,:).*c)'+B*(gamm_s(q,:).*(-c))';
end
r=(1/sqrt((2*M)+1))*((r_I)+j*(r_Q)); % Low-pass equivalent time correlated multipath samples
% disp(size(r)):
% rrr=r(:,1).';
% new_AC=TimeAC(rrr);%/var(rrr);
% % disp(abs(new_AC));
% 
% kkk=fm*(0:P-1)/fs;
% figure,plot(kkk,real(new_AC),kkk,ideal);
% legend('Simulated ','Ideal');
% axis([0 14 -1 1.1]);
% % 
% figure,semilogy(abs(r(:,1)));
% figure,semilogy(abs(r(:,2)));
% figure,semilogy(abs(r(:,3)));
% figure,semilogy(abs(r(:,4)));
% figure,semilogy(abs(r(:,5)));
% figure,semilogy(abs(r(:,6)));
