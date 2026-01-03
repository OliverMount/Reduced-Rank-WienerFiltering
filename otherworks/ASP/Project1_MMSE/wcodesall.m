clear all;
clc;

[l h]=wfilters('bior1.3','r');
% N=input('Enter  the Number of Stages');
N=5;
for loop=1:32,
    
user=loop;
data=zeros(1,2^N);
data(user)=1;
% data=2*data-1;
% figure,plot(data);

p=2^N; ph=p/2;
% First stage computation
for j=1:ph,
    
    k(j,:)=cconv(upsample(data(2*j-1),2)',l);
    m(j,:)=cconv(upsample(data(2*j),2)',h);
    op(j,:)=k(j,:)+m(j,:);
end
clear k m j ph;

% % Second stage computation
ph=p/4;
for j=1:ph,
    
    k(j,:)=cconv(upsample(op(2*j-1,:),2),l);
    m(j,:)=cconv(upsample(op(2*j,:),2),h);
    op1(j,:)=k(j,:)+m(j,:);
end
clear k m j ph;
% % 
% Third Stage Computation
ph=p/8;
for j=1:ph,
    
    k(j,:)=cconv(upsample(op1(2*j-1,:),2),l);
    m(j,:)=cconv(upsample(op1(2*j,:),2),h);
    op2(j,:)=k(j,:)+m(j,:);
end
% 
clear k m j ph;
% 
% % Fourth stage computation
ph=p/16;
for j=1:ph,
    
    k(j,:)=cconv(upsample(op2(2*j-1,:),2),l);
    m(j,:)=cconv(upsample(op2(2*j,:),2),h);
    op3(j,:)=k(j,:)+m(j,:);
end

clear k m j ph;
% % Last stage computation
ph=p/32;
for j=1:ph,
    
    k(j,:)=cconv(upsample(op3(2*j-1,:),2),l);
    m(j,:)=cconv(upsample(op3(2*j,:),2),h);
    op4(j,:)=k(j,:)+m(j,:);
end

[val loc]=max(data);
ws(loop,:)=real(op4);
clear k m j ph;
clear data;
end
save wavcodes.mat ws;

% figure,plot(gen);
% axis([0 64 -1 1]);
% str=sprintf(' Wavelet code for the data bit at %dth   filter input', loc);
% hnd=title(str);
% set(hnd,'fontsize',15);
% hnd=xlabel(' Time in Seconds ');
% set(hnd,'fontsize',15);
% hnd=ylabel(' Amplitude in Volts ');
% set(hnd,'fontsize',15);

% figure,periodogram(gen);