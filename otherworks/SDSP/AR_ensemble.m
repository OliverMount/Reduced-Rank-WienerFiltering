% An ensemble of AR(1) processes
clc;clear all;close all;
k=200;iter=100;
% figure,plot(d);

for i=1:iter,

    noi=randn(1,k);
    v1(i,:)=AutoRnoise([ 1 -.8]',noi,k); % AR(1) process

end