% Program for Finding Autocorrletaion of an  MA process
%-------------------------------------------------------------
clear all;clc;close all;
% Exmaple 1
h=[1 .5];  % Process parameters
noisevar=1;
% Theretical  Autocorrelation by Yule walker method 
len=1000;
h_th=[1.25 .5];
x=MovA(h,noisevar,len);
rx=TimeAC(x);
figure,stem(rx(1:50),'r*');hold;
theo=[h_th zeros(1,len-length(h_th))];
stem(theo(1:50));
title(' MA correlation sequence ');

clear all;
%_______________________________________________________________

% Example 2

h=[-.5 1 2 -2.5];  % Process parameters
len=1000;noisevar=1;
% Theretical  Autocorrelation by Yule walker method 
h_th=[11.5 -3.5 -3.5 1.25];
x=MovA(h,noisevar,len);
rx=TimeAC(x);
figure,stem(rx(1:50),'r*');hold;
theo=[h_th zeros(1,len-length(h_th))];
stem(theo(1:50));
title(' MA correlation sequence ');