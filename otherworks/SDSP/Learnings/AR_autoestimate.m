% Program for Finding Autocorrletaion of an  AR process
%-------------------------------------------------------------
clear all;clc;close all;
% Exmaple 1 AR(2) process
a=[1 0 .81];len=1000;
noisevar=1;
x=AutoR(a,noisevar,len);
r=TimeAC(x);

% r_th=[1 0 -.81 0 .6561 0 -0.5314 0   0.4305 0  -0.3487 0   0.2824 0  -0.2288 0   0.1853 0  -0.1501 0 0.1216 0  -0.0985 0];
r_th=[2.9078 0 -2.355 0  1.9076 0 -1.5451 0   1.2515 0   -1.0138  0  0.8211  0  -0.6651 0  0.5387 0 -0.4364 0  0.3535 0  -0.2863 0];

figure,stem(r(1:24),'r*');hold;
% theo=[r_th zeros(1,len-length(r_th))];
theo=r_th;
stem(theo(1:24),'kp');
title(' AR correlation sequence ');