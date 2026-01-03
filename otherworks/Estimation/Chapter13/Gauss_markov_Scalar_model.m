% Typical realization of Gauss Markov Model
clear all;clc;close all;
len=100; % Length of the input sequence
a=.98;
var_u=0.01; % Process Noise variance

s=AutoR_Kalman([1 -a],var_u,len); % Process Model (AR(1))

figure,plot(s);
axis([0 100 0 10]);