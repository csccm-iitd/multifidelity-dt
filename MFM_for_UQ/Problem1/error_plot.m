clear all
clc
s=load('50_9.mat')
e(1)=s.RMS
clear s
s=load('50_10.mat')
e(2)=s.RMS
clear s

s=load('50_12.mat')
e(3)=s.RMS
clear s


s=load('50_16.mat')
e(4)=s.RMS
clear s
x=[9 10 12 16]
plot(x,e)



