%Root mean square error graph vs nhf
%AARYA DESAI
clear all;
clc
Error_multi=zeros(5,1);
Error_high=zeros(5,1);
No_of_hf=zeros(5,1);
coeff_determination_h=zeros(5,1);
coeff_determination_m=zeros(5,1);
for i=1:5
    if i==1
s=load('freq_mass_5')
    elseif i==2
s=load('freq_mass_7')
    elseif i==3
s=load('freq_mass_9')
    elseif i==4
s=load('freq_mass_10')
    elseif i==5
s=load('freq_mass_12')
    end        
     
Kp_LF = predict(s.model_LF,[s.TimeSampling]');
Kp_MF = predict(s.model_HF,[[s.TimeSampling]',Kp_LF(:)]);
Kp_HF = predict(s.model_HF2,[s.TimeSampling]');

%error calculation

no=length(s.TimeSampling) %no of data points

Error_m=s.Identified_mass_array+1-Kp_MF;
sq_m=Error_m.^2;              
MSE_m=sum(sq_m)/no;        %sample_variance
RMSE_m=sqrt(MSE_m)

Error_h=s.Identified_mass_array+1-Kp_HF;
sq_h=Error_h.^2;
MSE_h=sum(sq_h)/no;
RMSE_h=sqrt(MSE_h)

% calculation of coefficient of determination
mean_original_data=sum(s.Identified_mass_array+1)/no; %mean of actual system
deviation=s.Identified_mass_array+1-ones(no,1)*mean_original_data;
sq_deviation=deviation.^2;

coeff_determination_m(i,1)=1-sum(sq_m)/sum(sq_deviation)
coeff_determination_h(i,1)=1-sum(sq_h)/sum(sq_deviation)

if coeff_determination_h(i,1)<0
    coeff_determination_h(i,1)=0
end

if coeff_determination_m(i,1)<0
    coeff_determination_m(i,1)=0
end


No_of_hf(i,1)=s.nHF
Error_multi(i,1)=RMSE_m
Error_high(i,1)=RMSE_h
clear s;

    end

plot(No_of_hf,Error_multi)
xlabel('No of high fidility')
ylabel('RMSE OF multifidility')
xlim([5,12])

%figure;
%{
plot(No_of_hf,Error_high)
xlabel('No of high fidility')
ylabel('RMSE of only high fidility')
xlim([5,12])
%}



