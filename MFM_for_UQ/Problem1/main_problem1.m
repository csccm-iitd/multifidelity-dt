% Aarya Sheetal Desai
% Problem 1 - Linear correlation

clc;
clear;

nsamp_LF = 50;
nsamp_HF = 16;

lb = 0.0;
ub = 1.0;

xsamp_LF = linspace(lb,ub,nsamp_LF);
ysamp_LF = sin(8*pi*xsamp_LF);

xsamp_HF = linspace(lb,ub,nsamp_HF);
ysamp_HF = (xsamp_HF - sqrt(2)).*(sin(8*pi*xsamp_HF)).^2;

[dmodel_lf, dmodel_hf] = MF_HPCFE(xsamp_LF', ysamp_LF', xsamp_HF', ysamp_HF');

xsimu = linspace(lb,ub,1000);

ysimu = MF_predict(xsimu',dmodel_lf,dmodel_hf);

ysimu_act = (xsimu - sqrt(2)).*(sin(8*pi*xsimu)).^2;

% Train from HF data

theta0 = 10;
lob = 0.1;
upb = 20;

cd('./H-PCFE_p')

dmodel_pk_lf = PK_fit(xsamp_HF', ysamp_HF, @PCFE, @corrgauss, theta0, lob, upb);
ysimu_HF = predictor1(xsimu',dmodel_pk_lf);

cd('../')


theta0 = 10;
lob = 0.1;
upb = 20;

cd('./H-PCFE_p')

dmodel_pk_hf = PK_fit(xsamp_LF', ysamp_LF, @PCFE, @corrgauss, theta0, lob, upb);
ysimu_LF = predictor1(xsimu',dmodel_pk_hf);

cd('../')

figure(1)
plot(xsimu', ysimu_act,'r*', xsimu', ysimu','g', xsimu', ysimu_LF,'b', xsimu', ysimu_HF,'c')
legend('Actual system','MF-HPCFE','LF-HPCFE','HF-HPCFE')


figure(2)
x=linspace(0,1,1000);
y = sin(8*pi*x);
z = (x - sqrt(2)).*(sin(8*pi*x)).^2;

plot(x,y,'r',x,z,'g',xsamp_LF,ysamp_LF,'o',xsamp_HF,ysamp_HF,'*')
legend('g_{low}','g_{high}','low fidelity training data (50 points)','high fidelity training data (16 points)')
ylim([-2 2])

RMS=sqrt(mean((ysimu_act'-ysimu).^2));

xsimu=rand(1,10000);
ysimu = MF_predict(xsimu',dmodel_lf,dmodel_hf);
ysimu_LF = predictor1(xsimu',dmodel_pk_lf);
ysimu_HF = predictor1(xsimu',dmodel_pk_hf);

ysimu_act = (xsimu - sqrt(2)).*(sin(8*pi*xsimu)).^2;

figure(3)
[f,xi]=ksdensity(ysimu_act);
plot(xi,f)
hold on;
[f,xi]=ksdensity(ysimu);
plot(xi,f)
hold on;
[f,xi]=ksdensity(ysimu_HF);
plot(xi,f)
hold on;
[f,xi]=ksdensity(ysimu_LF);
plot(xi,f)
legend('Actual','MF-DT','GP-DT_(high)','GP-DT_(low)')
hold off;
