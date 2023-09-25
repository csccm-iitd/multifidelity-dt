clear all; close all; clc; clf
fname='times';
fsize=16;


nf=3001;
f_max=2;
omega=linspace(0,f_max,nf);
zeta_0=0.02

X0_om=1./sqrt( (1-omega.^2).^2 + 4*zeta_0^2*omega.^2);
Delta_k=-0.25;
Xk_om=1./sqrt( (1+Delta_k-omega.^2).^2 + 4*zeta_0^2*omega.^2);

% Maximum points
Delta1=0
Hmax_p0= sqrt(1+Delta1-2*zeta_0^2)
H0_max=1/ (2 * zeta_0 * sqrt(1-zeta_0^2+Delta1))

Delta1=Delta_k
Hmax_pk= sqrt(1+Delta1-2*zeta_0^2)
Hk_max=1/ (2 * zeta_0 * sqrt(1-zeta_0^2+Delta1))



%% Maxima of the freqnecy response polt and analysis
figure(1);clf;hold all;
plot(omega,X0_om,'-k','linewidth',1.5);
x_p=[ Hmax_p0, Hmax_p0]
y_p=[0,H0_max]
plot(x_p,y_p,'*:k','MarkerSize',10,'HandleVisibility','off')
    
plot(omega,Xk_om,'--r','linewidth',1.5);
x_p=[ Hmax_pk, Hmax_pk]
y_p=[0,Hk_max]
plot(x_p,y_p,'*:r','MarkerSize',10,'HandleVisibility','off')

h1=legend(['Nominal',sprintf('\n'),'system'],['\Delta_k (t_s)='  num2str(Delta_k)]);
set(h1,'FontName',fname,'FontSize',fsize,'Box','off','Location','northwest');

% Texts
text(1.1,10,'frequency shift','FontName',fname,'fontsize',fsize)
text(1.13782505910165,26.3444227005871,'peak-response shift','FontName',fname,'fontsize',fsize)

text(1.05,20,'${\mathcal H}_0$','FontName',fname,'fontsize',20,'color','k','Interpreter','latex','HandleVisibility','off')
text(0.64,24,'${\mathcal H}_m$','FontName',fname,'fontsize',20,'color','r','Interpreter','latex','HandleVisibility','off')

% Create arrow
annotation('arrow',[0.586907449209932 0.489841986455982],...
    [0.402547671840355 0.321507760532151]);

% Create line
annotation('line',[0.4785 0.5056],[0.9147 0.9147]);

% Create arrow
annotation('arrow',[0.551918735891648 0.498871331828442],...
    [0.821616407982262 0.83370288248337]);

% Create line
annotation('line',[0.479721001221001 0.506821001221001],...
    [0.789397926634769 0.789397926634769]);

% Create doublearrow
annotation('doublearrow',[0.4943 0.493284493284493],...
    [0.9147 0.786283891547049]);

% Create doublearrow
annotation('doublearrow',[0.4654522803507 0.515262515262515],...
    [0.318 0.317384370015949]);

xlabel('Normalised frequency: \Omega = \omega/\omega_0','FontName',fname,'fontsize',fsize)
ylabel('Normalised respone: |U(\Omega)|/U_{st}','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');

%print -depsc figs/freqd_response_stiffness.eps
%print -djpeg figs/freqd_response_stiffness.jpeg




%% Define the stiffness chgane function
t1=linspace(0,1000,1000);
alpha1=0.4e-3;
beta1=2e-1;
eps1=0.01;
StiffnessDegradeFunction=exp(-alpha1*t1).*[1+eps1*(cos(beta1*t1))]/(1+eps1);


% Assume data is avialable at a lower smapling rate
nsamp=500;
TimeSampling=linspace(3,1000,nsamp);


%% Frequency-domain response
stiffness_array=[1,interp1(t1,StiffnessDegradeFunction,TimeSampling)];
nf=5001;
f_max=2;
omega=linspace(0.0,f_max,nf);
zeta_0=0.02

for j=1:nsamp+1
    X_k_om(:,j)=1./sqrt( (stiffness_array(j)-omega.^2).^2 + 4*zeta_0^2*omega.^2);  
end

% Plot the response cloud
figure(4);clf;hold all
plot(omega,X_k_om(:,1),'-k','linewidth',1.5)
plot(omega,X_k_om(:,2:3:nsamp+1),'-y','linewidth',0.5)
plot(omega,X_k_om(:,1),'-k','linewidth',1.5)
xlabel('Normalised frequency: \Omega = \omega/\omega_0','FontName',fname,'fontsize',fsize)
ylabel('Normalised respone: |U(\Omega)|/U_{st}','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');
h1=legend(['Nominal system'],['Responses at', sprintf('\n') , 'different slow times']);
set(h1,'FontName',fname,'FontSize',fsize,'Box','off','Location','northwest');
%print -depsc figs/freq_response_cloud_stiffness.eps
%print -djpeg figs/freq_response_cloud_stiffness.jpeg


%% Digital Twin derivation

j1=min(find( stiffness_array < 0.78))
delta_k1=stiffness_array(j1)

% original system - find the peak in the FRF
[H_0,id_max0]=max(X_k_om(:,1))   
f0=omega(id_max0)

% modified system: - find the peaks
[H_k,id_max_k]=max(X_k_om(:,j1))   
fk=omega(id_max_k)

R=H_0/H_k

fk^2-f0^2

% use the formula
Identified_stiffness=(1-zeta_0^2)*( R^2 -1 )
Original_stiffness=delta_k1-1
error_p=100*(Original_stiffness-Identified_stiffness)/Original_stiffness



% Identify at selected points in slow time with "exact" data
Identified_stiffness_array=zeros(nsamp+1,1);

for j1=2:nsamp+1   
    [H_k,id_max_k]=max(X_k_om(:,j1));   
    fk=omega(id_max_k);
    R=H_0/H_k;
    Identified_stiffness_array(j1)=(1-zeta_0^2)*( R^2 -1 );
    
end

% Plot original and DT
figure(6);clf;
plot(t1,StiffnessDegradeFunction,'-k','linewidth',1.5)
hold;
axis([0,1000,0,1.0]);
xlabel('Normalised slow time: t_s/T_0','FontName',fname,'fontsize',fsize)
ylabel('Normalised changes','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');
plot([0,TimeSampling],1+Identified_stiffness_array,'x:b','linewidth',1.5)
h1=legend('Actual system','Digital twin')
grid on
set(h1,'FontName',fname,'fontsize',fsize,'box','off','location','best')
axis([0,1000,0.5,1.0]);

nHF=10;

TimeSampling = [0,TimeSampling];
lt = length(TimeSampling);
tsel = round(linspace(10,lt-10,nHF),0);
tHF = TimeSampling(tsel);
tLF = TimeSampling;
KHF = Identified_stiffness_array(tsel);
KLF = 0.75.*Identified_stiffness_array(:) + 0.01.*sin(tLF(end) + tLF(:).*pi/10 .* Identified_stiffness_array(:));
% KLF = 0.55.*sin( tLF(:).*pi/10 .* Identified_stiffness_array(:));

close all

plot(TimeSampling,1+Identified_stiffness_array,'-r','Linewidth',2)
hold all
plot(tLF,1+KLF,'-g','Linewidth',2)

model_LF = fitrgp(tLF(:), 1+KLF(:), 'basis', 'pureQuadratic','KernelFunction','ardmatern32',...
    'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement','MaxObjectiveEvaluations',100),'Optimizer','lbfgs'); 

Ktilde_HF = predict(model_LF, tHF(:));

model_HF = fitrgp([tHF(:), Ktilde_HF(:)], 1+KHF(:), 'basis', 'pureQuadratic','KernelFunction','ardmatern32',...
    'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement','MaxObjectiveEvaluations',100),'Optimizer','lbfgs'); 


model_HF2 = fitrgp([tHF(:)], 1+KHF(:), 'basis', 'pureQuadratic','KernelFunction','ardmatern32',...
    'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement','MaxObjectiveEvaluations',100),'Optimizer','lbfgs'); 

tpred = linspace(0,1000,1001);
Kpred_LF = predict(model_LF,tpred(:));
Kpred_MF = predict(model_HF,[tpred(:),Kpred_LF(:)]);
Kpred_HF = predict(model_HF2,[tpred(:)]);

figure(10);
plot(tpred,Kpred_MF,'--r','Linewidth',2)
hold all;
plot(tpred,Kpred_HF,'-.g','Linewidth',2)
hold all;
plot(TimeSampling,1+Identified_stiffness_array,'-b','Linewidth',2);
title(['freq stiffness results for nhf=',num2str(nHF)])

%Indentified stiffness array is delta_k so KLF and KHF are also delta_k





