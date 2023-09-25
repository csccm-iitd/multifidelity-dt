clear all;clf;close all;

% fonts for the plot
fname='times';
fsize=14;
%% Normalised system
m0=1;
k0=1;
zeta0=0.02;
omega0=sqrt(k0/m0)
c0=2*zeta0*omega0


%% Define the mass chgane function
t1=linspace(0,1000,1000);
alpha2=0.15
eps2=0.35
MassFunction=1+eps2*sawtooth(alpha2*t1).*(sin(2.0*alpha2*t1)).^2;

% Plot changes
figure(1); hold all
plot(t1,MassFunction,'-k','linewidth',1.5)
xlabel('Normalised slow time: t_s/T_0','FontName',fname,'fontsize',fsize)
ylabel('Normalised mass \Delta_m(t_s)','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');


%% Changes in the natural frequency 
nsamp=500;

% Assume data is avialable at a lower smapling rate
TimeSampling=linspace(3,1000,nsamp);

K_t=k0;
M_t=m0*MassFunction;  % mass constant
omega_t=sqrt(K_t./M_t);
omega_i = interp1(t1,omega_t,TimeSampling); 

figure(2);clf
plot(t1,omega_t,'-k','linewidth',1.5);hold
plot(TimeSampling,omega_i,'*r');
axis([0,1000,0,1.5]);
grid on
xlabel('Normalised slow time: t_s/T_0','FontName',fname,'fontsize',fsize)
ylabel('Normalised change in the natural frequency ','FontName',fname,'fontsize',fsize);
h1=legend('Actual change','Samples available for digitial twin');
set(h1,'FontName',fname,'fontsize',fsize,'box','off','location','best')
ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');
%print -depsc figs/natural_frequency_changes_mass.eps
%print -djpeg figs/natural_frequency_changes_mass.jpeg


%% Time-domain response
mass_array=[1,interp1(t1,MassFunction,TimeSampling)];
tau_max=20;
n_tau=5001;
tau=linspace(0,tau_max,n_tau);
u0=1;

for j=1:nsamp+1

    omega_m=omega0/sqrt(mass_array(j));     % mass loaded natural frequency
    zeta_m=zeta0/sqrt(mass_array(j));           % mass loaded natural frequency
    omega_d_m=omega_m*sqrt(1-zeta_m^2);
    A_m=u0/sqrt(1-zeta_m^2);
    phi=atan(sqrt(1-zeta_m^2)/(zeta_m));

    f1=2*pi*zeta0*omega0/(mass_array(j));
    f2=2*pi*sqrt(mass_array(j)-zeta0^2)/mass_array(j);    
    %X_m_tau(:,j)=A_m*exp(-(zeta_n_m/sqrt(mass_array(j)))*2*pi*tau).*sin( (sqrt(1-zeta_n_m^2)/sqrt(mass_array(j)))*2*pi*tau+phi_m);
    X_m_tau(:,j)=A_m*exp(-f1*tau).*sin(f2*tau+phi);
    X_decay(:,j)=A_m*exp(-f1*tau);
    eigenvalue_all(:,j)=-zeta_m*omega_m + 1i * omega_d_m;

end

% Plot the response cloud
figure(4);clf;hold all
plot(tau,X_m_tau(:,1),'-k','linewidth',1.5)
plot(tau,X_m_tau(:,2:5:nsamp+1),'-y','linewidth',0.5)
plot(tau,X_m_tau(:,1),'-k','linewidth',1.5)
xlabel('Normalised time: \tau =(\omega_0/2\pi)t','FontName',fname,'fontsize',fsize)
ylabel('Normalised respone: u(\tau)/u_{0}','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');
h1=legend('Nominal system','Responses at different slow times');
set(h1,'FontName',fname,'FontSize',fsize,'Box','off','Location','southeast');
%print -depsc figs/time_response_cloud_mass.eps
%print -djpeg figs/time_response_cloud_mass.jpeg


j1=min(find( mass_array > 1.25))
delta_m1=mass_array(j1)

% Plot the identification method
figure(5);clf;hold all;

plot(tau,X_m_tau(:,1),'-k','linewidth',1.5)
plot(tau,X_m_tau(:,j1),'--r','linewidth',1.5)
plot(tau,X_decay(:,1),'-k','linewidth',1.5)
plot(tau,-X_decay(:,1),'-k','linewidth',1.5)
plot(tau,X_decay(:,j1),'-r','linewidth',1.5)
plot(tau,-X_decay(:,j1),'-r','linewidth',1.5)

    
% Plot x-axis
x_p=[0.0,tau_max];
y_p=[0,0];
plot(x_p,y_p,'-b','linewidth',1.5)

% Find two sucessive peaks and plot them
ut=X_m_tau(:,j1);
%plot(tau,ut,'-k')

[a1,b1]=find(tau==1)
[a2,b2]=max(ut(b1+1:n_tau))
x_p1=[tau(b1+b2),tau(b1+b2)];
y_p1=[0,a2];
plot(x_p1,y_p1,'-r','linewidth',1.5)

[a1,b1]=find(tau==2.2)
[a2,b2]=max(ut(b1+1:n_tau))
x_p1=[tau(b1+b2),tau(b1+b2)];
y_p1=[0,a2];
plot(x_p1,y_p1,'-r','linewidth',2.5)

% Annotations
text(5,0.65,'$e^{-\zeta_0 \omega_0 t}$','FontName',fname,'fontsize',24,'color','k','Interpreter','latex')
text(10,0.5,'$e^{-{\frac {\zeta_0 \omega_0} {1+\Delta} } t}$','FontName',fname,'fontsize',24,'color','r','Interpreter','latex')

text(0.880507841672891,1.05513784461153,'sucessive peaks','FontName',fname,'fontsize',fsize,'color','r')
% Create arrow
annotation('arrow',[0.203703703703704 0.180555555555555],...
    [0.788828337874659 0.731607629427793]);

% Create arrow
annotation('arrow',[0.342592592592593 0.350694444444444],...
    [0.682561307901907 0.628065395095368]);

% Create arrow
annotation('arrow',[0.209490740740741 0.225694444444444],...
    [0.795640326975477 0.690735694822888]);

% Create arrow
annotation('arrow',[0.549768518518518 0.570601851851852],...
    [0.640326975476839 0.599455040871935]);

xlabel('Normalised time: \tau =(\omega_0/2\pi)t','FontName',fname,'fontsize',fsize)
ylabel('Normalised respone: u(\tau)/u_{0}','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');

h1=legend(['Nominal system'],['Digital twin with: \Delta_m(t_s)='  num2str(delta_m1-1)]);
set(h1,'FontName',fname,'FontSize',fsize,'Box','off','Location','southeast');

%print -depsc figs/timed_response_details_mass.eps
%print -djpeg figs/timed_response_details_mass.jpeg
    
%% Digital Twin derivation


% check for one case (j1 - as above);

% original system - find the peaks
id_1=find(tau > 0.5 & tau < 1.5);
u1_0=max(X_m_tau(id_1,1))   

id_2=find(tau > 1.5 & tau < 2.75);
u2_0=max(X_m_tau(id_2,1))

delta_0=log(u1_0/u2_0)


% modified system: - find the peaks
id_1=find(tau > 0.5 & tau < 1.5);
u1_m=max(X_m_tau(id_1,j1))

id_2=find(tau > 1.5 & tau < 2.75);
u2_m=max(X_m_tau(id_2,j1))

delta_m=log(u1_m/u2_m)

% use the formula
Identified_mass=(1-zeta0^2)*( (delta_0/delta_m)^2 -1 )
Original_mass=delta_m1-1
error_p=100*(Original_mass-Identified_mass)/Original_mass



% Identify at selected points in slow time with "exact" data
Identified_mass_array=zeros(nsamp+1,1);

figure(1);
%plot([0,TimeSampling],mass_array,'*--r','linewidth',1.5)

for j1=2:nsamp+1
    
    u1_m=max(X_m_tau(id_1,j1));
    u2_m=max(X_m_tau(id_2,j1));
    delta_m=log(u1_m/u2_m);
    Identified_mass_array(j1)=(1-zeta0^2)*( (delta_0/delta_m)^2 -1 )
    

end

plot([0,TimeSampling],1+Identified_mass_array,'x:b','linewidth',1.5)
h1=legend('Actual system','Digital twin')
grid on
set(h1,'FontName',fname,'fontsize',fsize,'box','off','location','best')
axis([0,1000,0.5,1.5]);

nHF =19;

TimeSampling = [0,TimeSampling];
lt = length(TimeSampling);
tsel = round(linspace(10,lt-10,nHF),0);
tHF = TimeSampling(tsel);
tLF = TimeSampling;
KHF = Identified_mass_array(tsel);
%KLF=0.25*sawtooth(0.15*(tLF(:)-pi/0.15))
 KLF =0.75.*Identified_mass_array(:)+0.01*cos(tLF(:)*pi/10)+0.025
% KLF = 0.55.*sin( tLF(:).*pi/10 .* Identified_stiffness_array(:));
%KLF = 0.5.*Identified_mass_array(:)+0.1.*cos(tLF(:).*pi/10)
%KLF =(0.75.*Identified_mass_array(:)./sawtooth(alpha2*tLF(:))).*sawtooth(0.12*tLF(:))
%KLF =0.75.*Identified_mass_array(:)+0.01*cos(tLF(:)*pi/10)+0.025

close all;

plot(TimeSampling,1+Identified_mass_array,'-r','Linewidth',2)
hold all
plot(tLF,1+KLF,'-g','Linewidth',2)


model_LF = fitrgp(tLF(:),1+KLF(:), 'basis', 'pureQuadratic','KernelFunction','ardmatern32',...
    'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement','MaxObjectiveEvaluations',100),'Optimizer','lbfgs'); 

Ktilde_HF = predict(model_LF, tHF(:));

model_HF = fitrgp([tHF(:), Ktilde_HF(:)],1+KHF(:), 'basis', 'pureQuadratic','KernelFunction','ardmatern32',...
    'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement','MaxObjectiveEvaluations',100),'Optimizer','lbfgs'); 


model_HF2 = fitrgp([tHF(:)],1+KHF(:), 'basis', 'pureQuadratic','KernelFunction','ardmatern32',...
    'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement','MaxObjectiveEvaluations',100),'Optimizer','lbfgs'); 

tpred = linspace(0,1000,1001);
Kpred_LF = predict(model_LF,tpred(:));
Kpred_MF = predict(model_HF,[tpred(:),Kpred_LF(:)]);
Kpred_HF = predict(model_HF2,[tpred(:)]);

figure(10);
plot(tpred,Kpred_MF,'--r','Linewidth',2)
hold all;
plot(TimeSampling,1+Identified_mass_array,'-b','Linewidth',2)
hold all;
plot(tpred,Kpred_HF,'-.g','Linewidth',2)

h1=legend('Dig','Actual system')
grid on
set(h1,'FontName',fname,'fontsize',fsize,'box','off','location','best')
title(['time mass results for nhf=',num2str(nHF)])
%Identified_stiffness_array and KLF ,KHF are delta_k


