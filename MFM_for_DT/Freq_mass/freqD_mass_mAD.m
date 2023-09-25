clear all; close all; clc; clf
fname='times';
fsize=20;

omega_0=1;
u_0=1;

n_values=4;
Q_array=[5,10,20,50]
%zeta_array=1./sqrt(1+4.*Q_array.^2)

zeta_array=[0.1, 0.05, 0.025, 0.01]

nf=3001;
f_max=2;
omega=linspace(0,f_max,nf);


% with additional mass
nsamp=3;
%Delta=linspace(0,0.5,nsamp)

Delta_array=[0,0.1,0.5];


mass_array=1+Delta_array;

color_array='krb'
figure(1);clf;hold all;

linS = {'-','--','-.',':'};

eigenvalue_all=zeros(1,nsamp,n_values);

for k=1:n_values
zeta_0=zeta_array(k);

    for j=1:nsamp

        omega_m=omega_0/sqrt(mass_array(j));     % mass loaded natural frequency
        zeta_m=zeta_0/sqrt(mass_array(j));           % mass loaded natural frequency
        omega_d_m=omega_m*sqrt(1-zeta_m^2);
        
        X_m_om(:,j,k)=1./sqrt( (1-mass_array(j)*omega.^2).^2 + 4*zeta_0^2*omega.^2);
        X2_m_om2(:,j,k)=1./sqrt( (1-mass_array(j)*omega).^2 + 4*zeta_0^2*omega);

        Delta=mass_array(j)-1;
        Hmax_p(j,k)= sqrt(1+Delta-2*zeta_0^2)/(1+Delta);
        H_max(j,k)=(1/2*(1+Delta))*sqrt(1/(zeta_0^2*(Delta+1-zeta_0^2)));
        
        Half_amp(j,k)=H_max(j,k)/sqrt(2);
        B=mass_array(j);p=zeta_0;       
        Half_power_points_high(j,k)=sqrt(B-2*p^2+2*sqrt(B*p^2-p^4))/B;
        Half_power_points_low(j,k)=sqrt(B-2*p^2-2*sqrt(B*p^2-p^4))/B;
        Half_power_BandWidth(j,k)=Half_power_points_high(j,k)- Half_power_points_low(j,k);
        
        eigenvalue_all(:,j,k)=-zeta_m*omega_m + i * omega_d_m;

    end

end


omega2=linspace(0,f_max^2,nf^2);


%% Freq response plot
fsize=26;
for k=1:n_values
    figure(k);clf;hold all;
    for j=1:nsamp
        semilogy(omega,(X_m_om(:,j,k)),'LineStyle',linS{j},'linewidth',1.5,'color',color_array(j));
    end
    if k==1
        h1=legend(['Nominal systems'],...
          ['\Delta_m='  num2str(Delta_array(2))],['\Delta_m='  num2str(Delta_array(3))]);
        set(h1,'FontName',fname,'FontSize',fsize,'Box','off','Location','northeast')
    end

    xlabel('Normalised frequency: \Omega = \omega/\omega_0','FontName',fname,'fontsize',fsize)
    ylabel('Normalised respone: |U|/U_{st}','FontName',fname,'fontsize',fsize);
    %axis([0 , 2, 0, 4]);
    ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','off');

    %eval(['print -depsc figs/freqd_response_Q' num2str(Q_array(k)) '.eps']);
    %eval(['print -djpeg figs/freqd_response_Q' num2str(Q_array(k)) '.jpeg']);
        
end

%% Maxima of the freqnecy response polt and analysis
color_array='kbr'

linS = {'-','-.','--',':'};


fsize=20;
k=3;
figure(5);clf;hold all;
for j=1:2:3    
    plot(omega,X_m_om(:,j,k),'LineStyle',linS{j},'linewidth',1.5,'color',color_array(j));
end
h1=legend(['Nominal',sprintf('\n'),'system'],['\Delta_m (t_s)='  num2str(Delta_array(3))]);
set(h1,'FontName',fname,'FontSize',fsize,'Box','off','Location','northwest');
for j=1:2:3    
    x_p=[ Hmax_p(j,k), Hmax_p(j,k)]
    y_p=[0,H_max(j,k)]
    plot(x_p,y_p,'LineStyle',':','linewidth',1.5,'color',color_array(j),'Marker','*','MarkerSize',10,'HandleVisibility','off')
end

% Texts
text(1.1,10,'frequency shift','FontName',fname,'fontsize',fsize)
text(1.1,22,'peak-response shift','FontName',fname,'fontsize',fsize)

text(1.05,20,'${\mathcal H}_0$','FontName',fname,'fontsize',20,'color','k','Interpreter','latex','HandleVisibility','off')
text(0.64,24,'${\mathcal H}_m$','FontName',fname,'fontsize',20,'color','r','Interpreter','latex','HandleVisibility','off')

% Annotations
% Create doublearrow
annotation('doublearrow',[0.444695259593679 0.516930022573363],...
    [0.318 0.318]);

% Create arrow
annotation('arrow',[0.586907449209932 0.489841986455982],...
    [0.402547671840355 0.321507760532151]);

% Create line
annotation('line',[0.4785 0.5056],...
    [0.7575 0.7575]);

% Create line
annotation('line',[0.4785 0.5056],...
    [0.9147 0.9147]);

% Create doublearrow
annotation('doublearrow',[0.4943 0.4943],...
    [0.9147 0.7575]);

% Create arrow
annotation('arrow',[0.551918735891648 0.498871331828442],...
    [0.821616407982262 0.83370288248337]);

xlabel('Normalised frequency: \Omega = \omega/\omega_0','FontName',fname,'fontsize',fsize)
ylabel('Normalised respone: |U|/U_{st}','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');

%print -depsc figs/freqd_response_details.eps
%print -djpeg figs/freqd_response_details.jpeg


%%
clear all; close all; clc; clf
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


% Assume data is avialable at a lower smapling rate
nsamp=500;
TimeSampling=linspace(3,1000,nsamp);



%% Frequency-domain response
mass_array=[1,interp1(t1,MassFunction,TimeSampling)];
nf=5001;
f_max=2;
omega=linspace(0,f_max,nf);
zeta_0=0.02

for j=1:nsamp+1
    
    X_m_om(:,j)=1./sqrt( (1-mass_array(j)*omega.^2).^2 + 4*zeta_0^2*omega.^2);

end

% Plot the response cloud
figure(4);clf;hold all
plot(omega,X_m_om(:,1),'-k','linewidth',1.5)
plot(omega,X_m_om(:,2:3:nsamp+1),'-y','linewidth',0.5)
plot(omega,X_m_om(:,1),'-k','linewidth',1.5)
xlabel('Normalised frequency: \Omega = \omega/\omega_0','FontName',fname,'fontsize',fsize)
ylabel('Normalised respone: |U(\Omega)|/U_{st}','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');
h1=legend(['Nominal system'],['Responses at', sprintf('\n') , 'different slow times']);
set(h1,'FontName',fname,'FontSize',fsize,'Box','off','Location','northwest');
%print -depsc figs/freq_response_cloud_mass.eps
%print -djpeg figs/freq_response_cloud_mass.jpeg


%% Digital Twin derivation

j1=min(find( mass_array > 1.25))
delta_m1=mass_array(j1)

% original system - find the peak in the FRF
[H_0,id_max0]=max(X_m_om(:,1))   
f0=omega(id_max0)

% modified system: - find the peaks
[H_m,id_max_m]=max(X_m_om(:,j1))   
fm=omega(id_max_m)

R=H_m/H_0;

% use the formula
Identified_mass=(1+zeta0^2)*( R^2 -1 )
Original_mass=delta_m1-1
error_p=100*(Original_mass-Identified_mass)/Original_mass


% Identify at selected points in slow time with "exact" data
Identified_mass_array=zeros(nsamp+1,1);


for j1=2:nsamp+1
    
    [H_m,id_max_m]=max(X_m_om(:,j1))   
    fm=omega(id_max_m)
    
    Identified_mass_array(j1)=(1+zeta0^2)*( (H_m/H_0)^2 -1 )
    

end



% Plot original and DT
figure(6); hold all
plot(t1,MassFunction,'-k','linewidth',1.5)
xlabel('Normalised slow time: t_s/T_0','FontName',fname,'fontsize',fsize)
ylabel('Normalised mass \Delta_m(t_s)','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');
plot([0,TimeSampling],1+Identified_mass_array,'+:r','linewidth',1.5)
h1=legend('Actual system','Digital twin')
grid on
set(h1,'FontName',fname,'fontsize',fsize,'box','off','location','best')
axis([0,1000,0.5,1.5]);

nHF=12;

TimeSampling = [0,TimeSampling];
lt = length(TimeSampling);
tsel = round(linspace(10,lt-10,nHF),0);
tHF = TimeSampling(tsel);
tLF = TimeSampling;
KHF = Identified_mass_array(tsel);
%%KLF=0.24*sawtooth(0.15*(tLF(:)-pi/0.15))
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
title(['freq mass results for nhf =',num2str(nHF)])
%Identified_stiffness_array and KLF ,KHF are delta_k



