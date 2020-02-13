clearvars

%create scenario
%if singe rune is required then:
seasonal_switch=0; %set to 0 for the steady environment, set to 1 for the seasoanl scenario
linear_mortality_switch=1; % 1=Linear mortality, 0--> remove linear mortality
dd_mortality_switch=1; % 1=density dependent mortality on largest size classes, 0--> remove dd mortality
%fro both constant and seasonal environment:
yearsrun=10; %for how many years to we want the run
N_input_rate=10^-1.3; %parameter of input of nutrients to the system (only works for constan environment, seasonal is defined in "func_seasonal_forcing")
No_deep=140; % Nitrogen concentration of the deep layer
latitude=55; %latitude (only for the seasonal scenario)

%if parameter sweep is required (this can take quite long, unless you run in parallel, can be set inside "function_call_sensitivity_analysis_diags")
sensit_switch=1; %set to 0 for a simple run, set to 1 to run the parameter sweep (seasonal_switch must be =0)
if sensit_switch==1
   seasonal_switch=0; 
end
step_nb=60; %number of iterations
min_val=-3; %min values of parameter sweep. note that here we do in log, so the min values here is in fcat 10^-3
max_val=-2; %max values of parameter sweep.
param_sweep=logspace(min_val,max_val,step_nb); %if not in log scale, change to "linspace"
parallerun_mode=0;%set to 1 if you want to run in parallel (parfol toolbox needed)
%---------
ontogeny_switch=1; %do not change
diagnostics_switch=0; %do not swith
sensit_param=0; %do not change
%----------


%Load parameters
% param=parameters_copepod_model(linear_mortality_switch,dd_mortality_switch,ontogeny_switch,sensit_switch,sensit_param);
months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];

                                                           
if sensit_switch==0 % perform simple run

    disp('Normal run')
    
    if seasonal_switch==1  
        disp('Seasonal environment') 
    elseif seasonal_switch==0
        disp('Steady environment') 
    end

param=parameters_copepod_model(linear_mortality_switch,dd_mortality_switch,ontogeny_switch,sensit_switch,sensit_param);
param.D=N_input_rate; %parameter of input of nutrients to the system (only for constan environment)
param.No=No_deep; % Nitrogen concentration of the deep layer
param.lat=latitude; %latitude (only for the seasonal scenario)

tic

% Load environmental forcings
[temp,Diff,zmld,Ls,dzdt]=func_seasonal_forcing(param, yearsrun, seasonal_switch);

%load detritus sizes (it is a correction from notations of old versions)
param.WD=param.WF;
param.nbr_D=param.nbr_fp;
param.D_dw=param.F_dw;
param.D_up=param.F_up;

%load feeding kernels
[theta_P_P,theta_P_D, theta_P, theta_cop_P, theta_cop_cop, theta_cop_D, theta_cop]=func_feeding_kernels(param);


%Initial conditions
No=5;
Po=zeros(1,param.nbr_P);
Po(:)=1;
Co=zeros(1,length(param.Wvec));
Co(:)=1;
Do=zeros(1,param.nbr_D);
Do(:)=0;
Fo=zeros(1,param.nbr_fp);
Fo(:)=0;

%call ODE solver
    options1 = odeset('Refine',1);
    options2 = odeset(options1,'NonNegative',1:1+param.nbr_P+param.nbr_cops*param.nbr_stages+param.nbr_fp);
    [t,x]= ode23(@(t,x) function_model_NPZF(t,x,param,theta_P_P,theta_cop,theta_cop_P,theta_cop_cop,...
        theta_cop_D,Ls,zmld,dzdt,Diff,temp,seasonal_switch,diagnostics_switch),0:1:365*yearsrun, [No, Po, Co, Fo],options2); 

%Rename output
    N=x(:,1);
    P=x(:,2:param.nbr_P+1);
    C=x(:,param.nbr_P+2:param.nbr_P+param.nbr_cops*param.nbr_stages+1);
    D=x(:,param.nbr_P+param.nbr_cops*param.nbr_stages+2:end);


toc

% call function for plots 
plot_fncs= function_plots();
fsize=10;
plot_fncs.plot_full_run(t,N,P,C,param,fsize)

%% some diagnostics

diagnostics_switch=1;

%create empty matrices to fill in the loop
Flvl_diag=zeros(length(t),length(param.Wvec)); %feeding level of copepods
total_mortality_diag=zeros(length(t),length(param.Wvec)); %total mortality of copepods
dCdt_diag=zeros(length(t),length(param.Wvec)); %copepods ode
total_mortC_diag=zeros(length(t),length(param.Wvec)); %total mortality of copepods
pred_C_C_diag=zeros(length(t),length(param.Wvec)); %predation by copepods on copepods
nu_diag=zeros(length(t),length(param.Wvec)); %net energy gain copepods
dd_mortC_diag=zeros(length(t),length(param.Wvec)); %mortality by HTL on copepods
starvation_diag=zeros(length(t),length(param.Wvec)); %starvation ortality copepods
reproduction_diag=zeros(length(t),length(param.ind_a)); %reproduction copepods
Pfrac_diag=zeros(length(t),length(param.Wvec)); %fraction of the feeding level of cops. that comes from eating protists
Cfrac_diag=zeros(length(t),length(param.Wvec)); %fraction of the feeding level of cops. that comes from eating copepods
Dfrac_diag=zeros(length(t),length(param.Wvec)); %fraction of the feeding level of cops. that comes from eating fecal pellets

fg_diag=zeros(length(t),length(param.V)); %feeding level protists
pred_P_P_diag=zeros(length(t),length(param.V)); %predation by protists on protists
mu_diag=zeros(length(t),length(param.V)); %net growth protists
pred_C_P_diag=zeros(length(t),length(param.V)); %predation by copepods on protists
dPdt_diag=zeros(length(t),length(param.V)); %protists ode
NPP_diag=zeros(length(t),length(param.V)); %Net primary production
mback_diag=zeros(length(t),length(param.V)); %background mortality protists
JN_diag=zeros(length(t),length(param.V)); %uptake of nitrogen by protists
JL_diag=zeros(length(t),length(param.V)); %photosynthetic rate
JF_diag=zeros(length(t),length(param.V)); %food uptake by protists

massbalance_diag=zeros(1,length(t)); %light experienced by protists (takes into account MLD and shading)
detrit_flux_diag=zeros(length(t),param.nbr_fp); %flux of detritus out of the ML
pred_C_on_D_diag=zeros(length(t),param.nbr_fp); %consumption of fecal pellets by copepods
FPP_prod_diag=zeros(length(t),1); %fecal pellets production by copepods

%start diagnostics calculations
for i=1:length(t)
    
    %call function with systems of ode
        [dxdt, output_diagnostic, Diagnostics_vec_protists, reproduction,detrit_flux,massbalance,pred_C_on_D,FPP_production]= function_model_NPZF(t(i),x(i,:)',param,theta_P_P,...
            theta_cop,theta_cop_P,theta_cop_cop, theta_cop_D,Ls,zmld,dzdt,Diff,temp,seasonal_switch,diagnostics_switch); 

   %fill matrices at each time-step
    Flvl_diag(i,:)=output_diagnostic(:,1);
    total_mortality_diag(i,:)=output_diagnostic(:,2);
    dCdt_diag(i,:)=output_diagnostic(:,3);
    total_mortC_diag(i,:)=output_diagnostic(:,4);
    pred_C_C_diag(i,:)=output_diagnostic(:,5);
    nu_diag(i,:)=output_diagnostic(:,6);
    dd_mortC_diag(i,:)=output_diagnostic(:,7);
    starvation_diag(i,:)=output_diagnostic(:,8);
    if seasonal_switch==1
    Pfrac_diag(i,:)=output_diagnostic(:,9);
    Cfrac_diag(i,:)=output_diagnostic(:,10);
    Dfrac_diag(i,:)=output_diagnostic(:,11);
    massbalance_diag(i)=massbalance;
    end
    massbalance_diag(i)=massbalance;
    fg_diag(i,:)=Diagnostics_vec_protists(:,1);
    pred_P_P_diag(i,:)=Diagnostics_vec_protists(:,2);
    mu_diag(i,:)=Diagnostics_vec_protists(:,3);
    pred_C_P_diag(i,:)=Diagnostics_vec_protists(:,4);
    dPdt_diag(i,:)=Diagnostics_vec_protists(:,5);
    NPP_diag(i,:)=Diagnostics_vec_protists(:,6);
    mback_diag(i,:)=Diagnostics_vec_protists(:,7);
    JN_diag(i,:)=Diagnostics_vec_protists(:,8);
    JL_diag(i,:)=Diagnostics_vec_protists(:,9);
    JF_diag(i,:)=Diagnostics_vec_protists(:,10);
    pred_C_on_D_diag(i,:)=pred_C_on_D;
    FPP_prod_diag(i)=FPP_production;
    
    detrit_flux_diag(i,:)=detrit_flux';
    reproduction_diag(i,:)=reproduction;
    
end

total_mortality_diag=total_mortality_diag-starvation_diag;

%check mass balance
if seasonal_switch==0
    figure
    plot(t,massbalance_diag,'k')
        if abs(massbalance_diag(end))>1e-5
            disp('no mass balance!')   
        end
end
%% PLOTS
plot_fncs= function_plots();
fsize=10;

if seasonal_switch==0

    plot_fncs.plot_numberspec(t,N,P,C,param,fsize,fg_diag,Flvl_diag,total_mortality_diag,pred_C_C_diag,dd_mortC_diag,...
                                    pred_P_P_diag,pred_C_P_diag,ontogeny_switch,seasonal_switch,mback_diag)

    plot_fncs.plot_diagnostics_protists(param,fsize,JL_diag,JN_diag,JF_diag,mu_diag,pred_P_P_diag,...
                                    pred_C_P_diag,mback_diag)

                                
elseif seasonal_switch>0
yr=1;
yrend=0*365;

    plot_fncs.plot_annual(t,N,P,C,param,fsize,yr,yrend)
    plot_fncs.plot_annual_biomass(t,N,P,C,param,fsize,yr,yrend)
    plot_fncs.plot_cohorts(t,N,P,C,param,fsize,yr,yrend)
    plot_fncs.plot_mortalities(t,N,P,C,param,fsize,starvation_diag,pred_C_C_diag, dd_mortC_diag,yr,yrend)
    plot_fncs.plot_reproduction(t,N,P,C,param,fsize,reproduction_diag,yr,yrend,mu_diag)
    plot_fncs.plot_all_stacked(t,N,P,C,param,fsize,yr,yrend,zmld)
    plot_fncs.plot_predations(t,N,P,C,param,fsize,pred_C_P_diag,pred_C_C_diag,pred_P_P_diag,dd_mortC_diag,yr,yrend)
    plot_fncs.plot_all(t,N,P,C,param,fsize,yr,yrend)
    plot_fncs.plot_Flvlsources(param,fsize,t,Pfrac_diag,Cfrac_diag,Dfrac_diag,Flvl_diag,yr,yrend)

end
                            
elseif sensit_switch==1
    disp('Sensitivity analysis mode')
    
    sensit_param=param_sweep;

      [Ncell,Pcell,Ccell,Dcell,NPPcell,fluxcell]=function_call_sensitivity_analysis_diags(seasonal_switch,linear_mortality_switch,...
        dd_mortality_switch,ontogeny_switch,sensit_switch,sensit_param,diagnostics_switch,parallerun_mode);
    
    
    
%%    

for i=1:length(sensit_param)
Nmean(i)=mean(Ncell{i});
Pmean(i,:)=mean(Pcell{i});
Cmean(i,:)=mean(Ccell{i});
Dmean(i,:)=mean(Dcell{i});
end


for i=1:length(sensit_param)
C=Ccell{i};
Cmin(i,:)=min(C);
Cmax(i,:)=max(C);

end

%%

sensit_param_log=(sensit_param);
param=parameters_copepod_model(linear_mortality_switch,dd_mortality_switch,ontogeny_switch,sensit_switch,sensit_param);

fillcolor='k';

figure
subplot(4,1,1)
plot(sensit_param_log,Nmean,'ko','markersize',3,'markerfacecolor','k')
ylabel('[\mugN L^{-1}]')
set(gca,'xscale','log')

subplot(4,1,2)
cmap=flip(inferno(param.nbr_P+1));
h1=plot(sensit_param_log,Pmean,'linewidth',1.5);
for i=1:length(cmap)-1
  set(h1(i) ,'Color', cmap(1+i,:));
end
ylim([1e-2 1e9])
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([min(sensit_param) max(sensit_param)])
ylabel('[# L^{-1}]')

subplot(4,1,3)
cmap=flip(inferno(param.C_sp/2+1));
h1=plot(sensit_param_log,Cmean(:,param.ind_a(1:param.C_sp_act)),'o-','markersize',5); %./param.Wvec(param.ind_a(1:param.C_sp_act))'*1000
hold on
for i=1:length(cmap)-1
  set(h1(i) ,'markeredgecolor', cmap(1+i,:));
  set(h1(i) ,'markerfacecolor', cmap(1+i,:));
  legendInfo{i} = ['m_a= ' num2str(param.Wvec(param.ind_a(i))) '\mugC']; 
end
% funtion_plots_shadings(sensit_param_log,Cmin(:,param.ind_a(1:param.C_sp_act))./param.Wvec(param.ind_a(1:param.C_sp_act))'*1000,...
%     Cmax(:,param.ind_a(1:param.C_sp_act))./param.Wvec(param.ind_a(1:param.C_sp_act))'*1000,fillcolor) 
funtion_plots_shadings(sensit_param_log,Cmin(:,param.ind_a(1:param.C_sp_act)),...
    Cmax(:,param.ind_a(1:param.C_sp_act)),fillcolor)
ylabel('[# m^{-3}]')
% legend(legendInfo)
% legend('Location','northwest')
set(gca,'xscale','log')
set(gca,'yscale','log')
ylim([1e-3 1e7])
xlim([min(sensit_param) max(sensit_param)])

subplot(4,1,4)
cmap=flip(inferno(param.C_sp/2+1));
h2=plot(sensit_param_log,Cmean(:,param.ind_a(param.C_sp_act+1:end)),'o-','markersize',5); %./param.Wvec(param.ind_a(param.C_sp_act+1:end))'*1000
for i=1:length(cmap)-1
  set(h2(i) ,'markeredgecolor', cmap(1+i,:));
  set(h2(i) ,'markerfacecolor', cmap(1+i,:));
%   set(h2(i),'MarkerFaceAlpha',0.5)
  legendInfo{i} = ['m_a= ' num2str(param.Wvec(param.ind_a(i))) '\mugC']; 
end
hold on
fillcolor='k';
funtion_plots_shadings(sensit_param_log,Cmin(:,param.ind_a(param.C_sp_act+1:end)),...
    Cmax(:,param.ind_a(param.C_sp_act+1:end)),fillcolor)   
legend(legendInfo)
legend('Location','northwest')
set(gca,'Yscale','log')
set(gca,'xscale','log')
xlabel('Diffusivity term [d^{-1}]')
% xlabel('Irradiance term [W m^{-2}]')
ylabel('[# m^{-3}]')
ylim([1e-3 1e7])
xlim([min(sensit_param) max(sensit_param)])
fsize=10;
plot_settings(fsize)
    

% filename = 'ws_cte_3.mat';
% save(filename)
end