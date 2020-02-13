clearvars

%create scenario
seasonal_switch=1; %0 if we do not want seasonality; 1 if we want seasonality
linear_mortality_switch=1; % 1=Linear mortality, 0--> remove linear mortality
dd_mortality_switch=1; % 1=density dependent mortality on largest size classes, 0--> remove dd mortality
ontogeny_switch=1;
sensit_switch=0;
sensit_param=0;
diagnostics_switch=0;

param=parameters_copepod_model(linear_mortality_switch,dd_mortality_switch,ontogeny_switch,sensit_switch,sensit_param);

                                                           
if sensit_switch==0

    disp('Normal run')
    
    if seasonal_switch==1  
        disp('Seasonal environment') 
    elseif seasonal_switch==0
        disp('Steady environment') 
    end

months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];

param=parameters_copepod_model(linear_mortality_switch,dd_mortality_switch,ontogeny_switch,sensit_switch,sensit_param);

yearsrun=20;
param.D=10^-2.2;%-1.84;%0.01;
param.No=140;%200;%50; 
param.lat=55;

% if (param.lat==55 && seasonal_switch==1)
%    param.No=80; 
% end
% 
% if (param.lat==60 && seasonal_switch==1)
%    param.d_c=param.d_c.*0.01; %0.5
%    param.No=140; 
% end

tic

[temp,Diff,zmld,Ls,dzdt]=func_seasonal_forcing(param, yearsrun, seasonal_switch);
% Diff=load('Diff_season.mat');
% Diff=Diff.alk3;
% Diff=Diff';

%fix
    param.WD=param.WF;
    param.nbr_D=param.nbr_fp;
    param.D_dw=param.F_dw;
    param.D_up=param.F_up;


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


    


    options1 = odeset('Refine',1);
    options2 = odeset(options1,'NonNegative',1:1+param.nbr_P+param.nbr_cops*param.nbr_stages+param.nbr_fp);
    [t,x]= ode23(@(t,x) function_model_NPZF(t,x,param,theta_P_P,theta_cop,theta_cop_P,theta_cop_cop,...
        theta_cop_D,Ls,zmld,dzdt,Diff,temp,seasonal_switch,diagnostics_switch),0:1:365*yearsrun, [No, Po, Co, Fo],options2); 
%     options2 = odeset(options1,'NonNegative',1:1+param.nbr_P+param.nbr_cops*param.nbr_stages+param.nbr_fp*2);
%     [t,x]= ode23(@(t,x) function_model_NPZFFout(t,x,param,theta_P_P,theta_cop,theta_cop_P,theta_cop_cop,...
%         theta_cop_D,Ls,zmld,dzdt,Diff,temp,seasonal_switch,diagnostics_switch),0:1:365*yearsrun, [No, Po, Co, Fo, Fo],options2); 
%     
    N=x(:,1);
    P=x(:,2:param.nbr_P+1);
    C=x(:,param.nbr_P+2:param.nbr_P+param.nbr_cops*param.nbr_stages+1);
    D=x(:,param.nbr_P+param.nbr_cops*param.nbr_stages+2:end);
%     D=x(:,param.nbr_P+param.nbr_cops*param.nbr_stages+2:param.nbr_P+param.nbr_cops*param.nbr_stages+2+param.nbr_fp);
%     F_out=x(:,param.nbr_P+param.nbr_cops*param.nbr_stages+2+param.nbr_fp+1:end);
        %mass conservation?



toc

plot_fncs= function_plots();
fsize=10;

plot_fncs.plot_full_run(t,N,P,C,param,fsize)

%% some diagnostics

diagnostics_switch=1;

Flvl_diag=zeros(length(t),length(param.Wvec));
total_mortality_diag=zeros(length(t),length(param.Wvec));
dCdt_diag=zeros(length(t),length(param.Wvec));
total_mortC_diag=zeros(length(t),length(param.Wvec));
pred_C_C_diag=zeros(length(t),length(param.Wvec));
nu_diag=zeros(length(t),length(param.Wvec));
dd_mortC_diag=zeros(length(t),length(param.Wvec));
starvation_diag=zeros(length(t),length(param.Wvec));
reproduction_diag=zeros(length(t),length(param.ind_a));
Pfrac_diag=zeros(length(t),length(param.Wvec));
Cfrac_diag=zeros(length(t),length(param.Wvec));
Dfrac_diag=zeros(length(t),length(param.Wvec));

fg_diag=zeros(length(t),length(param.V));
pred_P_P_diag=zeros(length(t),length(param.V));
mu_diag=zeros(length(t),length(param.V));
pred_C_P_diag=zeros(length(t),length(param.V));
dPdt_diag=zeros(length(t),length(param.V));
NPP_diag=zeros(length(t),length(param.V));
mback_diag=zeros(length(t),length(param.V));
JN_diag=zeros(length(t),length(param.V));
JL_diag=zeros(length(t),length(param.V));
JF_diag=zeros(length(t),length(param.V));
light_diag=zeros(1,length(t));
detrit_flux_diag=zeros(length(t),param.nbr_fp);
pred_C_on_D_diag=zeros(length(t),param.nbr_fp);
FPP_prod_diag=zeros(length(t),1);

for i=1:length(t)
    

        [dxdt, output_diagnostic, Diagnostics_vec_protists, reproduction,detrit_flux,light,pred_C_on_D,FPP_production]= function_model_NPZF(t(i),x(i,:)',param,theta_P_P,...
            theta_cop,theta_cop_P,theta_cop_cop, theta_cop_D,Ls,zmld,dzdt,Diff,temp,seasonal_switch,diagnostics_switch); 
        

   
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
    light_diag(i)=light;
    end
    light_diag(i)=light;%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
if seasonal_switch>0
figure
plot(t(end-365:end),light_diag(end-365:end),'k')
end

figure
plot(t,light_diag,'k')
if abs(light_diag(end))>1e-5
    disp('no mass balance!')
    
end
%%
plot_fncs= function_plots();
fsize=10;

if seasonal_switch==0

    plot_fncs.plot_numberspec(t,N,P,C,param,fsize,fg_diag,Flvl_diag,total_mortality_diag,pred_C_C_diag,dd_mortC_diag,...
                                    pred_P_P_diag,pred_C_P_diag,ontogeny_switch,seasonal_switch,mback_diag)

    plot_fncs.plot_diagnostics_protists(param,fsize,JL_diag,JN_diag,JF_diag,mu_diag,pred_P_P_diag,...
                                    pred_C_P_diag,mback_diag)
       
         mean_mort=mean(total_mortC_diag(ceil(t(end)*3/4):end,:));
         mean_nu=mean(nu_diag(ceil(t(end)*3/4):end,:));
         meanflux=mean(detrit_flux_diag(ceil(t(end)*3/4):end,:));
         flux1000=sum(meanflux.*exp(-param.remin./param.sink_fp.*900),2);
     figure; semilogx(param.Wvec,mean_nu,'ko'); hold on; semilogx(param.Wvec,mean_mort,'ro')
     figure; semilogx(flux1000,'ko')
                                
elseif seasonal_switch>0
yr=9;%30
yrend=8*365;
    plot_fncs.plot_annual(t,N,P,C,param,fsize,yr,yrend)
    plot_fncs.plot_annual_biomass(t,N,P,C,param,fsize,yr,yrend)
    if param.lat==55
        plot_fncs.plot_lat55(t,N,P,C,param,yearsrun,fsize)
    elseif param.lat==60
        plot_fncs.plot_lat60(t,N,P,C,param,fsize)
    end
%     plot_fncs.plot_protists(t,N,P,C,param,fsize)
%     plot_fncs.plot_chla(t,N,P,C,param,fsize)
    plot_fncs.plot_cohorts(t,N,P,C,param,fsize,yr,yrend)
%     
%     plot_fncs.plot_diagnostics(t,N,P,C,param,fsize,fg_diag,Flvl_diag,pred_C_P_diag,pred_C_C_diag,pred_P_P_diag,dd_mortC_diag)
%     plot_fncs.plot_growth_rates(t,N,P,C,param,fsize,mu_diag,dPdt_diag,nu_diag,dCdt_diag)
%     yr=9;
%     yrend=(yr-1)*365;
    plot_fncs.plot_mortalities(t,N,P,C,param,fsize,starvation_diag,pred_C_C_diag, dd_mortC_diag,yr,yrend)
%     plot_fncs.plot_nu_cohorts(t,N,P,C,param,fsize,nu_diag)
yr=18;% 18
yrend=10*365; %10
    plot_fncs.plot_reproduction(t,N,P,C,param,fsize,reproduction_diag,yr,yrend,mu_diag)
    plot_fncs.plot_all_stacked(t,N,P,C,param,fsize,yr,yrend,zmld)
    plot_fncs.plot_predations(t,N,P,C,param,fsize,pred_C_P_diag,pred_C_C_diag,pred_P_P_diag,dd_mortC_diag,yr,yrend)
    plot_fncs.plot_all(t,N,P,C,param,fsize,yr,yrend)
    plot_fncs.plot_Flvlsources(param,fsize,t,Pfrac_diag,Cfrac_diag,Dfrac_diag,Flvl_diag,yr,yrend)
    
     figure; scatter(Diff(end-365:end),light_diag(end-365:end),10,1:366,'filled'); colormap(viridis(12))
     colorbar
     set(gca,'xscale','log')

end
                            
elseif sensit_switch==1
    disp('Sensitivity analysis mode')
    
    sensit_param=logspace(-2.3,-2.3,1);
% %     sensit_param=linspace(1,15,40);
%     sensit_param=linspace(0,30,30);
    
%     [N,P,C,D,Pcell,Ccell,Nmean,Pmean,Cmean,Dmean]=function_call_sensitivity_analysis(seasonal_switch,linear_mortality_switch,...
%         dd_mortality_switch,ontogeny_switch,sensit_switch,sensit_param,diagnostics_switch,detritus_switch);
% 
%       [Ncell,Pcell,Ccell,Dcell]=function_call_sensitivity_analysis_parfor(seasonal_switch,linear_mortality_switch,...
%         dd_mortality_switch,ontogeny_switch,sensit_switch,sensit_param,diagnostics_switch,detritus_switch);
    
      [Ncell,Pcell,Ccell,Dcell,NPPcell,fluxcell]=function_call_sensitivity_analysis_diags(seasonal_switch,linear_mortality_switch,...
        dd_mortality_switch,ontogeny_switch,sensit_switch,sensit_param,diagnostics_switch,detritus_switch);
    
    
    
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
    

% filename = 'ws_cte_2.mat';
% save(filename)
end