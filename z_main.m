% Code for the paper "A generic size- and trait-based model of plankton
% communities" 2020
%
%This code simulates the runs for the CONSTANT ENVIRONEMNT only 
%(i.e. environmental forcing is constant over time)
% The only things that need to be defined before runnning is the run-time
% and the number of copepod populations etc (all between lines 33 to 54)
%if nothing is changed the code will run as it has been run to make the
%figures for the paper
%
% Things to take into account:
%   the model needs to have at least one population of active feeders and one
%   of passive feeders.
%   Most figures have been adapted to the paper, and if the number of state
%   variables is changed some figures might not adapt to it.
%
% Structure of the code:
    % 1- define number of state variables and the environmental forcing
    % 2- the ODE system is solved
    % 3- Some initial figures are created
    % 4- The diagnostics are obtained by re-running the ode function
    % 5- The rest of plots are created
%
% Any questions related to this code can be addressed to me:
% cam.serra90@gmail.com
%
% Camila Serra-Pompei 25/08/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars

% Number of days to run the model
nbdays=20000; %run-days

%Environemntal parameters
Diff=0.005; % Input-rate of nutrients. Recommended range [10^-3 0.5].
Nmax=140; % Maximum N concentration
light=100; % Light
T=15; % Temperature

%Parameters setup
nbrP=14; %number of protists size classes
nbrC_act=8; %number of copepod active feeders populations
nbrC_pass=3; %number of copepods passive feeders populations
nbrFP=3; %number of fecal pellets size classes

P_min=1e-7; %minimum size for protists
P_max=0.1; %maximum size for protists
Cact_min=0.2; %minimum size for active feeders
Cact_max=1000; %maximum size for active feeders
Cpass_min=0.2; %minimum size for passive feeders
Cpass_max=5; %maximum size for passive feeders
C_size_classes=8; %Number of size-classes within each copepod population
nbrC=(nbrC_act+nbrC_pass)*C_size_classes; %total number of copepod state variables

% call parameters function
param=z_function_parameters(nbrP, P_min, P_max, nbrC_act, Cact_min, Cact_max,...
    nbrC_pass, Cpass_min, Cpass_max, C_size_classes, nbrFP);

% call function with predator-prey preferences
[theta_P_P, theta_cop_P, theta_cop_cop, theta_cop_F, theta_cop]=...
    z_function_feeding_kernels(param);

% Initical conditions
N0 = 1;
P0 = zeros(1,nbrP);
C0 = zeros(1,nbrC);
F0 = zeros(1,nbrFP);

P0(:) = 5;
C0(:) =5;

%------------------------------------------------
% Other things needed to run the model - Do not modify

% Diagnostics function is not activated
diags=0;

global time_counter
time_counter=1;

% indexing for mortality of higher trophic levels        
for i=1:length(param.Wvec)
    idx_bcop{i}=find(param.Wvec>=param.Wvec(i)/10^0.5 & param.Wvec<=param.Wvec(i)*10^0.5);
end

%options for the ODE solver
options1 = odeset('Refine',1);
options2 = odeset(options1,'NonNegative',1:1+param.nbr_P+param.nbr_cops*param.nbr_stages+param.nbr_fp);

%------------------------------------------------

% Solve ODE
[t, Y] = ode23(@z_ode_copepod_model, 0:1:nbdays, [N0, P0, C0, F0], options2, ...
                param, theta_P_P, theta_cop_P, theta_cop_cop , theta_cop_F,...
                light,T,diags,Diff,Nmax,idx_bcop);

% Rename output (Y)
N = Y(:,1);
P = Y(:,2:1+param.nbr_P);
C = Y(:,1+param.nbr_P+1:1+param.nbr_P+param.nbr_Ctot);
F = Y(:,1+param.nbr_P+param.nbr_Ctot+1:end);

%% Some initial figures     

figure
subplot(2,3,1)
plot(t,N,'k')
hold on
plot(t,F,'k--')
legend('N','F')
xlabel('Day')
ylabel('\mugC L^{-1}')
title('N and F\newline[\mugC L^{-1}]')

subplot(2,3,4)
plot(t,P)
xlabel('Day')
ylabel('\mugC L^{-1}')
title('Protists\newline[\mugC L^{-1}]')

subplot(2,3,2)
plot(t,C(:,param.ind_a(1:nbrC_act)))
for i=1:nbrC_act
    legend_info{i}=(['m_a=' num2str(param.Wa(i))]);
end
legend(legend_info)
xlabel('Day')
ylabel('\mugC L^{-1}')
title('Adult active cops.\newline[\mugC L^{-1}]')

subplot(2,3,3)
semilogy(t,C(:,param.ind_a(1:nbrC_act))*1000./param.Wa(1:nbrC_act))
xlabel('Day')
ylabel('# m^{-3}')
title('Adult active cops.\newline[# m^{-3}]')
ylim([1e-4 1e5])

subplot(2,3,5)
plot(t,C(:,param.ind_a(nbrC_act+1:end)))
for i=1:nbrC_pass
    legend_info{i}=(['m_a=' num2str(param.Wa(nbrC_act+i))]);
end
legend(legend_info)
xlabel('Day')
ylabel('\mugC L^{-1}')
title('Adult passive cops.\newline[\mugC L^{-1}]')

subplot(2,3,6)
semilogy(t,C(:,param.ind_a(nbrC_act+1:end)).*1000./param.Wa(nbrC_act+1:end))
xlabel('Day')
ylabel('# m^{-3}')
title('Adult passive cops.\newline[# m^{-3}]')
ylim([1e-4 1e5])


%% Diagnostics -----------------------------------------------------------
diags=1;

% prepare matrices to be filles
Flvl_diag=zeros(size(C));
dc_diag=zeros(size(C));
nupos_diag=zeros(size(C));
pred_C_on_C_diag=zeros(size(C));
dd_mort_C_diag=zeros(size(C));
mortP_diag=zeros(size(P));
pred_P_diag=zeros(size(P));
pred_C_on_P_diag=zeros(size(P));
mu_diag=zeros(size(P));
fg_diag=zeros(size(P));
fc_diag=zeros(size(C));
detrit_flux_diag=zeros(size(F));
reproduction_diag=zeros(size(C(:,param.ind_a)));
fracP_diag=zeros(size(C));
fracC_diag=zeros(size(C));
fracF_diag=zeros(size(C));
FPP_diag=zeros(size(C));
pred_C_on_F_diag=zeros(size(F));

% initialise for-loop to obtain the diagnostics from the model
for i=1:length(t)
    
      y_init=[N(i), P(i,:), C(i,:), F(i,:)]';
      [dydt, Flvl, d_c, pred_C_on_C, dd_mort_C, nu_pos, mortP , pred_P, pred_C_on_P,mu,fg,fc,...
          detrit_flux,reproduction,fracP,fracC,fracF,FPP,pred_C_on_F] = z_ode_copepod_model(i, y_init,...
              param, theta_P_P, theta_cop_P, theta_cop_cop , theta_cop_F,...
               light,T,diags,Diff,Nmax,idx_bcop); 
           
      Flvl_diag(i,:)=Flvl;     
      dc_diag(i,:)=d_c;   
      nupos_diag(i,:)=nu_pos;   
      pred_C_on_C_diag(i,:)=pred_C_on_C;
      dd_mort_C_diag(i,:)=dd_mort_C;
      mortP_diag(i,:)=mortP;
      pred_P_diag(i,:)=pred_P;
      pred_C_on_P_diag(i,:)=pred_C_on_P;
      mu_diag(i,:)=mu;
      fg_diag(i,:)=fg;
      fc_diag(i,:)=fc;
      detrit_flux_diag(i,:)=detrit_flux;
      reproduction_diag(i,:)=reproduction;
      fracP_diag(i,:)=fracP; 
      fracC_diag(i,:)=fracC;
      fracF_diag(i,:)=fracF;
      FPP_diag(i,:)=FPP;
      pred_C_on_F_diag(i,:)=pred_C_on_F;
    
end

% do the means of the last quarter of the run of the model
Pmean=mean(P(ceil(t(end)*3/4):end,:),1);
Cmean=mean(C(ceil(t(end)*3/4):end,:),1);
deltaC=param.deltaC;
deltaC(end,:)=deltaC(end-1,:);

Flvl_mean=mean(Flvl_diag(ceil(t(end)*3/4):end,:),1);
Flvl_max=max(Flvl_diag(ceil(t(end)*3/4):end,:));
Flvl_min=min(Flvl_diag(ceil(t(end)*3/4):end,:));

nu_mean=mean(nupos_diag(ceil(t(end)*3/4):end,:),1);
dc_mean=mean(dc_diag(ceil(t(end)*3/4):end,:),1);
predCC_mean=mean(pred_C_on_C_diag(ceil(t(end)*3/4):end,:),1);
dd_mort_C_mean=mean(dd_mort_C_diag(ceil(t(end)*3/4):end,:),1);
mortP_mean=mean(mortP_diag(ceil(t(end)*3/4):end,:),1);
pred_P_mean=mean(pred_P_diag(ceil(t(end)*3/4):end,:),1);
pred_C_on_P_mean=mean(pred_C_on_P_diag(ceil(t(end)*3/4):end,:),1);
mu_mean=mean(mu_diag(ceil(t(end)*3/4):end,:),1);
fg_mean=mean(fg_diag(ceil(t(end)*3/4):end,:),1);
fc_mean=mean(fc_diag(ceil(t(end)*3/4):end,:),1);


%% Plots for biomass spectrum and some diagnostics

%things neede to make the plots
[Wsort_a, idx_sort_a]=sort(param.Wvec(param.ind_act));
[Wsort_p, idx_sort_p]=sort(param.Wvec(param.ind_act(end)+1:end));

% Plots configurations
%Colors
c_ca=[29 41 81]./255;
c_cp=[87 160 211]./255;
c_p= [236 197 68]./255;

lwid=1.5;%linewidth

%plot size
x0=0;
y0=0;
width=17;
height=15;
grey=[0.5, 0.5, 0.5];
yellow=[244, 185, 66]./255;

%start figure
figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');

%use function to make figure
ha = tight_subplot(3,1,[.01 .01],[.1 .05],[.15 0.25])

% Subplot 1 - Biomass spectra %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ha(1))
scatter(0,0)
hold on

% Analytical solutions first ---------------------------
Flvl_mean_act=Flvl_mean(param.ind_act);
Flvl_mean_pass=Flvl_mean(param.ind_pass);
f0=mean([Flvl_mean_act, Flvl_mean_pass]);
%     f0=0.3;
    beta=mean([param.beta_act, param.beta_pass]);
    sigma=mean([param.sigma_act,param.sigma_pass]);
    m=logspace(-7,3);
    q=-0.25;
    n=-0.25;
    v=mean([0.011, 0.0052]);
    h=mean([1.37, 0.4]).*2.^((T(1)-15)/10);

    fncs_analytical=z_function_analytical_solutions();
    Nc=fncs_analytical.community_spectrum(f0,beta,sigma,m,h,v,q,n);
    mu=fncs_analytical.predation_mortality(f0,beta,sigma,m,h,n);

    h11=plot(m,Nc.*m,'k-.','linewidth',1); %make plot analytical solution
%---------------------------

% Numerical solutions
Pspec=Pmean./param.delta_V; %Protists biomass spectrum
Pspec(Pspec<1e-20)=NaN; %remove anything below 10^-20
h12=plot(param.V,Pspec,'color',c_p,'linewidth',lwid);
% make copepod population spectrum
st1=1;
st2=C_size_classes;
for i=1:nbrC_act+nbrC_pass
    biom_spec=Cmean(st1:st2)./deltaC(:,i)';

    if nu_mean(param.ind_a(i))<=0 || any(biom_spec<1e-40) %do not plot if these conditions are met
        
    else
        if i<=nbrC_act
            h13=plot(param.W(:,i),Cmean(st1:st2)./deltaC(:,i)','color',c_ca,'linewidth',lwid);
        else
            h14=plot(param.W(:,i),Cmean(st1:st2)./deltaC(:,i)','color',c_cp,'linewidth',lwid);
        end
    end

st1=st2+1;
st2=st2+C_size_classes;
end
ylim([1e-5 1e10])
xlim([1e-7 1e3])
ylabel('Biomass spectrum\newline[mgC m^{-3} \mugC^{-1}]')
set(gca,'xscale','log')
set(gca,'yscale','log')
legend([h12,h13,h14,h11],{'Protists','Active copepods','Passive copepods','Analytical solution'},'Fontsize',9)
legend boxoff



% Subplot 2 - Feedin levels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ha(2))
scatter(0,0)
hold on
plot(Wsort_a,Flvl_mean(param.ind_act(idx_sort_a))','color',c_ca,'linewidth',lwid)
plot(Wsort_p,Flvl_mean(param.ind_pass(idx_sort_p))','color',c_cp,'linewidth',lwid)
plot(param.V,fg','color',c_p,'linewidth',lwid)
h22=plot(logspace(log10(min(param.Wvec)),log10(max(param.Wvec)),10),fc_mean(1).*ones(1,10),'k--');
h21=plot(0,0,'k','linewidth',1.5);
xlim([1e-7 1e3])
ylim([0 0.8])
set(gca,'xscale','log')
ylabel('Feeding level [-]')
legend([h21,h22],{'Feeding level','Crtitcal feeding level'},'Fontsize',9)
legend boxoff



% Subplot 3 - Mortalities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ha(3))
scatter(0,0)
hold on
h34=plot(m,mu,'k-.','linewidth',1);
plot(Wsort_a,predCC_mean(param.ind_act(idx_sort_a))','color',c_ca,'linewidth',lwid)
plot(Wsort_p,predCC_mean(param.ind_pass(idx_sort_p))','color',c_cp,'linewidth',lwid)
st1=1;
st2=C_size_classes;
for i=1:nbrC_act+nbrC_pass
    if i<=nbrC_act
        plot(param.W(:,i),dd_mort_C_mean(st1:st2)',':','color',c_ca,'linewidth',lwid)
    else
        plot(param.W(:,i),dd_mort_C_mean(st1:st2)',':','color',c_cp,'linewidth',lwid)
    end
    
 st1=st2+1;   
 st2=st2+C_size_classes;
end
plot(param.V,mortP_mean',':','color',c_p,'linewidth',lwid)
plot(param.V,pred_P_mean','--','color',c_p,'linewidth',lwid)
plot(param.V,pred_C_on_P_mean','color',c_p,'linewidth',lwid)
h31=plot(0,0,'k-','linewidth',lwid);
h32=plot(0,0,'k--','linewidth',lwid);
h33=plot(0,0,'k:','linewidth',lwid);
legend([h31,h32,h33,h34],{'Predation by copepods','Predation by protists',...
    'Background or HTL mortality','Analytical solution'})
legend boxoff
ylim([1e-3 10^0.5])
xlim([1e-7 1e3])
set(gca,'xscale','log')
set(gca,'yscale','log')
ylabel('Mortality rates [d^{-1}]')
xlabel('Body mass [\mugC]')
set(gcf,'color','w');

