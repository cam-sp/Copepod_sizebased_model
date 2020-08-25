% Code for the paper "A generic size- and trait-based model of plankton
% communities" 2020 
%
%This code simulates the runs for the SEASONAL ENVIRONEMNT only 
%(i.e. environmental forcing varies with time)
% The only things that need to be defined before runnning is the run-time
% and the number of copepod populations etc (all between lines 33 to 47)
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

tic
yearsrun=50; % years of the run

%Parameters setup
nbrP=14; %number of protists size classes
nbrC_act=8; %number of copepod active feeders populations
nbrC_pass=3; %number of copepods passive feeders populations
nbrFP=10; %number of fecal pellets size classes

P_min=1e-7; %minimum size for protists
P_max=0.1; %maximum size for protists
Cact_min=0.2; %minimum size for active feeders
Cact_max=1000; %maximum size for active feeders
Cpass_min=0.2; %minimum size for passive feeders
Cpass_max=5; %maximum size for passive feeders
C_size_classes=8; %Number of size-classes within each copepod population
nbrC=(nbrC_act+nbrC_pass)*C_size_classes; %total number of copepod state variables

%call parameters function
param=z_function_parameters(nbrP, P_min, P_max, nbrC_act, Cact_min, Cact_max,...
    nbrC_pass, Cpass_min, Cpass_max, C_size_classes, nbrFP);

% call function with predator-prey preferences
[theta_P_P, theta_cop_P, theta_cop_cop, theta_cop_F, theta_cop]=...
    z_function_feeding_kernels(param);

%call function with the seasonal forcing
[T,Diff,zmld,light,dzdt]=z_func_seasonal_forcing(param, yearsrun);
Nmax=140; %deep N concentration

% Initical conditions
N0 = 1;
P0 = zeros(1,nbrP);
P0(:) = 5;
C0 = zeros(1,nbrC);
C0(:) =5;
F0 = zeros(1,nbrFP);

%------------------------------------------------
% Other things needed to run the model - Do not modify

% Diagnostics function is not activated
diags=0;
seasonal_switch=1;

global time_counter
time_counter=1;

% indexing for mortality of higher trophic levels      
for i=1:length(param.Wvec)
    idx_bcop{i}=find(param.Wvec>=param.Wvec(i)/10^0.5 & param.Wvec<=param.Wvec(i)*10^0.5);
end
        
tic

%options for the ODE solver
options1 = odeset('Refine',1);
options2 = odeset(options1,'NonNegative',1:1+param.nbr_P+param.nbr_cops*param.nbr_stages+param.nbr_fp);

%------------------------------------------------

% Solve ODE
[t, Y] = ode23(@z_ode_copepod_model, 1:1:365*yearsrun, [N0, P0, C0, F0], options2, ...
                param, theta_P_P, theta_cop_P, theta_cop_cop , theta_cop_F,...
                light,T,diags,Diff,Nmax,idx_bcop,seasonal_switch,zmld,dzdt);

toc 

% Rename output (Y)
N = Y(:,1);
P = Y(:,2:1+param.nbr_P);
C = Y(:,1+param.nbr_P+1:1+param.nbr_P+param.nbr_Ctot);
F = Y(:,1+param.nbr_P+param.nbr_Ctot+1:end);



figure
subplot(2,3,1)
plot(t,N,'k')

subplot(2,3,2)
plot(t,P)

subplot(2,3,3)
plot(t,C(:,param.ind_a(1:nbrC_act)))
legend

subplot(2,3,4)
plot(t,C(:,param.ind_a(nbrC_act+1:end)))
legend

subplot(2,3,5)
plot(t,F)
title('F')

figure
subplot(2,3,1)
plot(t(end-365:end),N(end-365:end),'k')

subplot(2,3,2)
plot(t(end-365:end),P(end-365:end,:))

subplot(2,3,3)
plot(t(end-365:end),C(end-365:end,param.ind_a(1:nbrC_act)))
legend

subplot(2,3,4)
plot(t(end-365:end),C(end-365:end,param.ind_a(nbrC_act+1:end)))
legend

subplot(2,3,5)
plot(t(end-365:end),F(end-365:end,:))
title('F')



%% Diagnostics -----------------------------------------------------------
diags=1;

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
reproduction_diag=zeros(length(t),length(param.Wa));
fracP_diag=zeros(size(C));
fracC_diag=zeros(size(C));
fracF_diag=zeros(size(C));
FPP_diag=zeros(size(C));
pred_C_on_F_diag=zeros(size(F));
FPremin_diag=zeros(size(F));

for i=1:length(t)
    
      y_init=[N(i), P(i,:), C(i,:), F(i,:)]';
      [dydt, Flvl, d_c, pred_C_on_C, dd_mort_C, nu_pos, mortP , pred_P, pred_C_on_P,mu,fg,fc,...
          detrit_flux,reproduction,fracP,fracC,fracF,FPP,pred_C_on_F,Fpremin] = z_ode_copepod_model(i, y_init,...
              param, theta_P_P, theta_cop_P, theta_cop_cop , theta_cop_F,...
               light,T,diags,Diff,Nmax,idx_bcop,seasonal_switch,zmld,dzdt); 
           
                      
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
      pred_C_on_F_diag(i,:)=pred_C_on_F.*F(i,:);
      Fpremin_diag=Fpremin;
end

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

Pmean=mean(P(ceil(t(end)*3/4):end,:),1);
Cmean=mean(C(ceil(t(end)*3/4):end,:),1);
deltaC=param.deltaC;
% deltaC(end,:)=deltaC(end-1,:);


%% Main plot seasonal environment

yr=1;
yrend=0*365; 

fluxmld=sum(detrit_flux_diag,2);
flux1000=sum(detrit_flux_diag.*exp(-param.remin./param.sink_fp.*(1000-zmld)),2);

fsize=10;
rownbs=1;

s1=1e-6;
s2=1e-5;
s3=1e-4;
s4=1e-3;
s5=1e-2;

tinterval=t(end-365*yr):30.5:t(end)-yrend;
Gyear=P(end-365*yr:end-yrend,:);
gg1=sum(Gyear(:,param.V<s1),2);
gg2=sum(Gyear(:,param.V>=s1 & param.V<s2),2);
gg3=sum(Gyear(:,param.V>=s2 & param.V<s3),2);
gg4=sum(Gyear(:,param.V>=s3 & param.V<s4),2);
gg5=sum(Gyear(:,param.V>=s4 & param.V<s5),2);
gg6=sum(Gyear(:,param.V>=s5),2);

gg=cat(2,gg1,gg2,gg3,gg4,gg5,gg6);%.*zmld(end-365*yr:end-yrend);%./(14*5.6);
months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
bluemap=(brewermap(6,'YlGnBu'));


x0=0;
y0=0;
width=17;
height=16;
grey=[0.5, 0.5, 0.5];
yellow=[244, 185, 66]./255;

cp= [0.5882    0.5882    0.5882];

fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');

ha = tight_subplot(4,2,[.01 .2],[.1 .02],[.11 .15]);

%subplot Protists
axes(ha(1))
yyaxis left
h =bar(t(end-365*yr:end-yrend),gg,1,'stacked','EdgeColor','none');
plots=get(gca, 'Children');
hold on
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
ylabel({'Protists'; '[mg m^{-3}]'})
set(h,{'FaceColor'},{bluemap(1,:);bluemap(2,:);bluemap(3,:);bluemap(4,:);bluemap(5,:);bluemap(6,:)});
xlabel('Month')
legend boxoff
set(gca,'Xticklabel',[])
yyaxis right
plot(t(end-365*yr:end-yrend),N(end-365*yr:end-yrend),':','color',cp,'linewidth',1.5)
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cp;
ylabel('Nitrogen [\mugN L^{-1}]')

legf=legend(plots([1 2 3 4 5 6]),flip({'           m < 10^{-6}','10^{-6} \leq m < 10^{-5}','10^{-5} \leq m < 10^{-4}','10^{-4} \leq m < 10^{-3}','10^{-3} \leq m < 10^{-2}','10^{-2} \leq m'}));
title(legf,{'Legend panel a';'Protiosts mass ranges'})
legf.FontSize = 8;
 
%size classes to make body-mass ranges of copepods
sc1=1e-2;
sc2=1e-1;
sc3=1e0;
sc4=1e1;
sc5=1e2;

%group copepods in the body-mass ranges
Cyear=C(end-365*yr:end-yrend,param.ind_act);
cc1_a=sum(Cyear(:,param.W(param.ind_act)<sc1),2);
cc2_a=sum(Cyear(:,param.W(param.ind_act)>=sc1 & param.W(param.ind_act)<sc2),2);
cc3_a=sum(Cyear(:,param.W(param.ind_act)>=sc2 & param.W(param.ind_act)<sc3),2);
cc4_a=sum(Cyear(:,param.W(param.ind_act)>=sc3 & param.W(param.ind_act)<sc4),2);
cc5_a=sum(Cyear(:,param.W(param.ind_act)>=sc4 & param.W(param.ind_act)<sc5),2);
cc6_a=sum(Cyear(:,param.W(param.ind_act)>=sc5),2);

Cyear_p=C(end-365*yr:end-yrend,param.ind_pass);
cc1_p=sum(Cyear_p(:,param.W(param.ind_pass)<sc1),2);
cc2_p=sum(Cyear_p(:,param.W(param.ind_pass)>=sc1 & param.W(param.ind_pass)<sc2),2);
cc3_p=sum(Cyear_p(:,param.W(param.ind_pass)>=sc2 & param.W(param.ind_pass)<sc3),2);
cc4_p=sum(Cyear_p(:,param.W(param.ind_pass)>=sc3 & param.W(param.ind_pass)<sc4),2);
cc5_p=sum(Cyear_p(:,param.W(param.ind_pass)>=sc4 & param.W(param.ind_pass)<sc5),2);
cc6_p=sum(Cyear_p(:,param.W(param.ind_pass)>=sc5),2);

cc_a=cat(2,cc1_a,cc2_a,cc3_a,cc4_a,cc5_a,cc6_a);
cc_p=cat(2,cc1_p,cc2_p,cc3_p,cc4_p,cc5_p,cc6_p);

%subplot active feeders biomass
axes(ha(3))
h3=bar(t(end-365*yr:end-yrend),cc_a,1,'stacked');
hold on
set(h3,{'FaceColor'},{bluemap(1,:);bluemap(2,:);bluemap(3,:);bluemap(4,:);bluemap(5,:);bluemap(6,:)});
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
for i=13:12:length(tinterval)-1
plot([tinterval(i), tinterval(i)], [1e-5 1e8],'k--')
end
set(gca,'XtickLabels',months)
ylabel({'Active copepods';'[mg m^{-3}]'})
plot_settings(fsize)
set(gca,'Xticklabel',[])


%subplot passive feeders biomass
axes(ha(5))
h5=bar(t(end-365*yr:end-yrend),cc_p,1,'stacked');
hold on
set(h5,{'FaceColor'},{bluemap(1,:);bluemap(2,:);bluemap(3,:);bluemap(4,:);bluemap(5,:);bluemap(6,:)});
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
for i=13:12:length(tinterval)-1
plot([tinterval(i), tinterval(i)], [1e-5 1e8],'k--')
end
set(gca,'XtickLabels',months)
ylabel({'Passive copepods';'[mg m^{-3}]'})
plots=get(gca, 'Children');
legf=legend(plots([1 2 3 4 5 6]),flip({'m < 10^{-2} ','10^{-2} \leq m < 10^{-1}','10^{-1} \leq m < 10^{0}','10^{0} \leq m < 10^{1}','10^{1} \leq m < 10^{2}','10^{2} \leq m'}));
legend boxoff
legf.FontSize = 8;
title(legf,{'Legend panels b and d';'Copepods mass ranges'})
set(gca,'Xticklabel',[])
plot_settings(fsize)

% subplot fecal pellets fluxes
axes(ha(7))
yyaxis left
plot(t(end-365*yr:end-yrend),fluxmld(end-365*yr:end-yrend),'k','linewidth',1.5) %./zmld(end-365*yr:end-yrend)
hold on
plot(t(end-365*yr:end-yrend),flux1000(end-365*yr:end-yrend),'k--','linewidth',1.5)
ylabel({'Fecal pellets flux';'[mg m^{-2} d^{-1}]'})
set(gca,'YTick',[0 20 40 60]);
yyaxis right
plot(t(end-365*yr:end-yrend),flux1000(end-365*yr:end-yrend)./fluxmld(end-365*yr:end-yrend),'color',cp,'linewidth',1.5)
ylabel({'Transfer efficiency [-]'})
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
plot_settings(fsize)
xlim([t(end-365*yr) t(end)-yrend])
legend({'MLD', '1000m', '1000m:MLD'},'FontSize',8)
legend boxoff
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cp;
plot_settings(fsize)
xlabel('Month')

deltaC=param.deltaC;
deltaC(end,:)=param.deltaC(param.ind_a-1);

%subplot biomass spectrum cohort
axes(ha(4))
i=6; %which cohort do we want. this is the index for the adult of the population we want to show
st1=(i-1)*param.nbr_stages+1;
st2=st1-1+param.nbr_stages;
contourf(t(end-365*yr:end-yrend),param.Wvec(st1:st2),log10(C(end-365*yr:end-yrend,st1:st2)'./(deltaC(:,i))),20,'edgecolor','none')
set(gca,'yscale','log')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
bmap=(brewermap(10,'YlGnBu'));
colormap((bmap))
yticks([1e0 1e1 1e2])
ylim([1e0 1e2])
cbh=colorbar;
cbh.Label.String = 'Biomass spectrum [\mugC \mugC^{-1} L^{-1}]';
ylabel({'Body mass';'[\mugC]'})
plot_settings(fsize)
set(gca,'Xticklabel',[])

% subplot for the reproduction rate of the population
axes(ha(6))
yyaxis right
plot(t(end-365*yr:end-yrend),reproduction_diag(end-365*yr:end-yrend,i),'color',cp,'linewidth',1.5)
ylim([0 0.05])
yticks([0 0.01 0.02 0.03 0.04])
ylabel({'Specific Rep. rate ';'[d^{-1}]'})
yyaxis left
plot(t(end-365*yr:end-yrend),reproduction_diag(end-365*yr:end-yrend,i).*C(end-365*yr:end-yrend,param.ind_a(i)),'k','linewidth',1.5)
xlim([t(end-365*yr) t(end)-yrend])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cp;
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
ylabel({'Rep. rate ';'[mgC m^{-3} d^{-1}]'})
plot_settings(fsize)
set(gca,'Xticklabel',[])

%we need to sort to make the contourfs (we could make it as a surface, 
%but when compiling the figure for the paper it looks weird)
[Wsort, sortind]=sort(param.Wvec(param.ind_act));

pred_C_C_diag=pred_C_on_C_diag;

%subplot for predation rates
axes(ha(8))
bmap=flip(viridis(20));
contourf(t(end-365*yr:end-yrend),Wsort,(pred_C_C_diag(end-365*yr:end-yrend,sortind)'),50,'edgecolor','none')
cbh=colorbar;
cbh.Label.String = 'Predation [d^{-1}]';
set(gca,'yscale','log')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
colormap(flip(bmap))
axis tight
ylabel({'Body mass';'[\mugC]'})
plot_settings(fsize)
xlabel('Month')
yticks([1e-2 1e0 1e2])


%% reproduction rates
%Note that this plot has been made for the paper only, and if the number of
%population changes it might give an error
figure

for i=1:nbrC_act+nbrC_pass
    plotsuborder=[flip([1 3 5 7 9 11 13 15]),[16 14 12]];
    subplot(8,2,plotsuborder(i))
yyaxis right
semilogy(t(end-365*yr:end-yrend),reproduction_diag(end-365*yr:end-yrend,i),'o','color',cp,'markerfacecolor',cp,'markersize',1)
ylim([1e-4 1e0])
if i==1 || i==9
    ylabel({'Specific Rep. rate ';'[d^{-1}]'})
end
yyaxis left
plot(t(end-365*yr:end-yrend),reproduction_diag(end-365*yr:end-yrend,i).*C(end-365*yr:end-yrend,param.ind_a(i)),'ko','markerfacecolor','k','markersize',2)
ylim([1e-4 0.1])
set(gca,'yscale','log')
xlim([t(end-365*yr) t(end)-yrend])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cp;
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
if i==1 || i==9
    ylabel({'Rep. rate ';'[mgC m^{-3} d^{-1}]'})
end
plot_settings(fsize)
set(gca,'Xticklabel',[])
alk=sprintf('m_a= %9.1E',param.Wa(i));
text(t(end-365*yr+10),0.075,alk)

if plotsuborder(i)==15 || plotsuborder(i)==16
    xlabel('Month')
    set(gca,'XtickLabels',months)
end
end


%%
cc_large=sum(cc_a(:,5:6),2)+sum(cc_p(:,5:6),2);


%% Plot feeding levels

[Wsort_act, ind_s_act]=sort(param.Wvec(param.ind_act));
Flvl_diag_act=Flvl_diag(:,param.ind_act);

[Wsort_pass, ind_s_pass]=sort(param.Wvec(param.ind_pass));
Flvl_diag_pass=Flvl_diag(:,param.ind_pass);

figure
subplot(1,2,1)
contourf(t(end-365*yr:end-yrend),Wsort_act,Flvl_diag_act(end-365*yr:end-yrend,ind_s_act)',10,'edgecolor','none')
set(gca,'yscale','log')
colorbar
caxis([0 0.8])
title('Flvl active feeders')
ylim([min(param.Wvec) max(param.Wvec)])
xlim([t(end-365*yr) t(end-yrend)])
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
ylabel('Body mass\newline    [\mugC]')
xlabel('Month')


subplot(1,2,2)
contourf(t(end-365*yr:end-yrend),Wsort_pass,Flvl_diag_pass(end-365*yr:end-yrend,ind_s_pass)',10,'edgecolor','none')
set(gca,'yscale','log')
colorbar
caxis([0 0.8])
colormap(viridis)
title('Flvl passive feeders')
ylim([min(param.Wvec) max(param.Wvec)])
xlim([t(end-365*yr) t(end-yrend)])
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
ylabel('Body mass\newline    [\mugC]')
xlabel('Month')

set(gcf,'color','w');

%% Plot predation rates


[Wsort_act, ind_s_act]=sort(param.Wvec(param.ind_act));
pred_C_on_C_diag_act=pred_C_on_C_diag(:,param.ind_act);

[Wsort_pass, ind_s_pass]=sort(param.Wvec(param.ind_pass));
pred_C_on_C_diag_pass=pred_C_on_C_diag(:,param.ind_pass);



figure
subplot(4,1,1)
contourf(t(end-365*yr:end-yrend),param.V,pred_P_diag(end-365*yr:end-yrend,:)',10,'edgecolor','none')
set(gca,'yscale','log')
colorbar
title('predation by protists on protists')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
ylabel('Cell mass [\mugC]')
xlabel('Month')
set(gca,'YTick',[1e-6 1e-4 1e-2]);
colormap(viridis(10))
axis tight


subplot(4,1,2)
contourf(t(end-365*yr:end-yrend),param.V,pred_C_on_P_diag(end-365*yr:end-yrend,:)',10,'edgecolor','none')
set(gca,'yscale','log')
colorbar
title('Predation by copepods on protists')
axis tight
set(gca,'YTick',[1e-6 1e-4 1e-2]);
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
ylabel('Body mass [\mugC]')
xlabel('Month')
colormap(viridis(10))


subplot(4,1,3)
contourf(t(end-365*yr:end-yrend),Wsort_act,pred_C_on_C_diag_act(end-365*yr:end-yrend,ind_s_act)',10,'edgecolor','none')
set(gca,'yscale','log')
colorbar
caxis([0 0.12])
colormap(viridis)
set(gca,'YTick',[1e-2 1e0 1e2]);
title('Predation by copepods on active feeders')
ylim([min(param.Wvec) max(param.Wvec)])
xlim([t(end-365*yr) t(end-yrend)])
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
ylabel('Body mass [\mugC]')
xlabel('Month')


subplot(4,1,4)
contourf(t(end-365*yr:end-yrend),Wsort_pass,pred_C_on_C_diag_pass(end-365*yr:end-yrend,ind_s_pass)',10,'edgecolor','none')
set(gca,'yscale','log')
colorbar
title('Predation by copepods on passive feeders')
ylim([min(param.Wvec) max(param.Wvec)])
xlim([t(end-365*yr) t(end-yrend)])
caxis([0 0.12])
set(gca,'YTick',[1e-2 1e0 1e2]);
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
ylabel('Body mass [\mugC]')
xlabel('Month')
colormap(viridis(10))

set(gcf,'color','w');




%% Plot fecal pellets fluxes

D=F;
FPPtot=sum(FPP_diag,2).*zmld;


figure
subplot(3,1,1)
plot(t(end-365*yr:end-yrend),flux1000(end-365*yr:end-yrend)./fluxmld(end-365*yr:end-yrend),'b','linewidth',1.5)
hold on
plot(t(end-365*yr:end-yrend),fluxmld(end-365*yr:end-yrend)./FPPtot(end-365*yr:end-yrend),'k','linewidth',1.5)
plot(t(end-365*yr:end-yrend),sum(pred_C_on_F_diag(end-365*yr:end-yrend,:),2).*zmld(end-365*yr:end-yrend)...
    ./FPPtot(end-365*yr:end-yrend),'r','linewidth',1.5)
legend('Flux_{1000}/Flux_{MLD}','Flux_{MLD}/FPP','\mu_{F}/FPP')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
xlim([t(end-365*yr) t(end)-yrend])
legend boxoff
ylabel('Ratio [-]')

subplot(3,1,2)
plot(t(end-365*yr:end-yrend),-zmld(end-365*yr:end-yrend),'k','linewidth',1.5)
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
xlim([t(end-365*yr) t(end)-yrend])
ylabel('MLD [m]')

subplot(3,1,3)
plot(t(end-365*yr:end-yrend),cc_large./(sum(cc_a,2)+sum(cc_p,2)),'k','linewidth',1.5)
title('Fraction of copepods \newlinelarger than 10\mugC')
ylabel('[-]')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
xlim([t(end-365*yr) t(end)-yrend])
xlabel('Month')

set(gcf,'color','w');
