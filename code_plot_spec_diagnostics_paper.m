%% code to plot figure 5 in the paper ()
% biomass spectrum + diagnostic of constant environment
clearvars

biom_spec=0;

for ji=1:2
    
    scen=ji;

    %load runs
if scen==1
    load('ws_cte_1.mat')
elseif scen==2
    load('ws_cte_2.mat')
end
%%
ba=[29 41 81]./255;
bp=[87 160 211]./255;
cp= [236 197 68]./255;

[Wsort_act, idxWsort_act]=sort(param.Wvec(param.ind_act));
[Wsort_pass, idxWsort_pass]=sort(param.Wvec(param.ind_pass));
                                
deltaC=param.deltaC;
deltaC(end,:)=param.deltaC(param.ind_a-1);

sheldon=1;

if sheldon==1
Pmean=mean(P(ceil(t(end)*3/4):end,:),1).*1000./param.V;
Pspec=Pmean.*param.V.^2./(param.delta_V);%sheldon spectrum
Cmean=mean(C(ceil(t(end)*3/4):end,:),1).*1000./param.Wvec';
Cspec=Cmean.*param.Wvec'.^2./(deltaC(:)');%sheldon spectrum

Pmax=max(P(ceil(t(end)*3/4):end,:)).*1000.*param.V./(param.delta_V);
Cmax=max(C(ceil(t(end)*3/4):end,:)).*1000.*param.Wvec'./(deltaC(:)');

Pmin=min(P(ceil(t(end)*3/4):end,:)).*1000.*param.V./(param.delta_V);
Cmin=min(C(ceil(t(end)*3/4):end,:)).*1000.*param.Wvec'./(deltaC(:)');
Cmin(Cmin==0)=1e-15;
Pmin(Pmin==0)=1e-15;
Pmax(Pmax==0)=1e-15;

else
    
Pmean=mean(P(ceil(t(end)*3/4):end,:),1);
Pspec=Pmean./(param.delta_V);%biomass spectrum
Cmean=mean(C(ceil(t(end)*3/4):end,:),1);
Cspec=Cmean./(deltaC(:)');%biomass spectrum

Pmax=max(P(ceil(t(end)*3/4):end,:))./(param.delta_V);
Cmax=max(C(ceil(t(end)*3/4):end,:))./(deltaC(:)');

Pmin=min(P(ceil(t(end)*3/4):end,:))./(param.delta_V);
Cmin=min(C(ceil(t(end)*3/4):end,:))./(deltaC(:)');
Cmin(Cmin==0)=1e-15;
Pmin(Pmin==0)=1e-15;
Pmax(Pmax==0)=1e-15;


end

if biom_spec==1
Pmean=mean(P(ceil(t(end)*3/4):end,:),1);
Cmean=mean(C(ceil(t(end)*3/4):end,:),1);

Pmax=max(P(ceil(t(end)*3/4):end,:)).*1000./param.V;
Cmax=max(C(ceil(t(end)*3/4):end,:)).*1000./param.Wvec';

Pmin=min(P(ceil(t(end)*3/4):end,:)).*1000./param.V;
Cmin=min(C(ceil(t(end)*3/4):end,:)).*1000./param.Wvec';
end

fg_mean=mean(fg_diag(ceil(t(end)*3/4):end,:),1);
mean_nu=mean(nu_diag(ceil(t(end)*3/4):end,:),1);
nu_mean_act=mean(nu_diag(ceil(t(end)*3/4):end,param.ind_act),1);
nu_mean_pass=mean(nu_diag(ceil(t(end)*3/4):end,param.ind_pass),1);
Flvl_mean_act=mean(Flvl_diag(ceil(t(end)*3/4):end,param.ind_act),1);
Flvl_mean_pass=mean(Flvl_diag(ceil(t(end)*3/4):end,param.ind_pass),1);
total_mortality_mean_act=mean(total_mortality_diag(ceil(t(end)*3/4):end,param.ind_act),1);
total_mortality_mean_pass=mean(total_mortality_diag(ceil(t(end)*3/4):end,param.ind_pass),1);
pred_C_C_mean_act=mean(pred_C_C_diag(ceil(t(end)*3/4):end,param.ind_act),1);
pred_C_C_mean_pass=mean(pred_C_C_diag(ceil(t(end)*3/4):end,param.ind_pass),1);
ddmort_mean_act=mean(dd_mortC_diag(ceil(t(end)*3/4):end,param.ind_act),1);
ddmort_mean_pass=mean(dd_mortC_diag(ceil(t(end)*3/4):end,param.ind_pass),1);
pred_P_P_mean=mean(pred_P_P_diag(ceil(t(end)*3/4):end,:),1);
pred_C_P_mean=mean(pred_C_P_diag(ceil(t(end)*3/4):end,:),1);
mu_mean=mean(mu_diag(ceil(t(end)*3/4):end,:),1);
% mback_mean=mean(mback(ceil(t(end)*3/4):end,:),1);

nu_min_act=min(nu_diag(ceil(t(end)*3/4):end,param.ind_act));
nu_min_act(nu_min_act<=0)=1e-15;

nu_min_pass=min(nu_diag(ceil(t(end)*3/4):end,param.ind_pass));
nu_min_pass(nu_min_pass<=0)=1e-15;

nu_max_act=max(nu_diag(ceil(t(end)*3/4):end,param.ind_act));
nu_max_act(nu_max_act<=0)=1e-15;
nu_max_pass=max(nu_diag(ceil(t(end)*3/4):end,param.ind_pass));
nu_max_pass(nu_max_pass<=0)=1e-15;

tot_mort_prtoists=pred_P_P_diag+pred_C_P_diag+mback_diag;
tot_mort_prtoists_mean=mean(tot_mort_prtoists(ceil(t(end)*3/4):end,:),1);

Flvl_max_act=max(Flvl_diag(ceil(t(end)*3/4):end,param.ind_act));
Flvl_max_pass=max(Flvl_diag(ceil(t(end)*3/4):end,param.ind_pass));

Flvl_min_act=min(Flvl_diag(ceil(t(end)*3/4):end,param.ind_act));
Flvl_min_pass=min(Flvl_diag(ceil(t(end)*3/4):end,param.ind_pass));




lowxlim=1e-7;
highxlim=1e4;
%%


cp= [236 197 68]./255;


ba=[29 41 81]./255;
bp=[87 160 211]./255;
cp= [236 197 68]./255;

if scen==1
x0=0;
y0=0;
width=17;
height=15;
grey=[0.5, 0.5, 0.5];
yellow=[244, 185, 66]./255;

figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');

ha = tight_subplot(3,2,[.01 .01],[.1 .05],[.15 0.25])

axes(ha(1))
    
elseif scen==2
    axes(ha(2))
end


%%
Pspec(Pspec<1e-6)=NaN;
scatter(0,0)
set(gca,'yscale','log')
set(gca,'xscale','log')
hold on
p1=loglog(param.V,Pspec,'o-','color',cp,'linewidth', 1.5,'markersize',2)
funtion_plots_shadings(param.V,Pmin,Pmax,'k')
if ontogeny_switch==1
st1=1;
st2=param.nbr_stages;
for i=1:param.C_sp_act
p2=loglog(param.Wvec(st1:st2),Cspec(st1:st2),'color',ba,'linewidth', 1.5)
funtion_plots_shadings(param.Wvec(st1:st2),Cmin(st1:st2),Cmax(st1:st2),ba)

st1=st1+param.nbr_stages;
st2=st2+param.nbr_stages;
end
st1=param.nbr_stages.*param.C_sp_act+1;
st2=param.nbr_stages.*param.C_sp_act+param.nbr_stages;
for i=1:param.C_sp_act
p3=loglog(param.Wvec(st1:st2),Cspec(st1:st2),'color',bp,'linewidth', 1.5)
funtion_plots_shadings(param.Wvec(st1:st2),Cmin(st1:st2),Cmax(st1:st2),bp)
st1=st1+param.nbr_stages;
st2=st2+param.nbr_stages;
end


else
    loglog(param.Wvec(param.ind_a(1:param.C_sp_act)),Cmean(param.ind_a(1:param.C_sp_act))./deltaC(param.ind_a(1:param.C_sp_act)),'k','linewidth', 1.5)
    loglog(param.Wvec(param.ind_a(param.C_sp_act+1:end)),Cmean(param.ind_a(param.C_sp_act+1:end))./deltaC(param.ind_a(param.C_sp_act+1:end)),'r','linewidth', 1.5)
    
end

ylim([1e-5 1e7])
xlim([lowxlim highxlim])

if scen==1
ylabel({'Biomass spectrum';'[\mugC L^{-1} \mugC^{-1}]'})
end


%%%%%%%%%%%%%%%%%%%%%%%%
f0=mean([Flvl_mean_act, Flvl_mean_pass]);
%     f0=0.3;
    beta=mean([param.beta_act, param.beta_pass]);
    sigma=mean([param.sigma_act,param.sigma_pass]);
    m=logspace(-8,4);
    q=-0.25;
    n=-0.25;
    v=mean([0.011, 0.0052]);
    h=mean([1.37, 0.4]).*2.^((temp(1)-15)/10);

    fncs_analytical=function_analytical_solutions();
    Nc=fncs_analytical.community_spectrum(f0,beta,sigma,m,h,v,q,n);
    mu=fncs_analytical.predation_mortality(f0,beta,sigma,m,h,n);

%%%%%%%%%%%%%%%%%%%%%%%%%%

p4=loglog(m,Nc.*m,'k:','linewidth',1.5)
if scen==2
legend([p1(1),p2(1),p3(1),p4(1)],'Protists','Active cops.','Passive cops.','Analytical spectrum')
legend boxoff
end

if scen==1
axes(ha(3))
elseif scen==2
axes(ha(4))
end
scatter(0,0)
hold on
p1=plot(param.V, fg_mean,'color',cp,'linewidth', 1.5)
p2=plot(Wsort_act, Flvl_mean_act(idxWsort_act),'color',ba,'linewidth', 1.5)
p3=plot(Wsort_act, param.fc(param.ind_act),'--','color',ba,'linewidth', 1)
p4=plot(Wsort_pass, Flvl_mean_pass(idxWsort_pass),'color',bp,'linewidth', 1.5)
fc_pass=param.fc(param.ind_pass);
p5=plot(Wsort_pass, fc_pass(idxWsort_pass),'--','color',bp,'linewidth', 1)

funtion_plots_shadings(Wsort_act,Flvl_min_act(idxWsort_act),Flvl_max_act(idxWsort_act),ba)    
funtion_plots_shadings(Wsort_pass,Flvl_min_pass(idxWsort_pass),Flvl_max_pass(idxWsort_pass),bp)    

set(gca,'xscale','log')
ylim([0 0.8])
% set(gca,'yscale','log')
% ylim([0 1])
xlim([lowxlim highxlim])
% title('Net energy gain')
if scen==2
    legend([p1(1),p2(1),p3(1),p4(1),p5(1)],'f protists','f actives','f_c actives','f passives','f_c passives')
    legend boxoff
elseif scen==1
    ylabel({'Feeding level';'[-]'})
end

%%

if scen==1
axes(ha(5))
elseif scen==2
axes(ha(6))
end
scatter(0,0)
set(gca,'yscale','log')
set(gca,'xscale','log')
hold on
h1=semilogx(param.V,pred_C_P_mean,'color',cp,'linewidth', 1.5);
h2=semilogx(param.V,tot_mort_prtoists_mean,'--','color',cp,'linewidth', 1);
h3=semilogx(Wsort_act, pred_C_C_mean_act(idxWsort_act),'color',ba,'linewidth', 1.5);
h4=semilogx(Wsort_act, total_mortality_mean_act(idxWsort_act),'--','color',ba,'linewidth', 1);
h5=semilogx(Wsort_pass, pred_C_C_mean_pass(idxWsort_pass),'color',bp,'linewidth', 1.5);
h6=semilogx(Wsort_pass, total_mortality_mean_pass(idxWsort_pass),'--','color',bp,'linewidth', 1);


if scen==1
    ylabel({'Mortality rates';'[d^{-1}]'})
end

ylim([0 1])
xlim([lowxlim highxlim])
ylim([1e-5 1e1])
set(gca,'yscale','log')

xlabel('Body mass [\mugC]')
h7=loglog(m,mu,'k:')
if scen==2
legend([h1,h2,h3,h4,h5,h6,h7],{'pred. C on protists','total mortality protists','pred. C on act.',...
    'total mortality act.','pred. C on pass.','total mortality pass','Analytical solution'})
legend boxoff
end
plot_settings(fsize)

if scen==1
axes(ha(1))
title('\rho=0.005 d^{-1}')
xlim([1e-7 3e3])
set(gca,'Xticklabel',[])
yticks([1e-3 1e-1 1e1 10^3 1e5])


axes(ha(3))
xlim([1e-7 3e3])
set(gca,'Xticklabel',[])

axes(ha(5))
xlim([1e-7 3e3])
yticks([1e-4 10^-2 1e0])
yticklabels({'10^{-4}','10^{-2}','10^{0}'})
xticks([10^-5 1e-3 1e-1 1e1 1e3])


elseif scen==2
axes(ha(2))
title('\rho=0.05 d^{-1}')
xlim([1e-7 3e3])
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])

axes(ha(4))
xlim([1e-7 3e3])
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])

axes(ha(6))
xlim([1e-7 3e3])
set(gca,'Yticklabel',[])
xticks([1e-7 10^-5 1e-3 1e-1 1e1 1e3])

end

end