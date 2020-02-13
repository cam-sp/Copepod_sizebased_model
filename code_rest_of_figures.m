%% Figure R*
clearvars

m=logspace(-3,3);

ba=[29 41 81]./255;
bp=[87 160 211]./255;

effg=0.67;
% fc=0.18;
alpha_c=[0.011, 0.0052];
n=-0.25;
q=-0.25;
w=-0.25;
h=[1.37, 0.4];
d=0;%0.01;

PL=0.532*m.^0.32;
SR=1.801*PL-0.695;

sw=10^0.38.*PL.^0.93;
ratio=SR./sw;
tau=ratio;
tau(ratio<0)=0;
I=0.40*m.^(-0.25); %Ingestion rate non-calanoid
k_act=(0.16*m.^(-0.25));
k_pass=(0.048*m.^(-0.25)); %Metabolic cost thomas non-calanoids 0.0437
k=k_pass+tau.*(k_act-k_pass);

fc=0.18;
fc2=k./(effg.*I); 

ccolor=['k','r'];

cp= [236 197 68]./255;
% cp= [255 174 51]./255;

ba=[29 41 81]./255;
bp=[87 160 211]./255;

figure
for i=1:2

% Fstar=-(mumax.*(fc + (d.*m.^w)./(effg.*mumax)))./(alpha_c(i).*m.^q.*(fc + (d.*m.^w)./(effg.*mumax) - 1));

if i==1
Fstar=-(h(i).*m.^n.*(fc + (d.*m.^w)./(effg.*h(i).*m.^n)))./(alpha_c(i).*m.^q.*(fc + (d.*m.^w)./(effg.*h(i).*m.^n) - 1));
semilogx(m,Fstar,'color',ba,'linewidth',2)
elseif i==2
Fstar=-(h(i).*m.^n.*(fc2 + (d.*m.^w)./(effg.*h(i).*m.^n)))./(alpha_c(i).*m.^q.*(fc2 + (d.*m.^w)./(effg.*h(i).*m.^n) - 1));
semilogx(m,Fstar,'color',bp,'linewidth',2) 
end
hold on
end

Fs=-(h.*fc+d./effg)./(alpha_c.*(fc+d./(effg.*h)-1))


% clearvars

m=logspace(-7,0);

effg=1;
fc=0.2;
alpha_c=0.0024;
n=-0.25;
q=-0.25;
w=-0.25;
V=2.96e7.*m.^1.221;
rho=0.024.*V.^1.1;
Q=0.032.*V.^0.76;
mu_inf=4.7.*V.^-0.26;
mu_max=mu_inf.*rho./(mu_inf.*Q + rho);
mumax=mu_max;

fsize=10;
plot_settings(fsize)
% mumax=(0.11*m.^-0.25); 
d=0;%0.01;

ccolor=['k','r'];
% cp= [0.5882    0.5882    0.5882];



% figure
for i=1
Fstar=-(mumax.*(fc + (d.*m.^w)./(effg.*mumax)))./(alpha_c(i).*m.^q.*(fc + (d.*m.^w)./(effg.*mumax) - 1));

if i==1
semilogx(m,Fstar,'color',cp,'linewidth',1.5)
elseif i==2
semilogx(m,Fstar,'r') 
end
hold on
end

ylim([0 50])
xlabel('Organism body mass [\mugC]')
ylabel('E^* [\mugC L^{-1}]')
legend('Active copepods','Passive copepods','Protists')
legend boxoff
% saveas(gcf,'figure_Rstar','epsc')

Fs=-(h.*fc+d./effg)./(alpha_c.*(fc+d./(effg.*h)-1))


%% who wins Rmax
 clearvars
 

ba=[29 41 81]./255;
bp=[87 160 211]./255;
% bp=[128 188 247]./255;
 
m=logspace(-3,3,400);
F=linspace(10,100,500);
% f00=linspace(0.1,1,500);
% F=70;
a=linspace(0.1,0.9,500);
a=0.7;
af=[0.3,0.7];


x0=0;
y0=0;
width=18;
height=10;
grey=[0.5, 0.5, 0.5];
yellow=[244, 185, 66]./255;

cp= [0.5882    0.5882    0.5882];

fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');


st=1;
for j=1:2
    
    a=af(j);
    a3=a/4;%/3;

    
PL=0.532*m.^0.32;
SR=1.801*PL-0.695;

sw=10^0.38.*PL.^0.93;
ratio=SR./sw;
tau=ratio;
tau(ratio<0)=0;
h=1.37;
h2=0.4;
n=-0.25;
z=100;%62.*m.^0.34;%200;
m0=m/z;
m0(m0<1e-4)=1e-3;
% Recruitment efficiencies
eff_egg_a=0.5;%0.37;
eff_egg_p=0.5;%0.74;

a2=a3+tau.*(a-a3);

% Reproduction efficiencies
eff_r_a=0.5;%0.42;
eff_r_p=0.5;%0.58+tau.*(0.42-0.58);

erep=eff_r_a*eff_egg_a;
erep_p=eff_r_p*eff_egg_p;

I_act=h*m.^(-0.25); %Ingestion rate
k_act=(0.16*m.^(-0.25)); %Metabolic cost
clearance_act=(0.011*m.^(-0.25)); %Clearance rate
d_act=(0.00005.*m.^(-0.25)); %mortality rate 0.006

I_pass=h2*m.^(-0.25); %Ingestion rate non-calanoid
k_pass=(0.048*m.^(-0.25)); %Metabolic cost thomas non-calanoids 0.0437
clearance_pass=(0.0052*m.^(-0.25)); %Clearance rate
d_pass=(0.00005.*m.^(-0.25));%(0.007.*param.W(param.ind_act).^0); %mortality rate 0.003

tau2=tau;
tau2(:)=0;

Im=I_pass+tau2.*(I_act-I_pass);
km=k_pass+tau.*(k_act-k_pass);
alpha_cm=clearance_pass+tau2.*(clearance_act-clearance_pass);
d_cm=d_pass+tau.*(d_act-d_pass);

eff=0.67;
fc_act=k_act./(eff.*I_act); 
fc_pass=km./(eff.*Im); 

if length(F)>1
for i=1:length(F)
f0=clearance_act.*F(i)./(clearance_act.*F(i) + I_act);
% f0=f00(i);
% f02=f00(i);
f02=alpha_cm.*F(i)./(alpha_cm.*F(i) + Im);
A=0.67*h*(f0-fc_act);
A2=0.67*h2*(f02-fc_pass);
rmax(:,i)=A.*m.^n*n/(z.^n-1)*((1-a)*log(z)+log(erep./a));
rmax2(:,i)=A2.*m.^n.*n./(z.^n-1).*((1-a2).*log(z)+log(erep_p./(a2)));

flvl_a(:,i)=f0;
flvl_p(:,i)=f02;
end

% figure
% subplot(2,1,1)
% contourf(m,F,rmax')
% colorbar
% set(gca,'xscale','log')
% caxis([0 0.6])
% 
% subplot(2,1,2)
% contourf(m,F,rmax2')
% colorbar
% set(gca,'xscale','log')
% caxis([0 0.6])

rwin=zeros(size(rmax));
rwin(rmax>rmax2)=1;
rwin(rmax<0)=-1;
rwin(rmax2<0)=-1;

idx=find(F>25 & F<30);
lwid=1.5;

subplot(2,2,st)
semilogx(m,rmax(:,end),'color',ba,'linewidth',lwid)
hold on
semilogx(m,rmax(:,idx(end)),'--','color',ba,'linewidth',lwid)
semilogx(m,rmax2(:,end),'color',bp,'linewidth',lwid)
semilogx(m,rmax2(:,idx(end)),'--','color',bp,'linewidth',lwid)
ylim([1e-3 10^0])
yticks([1e-3 1e-2 1e-1 1e0])
yticklabels({'10^{-3}','10^{-2}','10^{-1}'})
xlim([1e-1 1e3])
xticks([1e-1 1e0 1e1 1e2 1e3])
xticklabels({'10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}'})

set(gca,'yscale','log')
xlabel('Copepod adult body-mass [\mugC]')
ylabel('r_{max} [d^{-1}]')
% legend('F=100\mugC L^{-1}', 'F=30\mugC L^{-1}','F=100\mugC L^{-1}', 'F=30\mugC L^{-1}')
columnlegend(2, {'Active cops. E=100\mugC L^{-1}', 'Active cops. E=30\mugC L^{-1}','Passive cops. E=100\mugC L^{-1}', 'Passive cops. E=30\mugC L^{-1}'})
legend boxoff

fsize=10;
plot_settings(fsize)

st=st+1;

subplot(2,2,st)
contourf(m,F,rwin')
set(gca,'xscale','log')
% colorbar
% colormap(inferno(3))
colormap([1 1 1; 1 0 0; 0 0 0])
colormap([1,1,1; bp; ba])
xlabel('Copepod adult body-mass [\mugC]')
ylabel('Prey conc. [\mugC L^{-1}]')
cbh = colorbar ; %Create Colorbar
 cbh.Ticks = -1:1;%linspace(-1,1,3) ; %Create 8 ticks from zero to 1
 cbh.TickLabels = {'no growth','Passive','Active'} ;  
xlim([1e-1 1e3])
st=st+1;
xlim([1e-1 1e3])
xticks([1e-1 1e0 1e1 1e2 1e3])
xticklabels({'10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}'})


else
for i=1:length(a)
f0=clearance_act.*F./(clearance_act.*F + I_act);
f02=alpha_cm.*F./(alpha_cm.*F + Im);
A=0.67*h*(f0-fc_act);
A2=0.67*h2*(f02-fc_pass);
rmax(:,i)=A.*m.^n*n/(z.^n-1)*((1-a(i))*log(z)+log(erep./a(i)));
rmax2(:,i)=A2.*m.^n.*n./(z.^n-1).*((1-a(i)).*log(z)+log(erep_p./a(i)));
end

figure
subplot(2,2,2)
contourf(m,a,rmax')
colorbar
set(gca,'xscale','log')
caxis([0 0.6])

subplot(2,2,4)
contourf(m,a,rmax2')
colorbar
set(gca,'xscale','log')
caxis([0 0.6])

rwin=zeros(size(rmax));
rwin(rmax<rmax2)=1;
rwin(rmax<0)=-1;
rwin(rmax2<0)=-1;


subplot(2,2,1)
contourf(m,a,rwin')
set(gca,'xscale','log')
% colorbar
colormap(viridis(3))
xlabel('Copepod adult body-mass [\mugC]')
ylabel('Prey concentration [\mugC L^{-1}]')
cbh = colorbar ; %Create Colorbar
 cbh.Ticks = linspace(-1,1,3) ; %Create 8 ticks from zero to 1
 cbh.TickLabels = {'no growth','Active','passive'} ;  

idx=find(F>25 & F<30);
lwid=1.5;

subplot(2,2,3)
semilogx(m,rmax(:,end),'k','linewidth',lwid)
hold on
semilogx(m,rmax(:,idx(end)),'k--','linewidth',lwid)
semilogx(m,rmax2(:,end),'r','linewidth',lwid)
semilogx(m,rmax2(:,idx(end)),'r--','linewidth',lwid)
% ylim([0 0.7])
% set(gca,'yscale','log')
xlabel('Copepod adult body-mass [\mugC]')
ylabel('r_{max} [d^{-1}]')
legend('F=100\mugC L^{-1}', 'F=30\mugC L^{-1}')
legend boxoff
xlim([1e-3 1e3])

fsize=10;
plot_settings(fsize)

end

end

% saveas(gcf,'figure_rmax_win','epsc')
subplot(2,2,1)
set(gca,'xTickLabel',[]);
subplot(2,2,2)
set(gca,'xTickLabel',[]);
fsize=10;
plot_settings(fsize)
set(findall(gcf,'-property','FontSize'),'FontSize',fsize)

%% development time
clearvars


ba=[29 41 81]./255;
bp=[87 160 211]./255;
% bp=[128 188 247]./255;

x0=0;
y0=0;
width=18;
height=8;
grey=[0.5, 0.5, 0.5];
yellow=[244, 185, 66]./255;

cp= [0.5882    0.5882    0.5882];

fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');




m=logspace(-1,3);
t=linspace(0,300);
z=100;%62.*m.^0.34;%200160;
m0=m./z;
h=1.37;
h2=0.40;
f0=linspace(0.3,1,2);
fc=0.17;
n=-1/4;



data1 = importdata('data_development_broadcast.csv'); %data from kiørboe and sabatini
data2 = importdata('data_development_sac.csv');

m_data=[data1(:,1); data2(:,1)];
dev_data=[data1(:,2); data2(:,2)];

m_data_a=data1(:,1);
m_data_p=data2(:,1);

dev_data_a=data1(:,2);
dev_data_p=data2(:,2);

% figure
subplot(1,2,1)
for i=1:length(f0)
    
    A=0.67*h*(f0(i)-fc);
    ta=(1-z.^(-n)).*m.^(-n)./(A.*n);
    
if i==1
    plot(m,ta,'--','color',ba,'linewidth',1.5)
else
    plot(m,ta,'color',ba,'linewidth',1.5)
end
    hold on
legendInfo{i} = ['f_0 = ' num2str(f0(i))]; 
end


for i=1:length(f0)

    A2=0.67*h2*(f0(i)-fc);
    ta2=(1-z.^(-n)).*m.^(-n)./(A2.*n);
if i==1    
    plot(m,ta2,'--','color',bp,'linewidth',1.5)
else
    plot(m,ta2,'color',bp,'linewidth',1.5)
end
    hold on
    legendInfo{length(f0)+i} = ['f_0 = ' num2str(f0(i))]; 
end

% scatter(m_data,dev_data,12,'markerfacecolor', 'k','markeredgecolor','k','MarkerFaceAlpha',0.5)
scatter(m_data_a,dev_data_a,12,'markerfacecolor', ba,'markeredgecolor',ba,'MarkerFaceAlpha',0.5)
scatter(m_data_p,dev_data_p,12,'markerfacecolor', bp,'markeredgecolor',bp,'MarkerFaceAlpha',0.5)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Adult copepods body mass [\mugC]')
ylabel('Development time [d]')
xlim([min(m) max(m)])
xticks([1e-1 1e0 1e1 1e2 1e3])
xticklabels({'10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}'})
fsize=10;
plot_settings(fsize)
% legend(legendInfo)
columnlegend(2,{'f_0=0.3 Active copepods','f_0=1 Active copepods','f_0=0.3 Passive copepods','f_0=1 Passive copepods'})
legend boxoff

% saveas(gcf,'figure_development','epsc')
% size at age
clearvars -except ba bp

m=100;
t=linspace(0,100);
z=100;%62.*m.^0.34;%200;
m0=m/z;
h=1.37;
h2=0.40;
f0=linspace(0.3,1,2);
fc=0.17;
n=-1/4;

subplot(1,2,2)
for i=1:length(f0)
    
    A=0.7*h*(f0(i)-fc);
    mt=(m0.^(-n)-A*n*t).^(-1/n);
    if i==1
        plot(t,mt,'--','color',ba,'linewidth',1.5)
    else
        plot(t,mt,'color',ba,'linewidth',1.5)
    end
    hold on
    legendInfo2{i} = ['f0 = ' num2str(f0(i))]; 
end

for i=1:length(f0)
    
    A2=0.67*h2*(f0(i)-fc);
    mt2=(m0.^(-n)-A2*n*t).^(-1/n);
    if i==1
        plot(t,mt2,'--','color',bp,'linewidth',1.5)
    else
        plot(t,mt2,'color',bp,'linewidth',1.5)
    end
    hold on
    legendInfo2{length(f0)+i} = ['f0 = ' num2str(f0(i))]; 
end
plot(t,mt*0+m,'k:')
legendInfo2{7} = 'Adult size'; 
% legend(legendInfo2)
set(gca,'yscale','log')
xlabel('Age [d]')
ylabel('Body mass [\mugC]')
ylim([1e0 1e2])

fsize=10;
plot_settings(fsize)


% %% sheldon appendix
% 
% % cp= [0.5882    0.5882    0.5882];
% cp= [236 197 68]./255;
% % cp= [255 174 51]./255;
% 
% ba=[29 41 81]./255;
% bp=[87 160 211]./255;
% cp= [236 197 68]./255;
% 
% if scen==1
% x0=0;
% y0=0;
% width=17;
% height=15;
% grey=[0.5, 0.5, 0.5];
% yellow=[244, 185, 66]./255;
% 
% figure('Units','centimeters',...
% 'Position',[x0 y0 width height],...
% 'PaperPositionMode','auto');
% 
% ha = tight_subplot(3,2,[.01 .01],[.1 .05],[.15 0.25])
% 
% axes(ha(1))
%     
% elseif scen==2
%     axes(ha(2))
% end
% 
% 
% %%
% if scen==1
% figure
% subplot(1,2,1)
% else
%     subplot(1,2,2)
% end
% Pspec(Pspec<1e-6)=NaN;
% scatter(0,0)
% set(gca,'yscale','log')
% set(gca,'xscale','log')
% hold on
% p1=loglog(param.V,Pspec,'o-','color',cp,'linewidth', 1.5,'markersize',2)
% funtion_plots_shadings(param.V,Pmin,Pmax,'k')
% if ontogeny_switch==1
% st1=1;
% st2=param.nbr_stages;
% for i=1:param.C_sp_act
% p2=loglog(param.Wvec(st1:st2),Cspec(st1:st2),'color',ba,'linewidth', 1.5)
% funtion_plots_shadings(param.Wvec(st1:st2),Cmin(st1:st2),Cmax(st1:st2),ba)
% 
% st1=st1+param.nbr_stages;
% st2=st2+param.nbr_stages;
% end
% st1=param.nbr_stages.*param.C_sp_act+1;
% st2=param.nbr_stages.*param.C_sp_act+param.nbr_stages;
% for i=1:param.C_sp_act
% p3=loglog(param.Wvec(st1:st2),Cspec(st1:st2),'color',bp,'linewidth', 1.5)
% funtion_plots_shadings(param.Wvec(st1:st2),Cmin(st1:st2),Cmax(st1:st2),bp)
% st1=st1+param.nbr_stages;
% st2=st2+param.nbr_stages;
% end
% 
% 
% else
%     loglog(param.Wvec(param.ind_a(1:param.C_sp_act)),Cmean(param.ind_a(1:param.C_sp_act))./deltaC(param.ind_a(1:param.C_sp_act)),'k','linewidth', 1.5)
%     loglog(param.Wvec(param.ind_a(param.C_sp_act+1:end)),Cmean(param.ind_a(param.C_sp_act+1:end))./deltaC(param.ind_a(param.C_sp_act+1:end)),'r','linewidth', 1.5)
%     
% end
% 
% % loglog(param.Wvec,1e6*param.D*param.Wvec.^(-2),'k--','linewidth',1.5)
% % loglog(logspace(-7,4),1e6*param.D*(logspace(-7,4)).^(-2),'k--','linewidth',1)
% ylim([1e-5 1e7])
% xlim([lowxlim highxlim])
% % title(['D= ',num2str(param.D),' '])
% ylabel({'Biomass spectrum';'[\mugC L^{-1} \mugC^{-1}]'})
% % title('Sheldon spectrum')
% % grid on
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%
% f0=mean([Flvl_mean_act, Flvl_mean_pass]);
% %     f0=0.3;
%     beta=mean([param.beta_act, param.beta_pass]);
%     sigma=mean([param.sigma_act,param.sigma_pass]);
%     m=logspace(-8,4);
%     q=-0.25;
%     n=-0.25;
%     v=mean([0.011, 0.0052]);
%     h=mean([1.37, 0.4]).*2.^((temp(1)-15)/10);
% 
%     fncs_analytical=function_analytical_solutions();
%     Nc=fncs_analytical.community_spectrum(f0,beta,sigma,m,h,v,q,n);
%     mu=fncs_analytical.predation_mortality(f0,beta,sigma,m,h,n);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% p4=loglog(m,Nc.*m.^2*1000,'k:')
% % p4=loglog(m,Nc.*m,'k:','linewidth',1.5)
% legend([p1(1),p2(1),p3(1),p4(1)],'Protists','Active cops.','Passive cops.','Analytical spectrum')
% legend boxoff
% ylim([1e0 1e6])
% 
% fsize=10;
% plot_settings(fsize)
% set(findall(gcf,'-property','FontSize'),'FontSize',fsize)