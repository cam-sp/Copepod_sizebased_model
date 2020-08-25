% Code for the paper "A generic size- and trait-based model of plankton
% communities" 2020 
%
%This code has the plots for the figures with the analytical solutions 
% Structure of the code:
    % 1- Plots for the R* (in appendix)
    % 2- Figure 4 in the paper
    % 3- Figure 3 in the paper

%
% Any questions related to this code can be addressed to me:
% cam.serra90@gmail.com
%
% Camila Serra-Pompei 25/08/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Figure E*
clearvars

%size range
m=logspace(-3,3);

%parameters
effg=0.67;
alpha_c=[0.011, 0.0052];
n=-0.25;
q=-0.25;
w=-0.25;
h=[1.37, 0.4];
d=0;
% PL=0.532*m.^0.32;
% SR=1.801*PL-0.695;
% sw=10^0.38.*PL.^0.93;
% ratio=SR./sw;
% tau=ratio;
% tau(ratio<0)=0;
I=0.40*m.^(-0.25); 
k_act=(0.16*m.^(-0.25));
k_pass=(0.048*m.^(-0.25)); 
k_pass(m>5)=NaN;
k=k_pass;
fc=0.18;
fc2=k./(effg.*I); 

%define colors
cp= [236 197 68]./255;
ba=[29 41 81]./255;
bp=[87 160 211]./255;

figure
for i=1:2

if i==1 %actives
    Fstar=-(h(i).*m.^n.*(fc + (d.*m.^w)./(effg.*h(i).*m.^n)))./(alpha_c(i).*m.^q.*(fc + (d.*m.^w)./(effg.*h(i).*m.^n) - 1));
    semilogx(m,Fstar,'color',ba,'linewidth',2)
elseif i==2 %passives
    Fstar=-(h(i).*m.^n.*(fc2 + (d.*m.^w)./(effg.*h(i).*m.^n)))./(alpha_c(i).*m.^q.*(fc2 + (d.*m.^w)./(effg.*h(i).*m.^n) - 1));
    semilogx(m,Fstar,'color',bp,'linewidth',2) 
end
hold on
end

% protists

m=logspace(-7,-1);

effg=1;
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
fcp=0.2;
Ip=10^-0.19.*24./1000.*(m./1000).^-0.33;
kp=0.2.*mumax;
alpha=alpha_c.*m.^q;

fsize=10;
plot_settings(fsize)
d=0;

ccolor=['k','r'];

Fstar=Ip./((alpha.*Ip./kp)-alpha);
semilogx(m,Fstar,'color',cp,'linewidth',1.5)

ylim([0 50])
xlabel('Organism body mass [\mugC]')
ylabel('E^* [\mugC L^{-1}]')
legend('Active copepods','Passive copepods','Protists')
legend boxoff


%% Figure 4 in the paper -  who wins Rmax
 clearvars
 

ba=[29 41 81]./255;
bp=[87 160 211]./255;
 
m=logspace(-3,3,400);
F=linspace(10,100,500);
% a=linspace(0.1,0.9,500);
a=0.7;
af=[0.3,0.7];


x0=0;
y0=0;
width=18;
height=10;


fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');


st=1;
for j=1:2
    
    a=af(j);
    a3=a/5;


h=1.37;
h2=0.4;
n=-0.25;
z=100;
m0=m/z;
m0(m0<1e-4)=1e-3;
a2=a/5;%+tau.*(a-a3);
erep=0.25;
erep_p=0.25;
I_act=h*m.^(-0.25); %Ingestion rate
k_act=(0.16*m.^(-0.25)); %Metabolic cost
clearance_act=(0.011*m.^(-0.25)); %Clearance rate
d_act=(0.00005.*m.^(-0.25)); %mortality rate 0.006
I_pass=h2*m.^(-0.25); %Ingestion rate non-calanoid
k_pass=(0.048*m.^(-0.25)); %Metabolic cost thomas non-calanoids 0.0437
clearance_pass=(0.0052*m.^(-0.25)); %Clearance rate
d_pass=(0.00005.*m.^(-0.25));%(0.007.*param.W(param.ind_act).^0); %mortality rate 0.003

eff=0.67;
fc_act=k_act./(eff.*I_act); 
fc_pass=k_pass./(eff.*I_pass); 


for i=1:length(F)
f0=clearance_act.*F(i)./(clearance_act.*F(i) + I_act);
f02=clearance_pass.*F(i)./(clearance_pass.*F(i) + I_pass);

A=0.67*h*(f0-fc_act);
A2=0.67*h2*(f02-fc_pass);
A2(m>5)=-10;
rmax(:,i)=A.*m.^n*n/(z.^n-1)*((1-a)*log(z)+log(erep./a));
rmax2(:,i)=A2.*m.^n.*n./(z.^n-1).*((1-a2).*log(z)+log(erep_p./(a2)));

flvl_a(:,i)=f0;
flvl_p(:,i)=f02;
end

rwin=zeros(size(rmax)); % 0 is passive feeders
rwin(rmax>rmax2)=1; % 1 is active feeders
rwin(rmax<0)=-1; % -1 is death

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
columnlegend(2, {'Active cops. E=100\mugC L^{-1}', 'Active cops. E=30\mugC L^{-1}','Passive cops. E=100\mugC L^{-1}', 'Passive cops. E=30\mugC L^{-1}'})
legend boxoff

fsize=10;
plot_settings(fsize)

st=st+1;

subplot(2,2,st)
contourf(m,F,rwin')
set(gca,'xscale','log')
colormap([1,1,1; bp; ba])
xlabel('Copepod adult body-mass [\mugC]')
ylabel('Prey conc. [\mugC L^{-1}]')
cbh = colorbar ; 
 cbh.Ticks = -1:1;
 cbh.TickLabels = {'no growth','Passive','Active'} ;  
 caxis([-1 1])
xlim([1e-1 1e3])
st=st+1;
xlim([1e-1 1e3])
xticks([1e-1 1e0 1e1 1e2 1e3])
xticklabels({'10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}'})

end

subplot(2,2,1)
set(gca,'xTickLabel',[]);
subplot(2,2,2)
set(gca,'xTickLabel',[]);
fsize=10;
plot_settings(fsize)
set(findall(gcf,'-property','FontSize'),'FontSize',fsize)

%% Figure 3 in the paper - development time
clearvars


ba=[29 41 81]./255;
bp=[87 160 211]./255;

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
z=100;
m0=m./z;
h=1.37;
h2=0.40;
f0=linspace(0.3,1,2);
fc=0.17;
n=-1/4;

data = importdata('data_development_times.xlsx'); %data from kiørboe and sabatini
dev_time=data.data(:,1);
ma=data.data(:,2);
act=data.data(:,3);
mix=data.data(:,4);
pass=data.data(:,5);

m_data_a=ma(act==1);
m_data_p=ma(pass==1);
m_data_mix=ma(mix==1);

dev_data_a=dev_time(act==1);
dev_data_p=dev_time(pass==1);
dev_data_mix=dev_time(mix==1);

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
    ta2(m>5)=NaN;
if i==1    
    plot(m,ta2,'--','color',bp,'linewidth',1.5)
else
    plot(m,ta2,'color',bp,'linewidth',1.5)
end
    hold on
    legendInfo{length(f0)+i} = ['f_0 = ' num2str(f0(i))]; 
end

scatter(m_data_a,dev_data_a,12,'markerfacecolor', ba,'markeredgecolor',ba,'MarkerFaceAlpha',0.5)
scatter(m_data_p,dev_data_p,12,'markerfacecolor', bp,'markeredgecolor',bp,'MarkerFaceAlpha',0.5)
scatter(m_data_mix,dev_data_mix,12,'markerfacecolor', 'none','markeredgecolor',ba,'MarkerFaceAlpha',0.5)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Adult copepods body mass [\mugC]')
ylabel('Development time [d]')
ylim([10^0.5 1e3])
xlim([min(m) max(m)])
xticks([1e-1 1e0 1e1 1e2 1e3])
xticklabels({'10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}'})
fsize=10;
plot_settings(fsize)
columnlegend(2,{'f_0=0.3 Active copepods','f_0=1 Active copepods','f_0=0.3 Passive copepods','f_0=1 Passive copepods'})
legend boxoff

clearvars -except ba bp

m=100;
t=linspace(0,80);
z=100;
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


plot(t,mt*0+m,'k:')
legendInfo2{7} = 'Adult size'; 
set(gca,'yscale','log')
xlabel('Age [d]')
ylabel('Body mass [\mugC]')
ylim([1e0 1e2])

fsize=10;
plot_settings(fsize)

