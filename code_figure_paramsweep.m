clearvars

%load run file
load('workspace_paramsweep.mat')

sensit_param_log=(sensit_param);
xlog=1;

%group protists biomass in size groups
Gyear=Pmean;

gg1=sum(Gyear(:,param.V<=1e-5 ),2);
gg2=sum(Gyear(:,param.V>1e-5 & param.V<=1e-3),2);
gg3=sum(Gyear(:,param.V>1e-3),2);
gg=cat(2,gg1,gg2,gg3);

%fecla pellets fluxes
idxstart=length(fluxcell{1})+1-length(Ncell{1});
flux=zeros(1,length(sensit_param_log));
flux1000=zeros(1,length(sensit_param_log));
NPPtot=zeros(1,length(sensit_param_log));
FPPtot=zeros(1,length(sensit_param_log));
for i=1:length(sensit_param_log)
    fluxx=fluxcell{i};
    flux(i)=mean(sum(fluxx(idxstart:end,:),2));
    fluxx1000=fluxx(idxstart:end,:).*exp(-param.remin./param.sink_fp.*960);
    flux1000(i)=mean(sum(fluxx1000,2));
    NPPtot(i)=mean(sum(NPPcell{i},2))*40;
    FPPtot(i)=mean(sum(NPPcell{i},2))*40;
end

%Copepods stuff
Cmean2=Cmean;
Cmean2(Cmean<1e-2)=NaN;

cnb=1;


%start figure
x0=0;
y0=0;
width=16;
height=15;
grey=[0.5, 0.5, 0.5];
yellow=[244, 185, 66]./255;

figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');

% ha = tight_subplot(3,2,[.01 .1],[.4 .05],[.1 .1])
ha = tight_subplot(4,1,[.01 .1],[.1 .02],[.15 0.25])

axes(ha(1))
 scatter(0,0)
 hold on

cmap=flip(inferno(length(gg(1,:))+1));
lwid=linspace(0.5,2.5,length(cmap));
lwid(:)=2;
h1=plot(sensit_param_log,gg,'-','linewidth',1.5);
for i=1:length(cmap)-1
  set(h1(i) ,'Color', cmap(1+i,:));
    set(h1(i) ,'linewidth', lwid(1+i));
end
hold on
ylim([0 120])
set(gca,'xscale','log')
xlim([min(sensit_param) max(sensit_param)])
ylabel('[mgC m^{-3}]')
set(gca,'Xticklabel',[])
legend boxoff
h2=plot(sensit_param_log,Nmean,'k:','linewidth',1.5);
set(gca,'Xticklabel',[])
legp=legend([h1;h2],{'m<10^{-5}','10^{-5}\leq m<10^{-3}','10^{-3}\leq m','N'});
ylim([0 150])


axes(ha(2))
 scatter(0,0)
 hold on
ifi=inferno(6);
ifi=ifi(2:end,:);

cmap=flip(inferno(param.C_sp/2+1));
lwid=linspace(0.5,2.5,length(cmap));
lwid(:)=2;
h1=plot(sensit_param_log,Cmean2(:,param.ind_a(1:param.C_sp_act)),'-');
hold on
for i=1:length(cmap)-1
  set(h1(i) ,'markeredgecolor', cmap(1+i,:));
  set(h1(i) ,'markerfacecolor', cmap(1+i,:));
  set(h1(i) ,'color', cmap(1+i,:));
  set(h1(i) ,'linewidth', lwid(1+i));
  set(h1(i) ,'markersize', 1);
  legendInfo{i}  = sprintf('m_a= %9.1E', param.Wvec(param.ind_a(i)));
end
funtion_plots_shadings(sensit_param_log,Cmin(:,param.ind_a(1:param.C_sp_act)),...
    Cmax(:,param.ind_a(1:param.C_sp_act)),fillcolor)
ylabel('[mg m^{-3}]')
if xlog==1
set(gca,'xscale','log')
end
set(gca,'yscale','log')
ylim([1e-2 1e2])
xlim([min(sensit_param) max(sensit_param)])
set(gca,'Xticklabel',[])
yticks([1e-2 1e-1 1e0 1e1])
legend(h1,legendInfo)
legend boxoff
set(gca,'Xticklabel',[])



axes(ha(3))
 scatter(0,0)
 hold on
cmap=flip(inferno(param.C_sp/2+1));
lwid=linspace(0.5,2.5,length(cmap));
lwid(:)=2;
h2=plot(sensit_param_log,Cmean2(:,param.ind_a(param.C_sp_act+1:end)),'-','markersize',5); %./param.Wvec(param.ind_a(param.C_sp_act+1:end))'*1000
for i=1:length(cmap)-1
  set(h2(i) ,'markeredgecolor', cmap(1+i,:));
  set(h2(i) ,'markerfacecolor', cmap(1+i,:));
   set(h2(i) ,'color', cmap(1+i,:));
  set(h2(i) ,'linewidth', lwid(1+i));
  legendInfo{i} = ['m_a= ' num2str(param.Wvec(param.ind_a(i)))]; 
end
hold on
fillcolor='k';
funtion_plots_shadings(sensit_param_log,Cmin(:,param.ind_a(param.C_sp_act+1:end)),...
    Cmax(:,param.ind_a(param.C_sp_act+1:end)),fillcolor)   
set(gca,'Xticklabel',[])
if xlog==1
set(gca,'xscale','log')
end
ylabel('[mg m^{-3}]')
set(gca,'yscale','log')
ylim([1e-2 1e2])
xlim([min(sensit_param) max(sensit_param)])
fsize=10;
plot_settings(fsize)
xlabel('Input of nitrogen \rho [d^{-1}]')
set(gca,'Xticklabel',[])
yticks([1e-2 1e-1 1e0 1e1])


cp= [0.5882    0.5882    0.5882];
axes(ha(4))
 scatter(0,0)
 hold on
  yyaxis left

h1=plot(sensit_param_log,flux,'k-','linewidth', 1.5)
hold on

set(gca,'xscale','log')

ylabel('mgC m^{-2} d^{-1}')
xlim([min(sensit_param) max(sensit_param)])
xlabel('Input of nitrogen \rho [d^{-1}]')
yticks([0 50 100])

fluxnan=flux;
fluxnan(sensit_param_log<0.003)=NaN;
yyaxis right
h2=plot(sensit_param_log,flux1000./fluxnan,'-','Color',cp,'linewidth', 1.5)
hold on
h3=plot(sensit_param_log,fluxnan./FPPtot,':','Color',cp,'linewidth', 1.5)

ylim([0 0.7])
ylabel('Fraction [-]')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cp;


legend([h1 h2 h3],{'Flux_{ML}','Flux_{1000}:Flux_{ML}','Flux_{ML}:FPP'})


legend boxoff

set(gca,'xTick',[1e-3, 1e-2, 1e-1]); 
set(gca,'xtickLabels',{'10^{-3}','10^{-2}','10^{-1}'})

fsize=10;
set(findall(gcf,'-property','FontSize'),'FontSize',fsize)