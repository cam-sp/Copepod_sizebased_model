%% Plots seasonality paper

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
% months = ['M';'J';'S';'D'];
% months = ['Ap';'Jl';'Oc';'Ja'];
bluemap=(brewermap(6,'YlGnBu'));
% bluemap=flip(viridis(6));
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

ha = tight_subplot(4,2,[.01 .2],[.1 .02],[.1 .15])

axes(ha(1))
yyaxis left
% subplot(4,rownbs,1)
h =bar(t(end-365*yr:end-yrend),gg,1,'stacked','EdgeColor','none');
plots=get(gca, 'Children');
hold on
% for i=13:12:length(tinterval)-1
% plot([tinterval(i), tinterval(i)], [1e-5 1e8],'k--')
% end
% ylim([0 max(sum(gg,2))])
% title('Protists')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
% set(gca,'XTick',t(end-365*yr)+30.5*3:30.5*3:t(end)-yrend); 
set(gca,'XtickLabels',months)
ylabel({'Protists'; '[mg m^{-3}]'})
set(h,{'FaceColor'},{bluemap(1,:);bluemap(2,:);bluemap(3,:);bluemap(4,:);bluemap(5,:);bluemap(6,:)});
xlabel('Month')
% legend(plots(2, 1), {'Second line', 'First line'});
% legf=legend(h,{'           m \leq 10^{-6}','10^{-6} < m \leq 10^{-5}','10^{-5} < m \leq 10^{-4}','10^{-4} < m \leq 10^{-3}','10^{-3} < m \leq 10^{-2}','10^{-2} < m'})
legend boxoff
% ylim([0 7000])
plot_settings(fsize)
set(gca,'Xticklabel',[])
set(gca,'YTick',[0 25 50 75 100]);
yyaxis right
plot(t(end-365*yr:end-yrend),N(end-365*yr:end-yrend),':','color',cp,'linewidth',1.5)
ax = gca;
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cp;
ylabel('Nitrogen [\mugN L^{-1}]')

legf=legend(plots([1 2 3 4 5 6]),flip({'           m < 10^{-6}','10^{-6} \leq m < 10^{-5}','10^{-5} \leq m < 10^{-4}','10^{-4} \leq m < 10^{-3}','10^{-3} \leq m < 10^{-2}','10^{-2} \leq m'}))
title(legf,{'Legend panel a';'Protiosts mass ranges'})
% legf=legend(h,flip({'           m < 10^{-6}','10^{-6} \leq m < 10^{-5}','10^{-5} \leq m < 10^{-4}','10^{-4} \leq m < 10^{-3}','10^{-3} \leq m < 10^{-2}','10^{-2} \leq m'}))
legf.FontSize = 8;
 

sc1=1e-2;
sc2=1e-1;
sc3=1e0;
sc4=1e1;
sc5=1e2;

Cyear=C(end-365*yr:end-yrend,param.ind_act);%*1000./param.Wvec(param.ind_act)';
cc1_a=sum(Cyear(:,param.W(param.ind_act)<sc1),2);
cc2_a=sum(Cyear(:,param.W(param.ind_act)>=sc1 & param.W(param.ind_act)<sc2),2);
cc3_a=sum(Cyear(:,param.W(param.ind_act)>=sc2 & param.W(param.ind_act)<sc3),2);
cc4_a=sum(Cyear(:,param.W(param.ind_act)>=sc3 & param.W(param.ind_act)<sc4),2);
cc5_a=sum(Cyear(:,param.W(param.ind_act)>=sc4 & param.W(param.ind_act)<sc5),2);
cc6_a=sum(Cyear(:,param.W(param.ind_act)>=sc5),2);

Cyear_p=C(end-365*yr:end-yrend,param.ind_pass);%*1000./param.W(param.ind_pass);
cc1_p=sum(Cyear_p(:,param.W(param.ind_pass)<sc1),2);
cc2_p=sum(Cyear_p(:,param.W(param.ind_pass)>=sc1 & param.W(param.ind_pass)<sc2),2);
cc3_p=sum(Cyear_p(:,param.W(param.ind_pass)>=sc2 & param.W(param.ind_pass)<sc3),2);
cc4_p=sum(Cyear_p(:,param.W(param.ind_pass)>=sc3 & param.W(param.ind_pass)<sc4),2);
cc5_p=sum(Cyear_p(:,param.W(param.ind_pass)>=sc4 & param.W(param.ind_pass)<sc5),2);
cc6_p=sum(Cyear_p(:,param.W(param.ind_pass)>=sc5),2);

cc_a=cat(2,cc1_a,cc2_a,cc3_a,cc4_a,cc5_a,cc6_a);%.*zmld(end-365*yr:end-yrend);%./sum(Cyear,2);
cc_p=cat(2,cc1_p,cc2_p,cc3_p,cc4_p,cc5_p,cc6_p);%.*zmld(end-365*yr:end-yrend);%./sum(Cyear_p,2);

% bluemap=flip(brewermap(6,'Accent'));
% bluemap=(viridis(6));
set(0,'DefaultAxesColorOrder',bluemap);
% subplot(3,2,1)
% axes(ha(2))
axes(ha(3))
% subplot(4,rownbs,2)
bar(t(end-365*yr:end-yrend),cc_a,1,'stacked')
hold on
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
% set(gca,'XTick',t(end-365*yr)+30.5*3:30.5*3:t(end)-yrend); 
for i=13:12:length(tinterval)-1
plot([tinterval(i), tinterval(i)], [1e-5 1e8],'k--')
end
set(gca,'XtickLabels',months)
set(gca,'YTick',[0 10 30 50]); 
ylabel({'Active copepods';'[mg m^{-3}]'})
% set(gca,'yscale','log')
% ylim([0 max(sum(cc_a,2))])
% grid on
plot_settings(fsize)
set(gca,'Xticklabel',[])

% bluemap=brewermap(6,'OrRd');
% set(0,'DefaultAxesColorOrder',bluemap);
% subplot(3,2,2)
% subplot(4,rownbs,3)
axes(ha(5))
% axes(ha(3))
bar(t(end-365*yr:end-yrend),cc_p,1,'stacked')
hold on
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
% set(gca,'XTick',t(end-365*yr)+30.5*3:30.5*3:t(end)-yrend); 
for i=13:12:length(tinterval)-1
plot([tinterval(i), tinterval(i)], [1e-5 1e8],'k--')
end
set(gca,'XtickLabels',months)
ylabel({'Passive copepods';'[mg m^{-3}]'})
plots=get(gca, 'Children');
% legend(plots(2, 1), {'Second line', 'First line'});
% legend('m \leq 10^{-1} ','10^{-1} < m \leq 10^{0}','10^{0} < m \leq 10^{1}','10^{1} < m \leq 10^{2}','10^{2} < m')
legf=legend(plots([1 2 3 4 5 6]),flip({'m < 10^{-2} ','10^{-2} \leq m < 10^{-1}','10^{-1} \leq m < 10^{0}','10^{0} \leq m < 10^{1}','10^{1} \leq m < 10^{2}','10^{2} \leq m'}));
% legf=legend({'m \leq 10^{-2} ','10^{-2} < m \leq 10^{-1}','10^{-1} < m \leq 10^{1}','10^{1} < m \leq 10^{2}','10^{2} < m \leq 10^{3}','10^{3} < m'});

legend boxoff
legf.FontSize = 8;
title(legf,{'Legend panels b and d';'Copepods mass ranges'})
% set(gca,'yscale','log')
set(gca,'Xticklabel',[])
set(gca,'YTick',[0 5 10]); 


% ylim([0 max(sum(cc_p,2))])
plot_settings(fsize)


% NPPtot=sum(NPP_diag,2).*zmld;
NPPtot=sum(FPP_prod_diag,2);

axes(ha(7))
% subplot(4,rownbs,4)
figure
yyaxis left
plot(t(end-365*yr:end-yrend),fluxmld(end-365*yr:end-yrend),'k','linewidth',1.5) %./zmld(end-365*yr:end-yrend)
hold on
plot(t(end-365*yr:end-yrend),flux1000(end-365*yr:end-yrend),'k--','linewidth',1.5)
ylabel({'Fecal pellets flux';'[mg m^{-2} d^{-2}]'})
set(gca,'YTick',[0 0.5 1]);
% semilogy(t(end-365*yr:end-yrend),N(end-365*yr:end-yrend),'k:')
yyaxis right
plot(t(end-365*yr:end-yrend),fluxmld(end-365*yr:end-yrend)./NPPtot(end-365*yr:end-yrend),'color',cp,'linewidth',1.5)
% hold on
% plot(t(end-365*yr:end-yrend),flux1000(end-365*yr:end-yrend)./NPPtot(end-365*yr:end-yrend),'--','color',cp,'linewidth',1.5)
% plot(t(end-365*yr:end-yrend),flux1000(end-365*yr:end-yrend)./fluxmld(end-365*yr:end-yrend),'color',cp,'linewidth',1.5)
ylim([0.05 0.35])
ylabel({'Fraction [-]'})
% ylim([0 0.4])
% set(gca,'yscale','log')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
plot_settings(fsize)
xlim([t(end-365*yr) t(end)-yrend])
legend({'MLD', '1000m', 'ML:1000m'},'FontSize',8)
legend boxoff
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cp;
plot_settings(fsize)
xlabel('Month')

deltaC=param.deltaC;
deltaC(end,:)=param.deltaC(param.ind_a-1);

axes(ha(4))
i=4;%2;
st1=(i-1)*param.nbr_stages+1;
st2=st1-1+param.nbr_stages;
contourf(t(end-365*yr:end-yrend),param.Wvec(st1:st2),log10(C(end-365*yr:end-yrend,st1:st2)'./(deltaC(:,i))),20,'edgecolor','none')
% surface(t(end-365*yr:end-yrend),param.Wvec(st1:st2),log10(C(end-365*yr:end-yrend,st1:st2)'./(deltaC(:,i))))
% shading interp
set(gca,'yscale','log')
yticks([1e-1 1e0])
% yticklabels({'10^{-2}','10^{-1}'})
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
bmap=(brewermap(10,'YlGnBu'));
% bmap=flip(brewermap(10,'Blues'));
colormap((bmap))
axis tight
% caxis([-2 1])
cbh=colorbar;
cbh.Ticks = [-6 -4 -2 0  4] ; %Create 8 ticks from zero to 1
cbh.TickLabels = {'10^{-6}','10^{-4}','10^{-2}','10^{0}','10^{2}','10^{4}'} ;
cbh.Label.String = 'Biom. spectrum [L^{-1}]';
ylabel({'Body mass';'[\mugC]'})
plot_settings(fsize)
set(gca,'Xticklabel',[])


axes(ha(6))
yyaxis right
semilogy(t(end-365*yr:end-yrend),reproduction_diag(end-365*yr:end-yrend,i),'o','color',cp,'markerfacecolor',cp,'markersize',1)
% ylim([1e-3 1e1])
yticks([1e-3 1e-1 1e1])
yticklabels({'10^{-3}','10^{-1}','10^{1}'})
ylabel({'Specific Rep. rate ';'[d^{-1}]'})
yyaxis left
plot(t(end-365*yr:end-yrend),reproduction_diag(end-365*yr:end-yrend,i).*C(end-365*yr:end-yrend,param.ind_a(i)),'ko','markerfacecolor','k','markersize',2)
% ylim([0 0.1])
xlim([t(end-365*yr) t(end)-yrend])
yticks([0.005 0.01 0.015])
% yticklabels({'0.05','0.1'})
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cp;
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
% set(gca,'XtickLabels',months)
ylabel({'Rep. rate ';'[mgC m^{-3} d^{-1}]'})
plot_settings(fsize)
set(gca,'Xticklabel',[])


[Wsort, sortind]=sort(param.Wvec(param.ind_act));

axes(ha(8))
bmap=(brewermap(10,'YlGnBu'));
% bmap(jet(10))
bmap=flip(viridis(20));
contourf(t(end-365*yr:end-yrend),Wsort,log10(pred_C_C_diag(end-365*yr:end-yrend,sortind)'),50,'edgecolor','none')
% surface(t(end-365*yr:end-yrend),param.Wvec(param.ind_act),log10(pred_C_C_diag(end-365*yr:end-yrend,param.ind_act)'))
% shading flat
cbh=colorbar;
cbh.Ticks = [-4 -3 -2 -1] ; %Create 8 ticks from zero to 1
cbh.TickLabels = {'10^{-4}','10^{-3}','10^{-2}','10^{-1}'} ;
cbh.Label.String = 'Predation [d^{-1}]';
set(gca,'yscale','log')
% title('Log_{10} Predation from Copepods on active cops.')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
set(gca,'yTick',[1e-3, 1e-1, 1e1, 1e3]); 
% set(gca,'ytickLabels',{'10^{-3}','10^{0}','10^{3}'})
xlim([t(end-365*yr) t(end)-yrend])
ylim([min(param.Wvec) max(param.Wvec)])
colormap(flip(bmap))
caxis([-4 -0.5])
ylabel({'Body mass';'[\mugC]'})
plot_settings(fsize)
xlabel('Month')

% print -depsc2 figure_seasonal.eps


%%


%% Plots seasonality paper

fluxmld=sum(detrit_flux_diag,2);
flux1000=sum(detrit_flux_diag.*exp(-param.remin./param.sink_fp.*(1000-zmld)),2);

fsize=9;
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
% months = ['M';'J';'S';'D'];
% months = ['Ap';'Jl';'Oc';'Ja'];
bluemap=(brewermap(6,'YlGnBu'));
% bluemap=flip(viridis(6));
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

ha = tight_subplot(3,1,[.01 .2],[.1 .02],[.1 .3])

axes(ha(1))
yyaxis left
% subplot(4,rownbs,1)
h =bar(t(end-365*yr:end-yrend),gg,1,'stacked','EdgeColor','none');
plots=get(gca, 'Children');
hold on
for i=13:12:length(tinterval)-1
plot([tinterval(i), tinterval(i)], [1e-5 1e8],'k--')
end
ylim([0 max(sum(gg,2))])
% title('Protists')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
% set(gca,'XTick',t(end-365*yr)+30.5*3:30.5*3:t(end)-yrend); 
set(gca,'XtickLabels',months)
ylabel('[mg m^{-2}]')
set(h,{'FaceColor'},{bluemap(1,:);bluemap(2,:);bluemap(3,:);bluemap(4,:);bluemap(5,:);bluemap(6,:)});
xlabel('Month')
% legend(plots(2, 1), {'Second line', 'First line'});
legf=legend({'           m \leq 10^{-6}','10^{-6} < m \leq 10^{-5}','10^{-5} < m \leq 10^{-4}','10^{-4} < m \leq 10^{-3}','10^{-3} < m \leq 10^{-2}','10^{-2} < m'})
title(legf,{'Legend panel A';'Protiosts mass ranges'})
legf.FontSize = 8;
legend boxoff
% ylim([0 7000])
plot_settings(fsize)
set(gca,'Xticklabel',[])
yyaxis right
plot(t(end-365*yr:end-yrend),N(end-365*yr:end-yrend),':','color',cp,'linewidth',1.5)
ax = gca;
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cp;
ylabel('Nitrogen [\mugN L^{-1}]')

sc1=1e-2;
sc2=1e-1;
sc3=1e0;
sc4=1e1;
sc5=1e2;

Cyear=C(end-365*yr:end-yrend,param.ind_act);%*1000./param.Wvec(param.ind_act)';
cc1_a=sum(Cyear(:,param.W(param.ind_act)<sc1),2);
cc2_a=sum(Cyear(:,param.W(param.ind_act)>=sc1 & param.W(param.ind_act)<sc2),2);
cc3_a=sum(Cyear(:,param.W(param.ind_act)>=sc2 & param.W(param.ind_act)<sc3),2);
cc4_a=sum(Cyear(:,param.W(param.ind_act)>=sc3 & param.W(param.ind_act)<sc4),2);
cc5_a=sum(Cyear(:,param.W(param.ind_act)>=sc4 & param.W(param.ind_act)<sc5),2);
cc6_a=sum(Cyear(:,param.W(param.ind_act)>=sc5),2);

Cyear_p=C(end-365*yr:end-yrend,param.ind_pass);%*1000./param.W(param.ind_pass);
cc1_p=sum(Cyear_p(:,param.W(param.ind_pass)<sc1),2);
cc2_p=sum(Cyear_p(:,param.W(param.ind_pass)>=sc1 & param.W(param.ind_pass)<sc2),2);
cc3_p=sum(Cyear_p(:,param.W(param.ind_pass)>=sc2 & param.W(param.ind_pass)<sc3),2);
cc4_p=sum(Cyear_p(:,param.W(param.ind_pass)>=sc3 & param.W(param.ind_pass)<sc4),2);
cc5_p=sum(Cyear_p(:,param.W(param.ind_pass)>=sc4 & param.W(param.ind_pass)<sc5),2);
cc6_p=sum(Cyear_p(:,param.W(param.ind_pass)>=sc5),2);

cc_a=cat(2,cc1_a,cc2_a,cc3_a,cc4_a,cc5_a,cc6_a);%.*zmld(end-365*yr:end-yrend);%./sum(Cyear,2);
cc_p=cat(2,cc1_p,cc2_p,cc3_p,cc4_p,cc5_p,cc6_p);%.*zmld(end-365*yr:end-yrend);%./sum(Cyear_p,2);

% bluemap=flip(brewermap(6,'Accent'));
% bluemap=(viridis(6));
set(0,'DefaultAxesColorOrder',bluemap);
% subplot(3,2,1)
% axes(ha(2))
axes(ha(2))
% subplot(4,rownbs,2)
bar(t(end-365*yr:end-yrend),cc_a,1,'stacked')
hold on
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
% set(gca,'XTick',t(end-365*yr)+30.5*3:30.5*3:t(end)-yrend); 
for i=13:12:length(tinterval)-1
plot([tinterval(i), tinterval(i)], [1e-5 1e8],'k--')
end
set(gca,'XtickLabels',months)
ylabel('[mg m^{-2}]')
% set(gca,'yscale','log')
ylim([0 max(sum(cc_a,2))])
% grid on
plot_settings(fsize)
set(gca,'Xticklabel',[])

% bluemap=brewermap(6,'OrRd');
% set(0,'DefaultAxesColorOrder',bluemap);
% subplot(3,2,2)
% subplot(4,rownbs,3)
axes(ha(3))
% axes(ha(3))
bar(t(end-365*yr:end-yrend),cc_p,1,'stacked')
hold on
% set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
plots=get(gca, 'Children');
% set(gca,'XTick',t(end-365*yr)+30.5*3:30.5*3:t(end)-yrend); 
set(gca,'XTick',t(end-365*yr):365:t(end)-yrend); 
for i=13:12:length(tinterval)-1
plot([tinterval(i), tinterval(i)], [1e-5 1e8],'k--')
end
ylim([0 max(sum(cc_p,2))])
set(gca,'XtickLabels',0:365:365*8)
ylabel('[mg m^{-2}]')
% legend(plots(2, 1), {'Second line', 'First line'});
% legend('m \leq 10^{-1} ','10^{-1} < m \leq 10^{0}','10^{0} < m \leq 10^{1}','10^{1} < m \leq 10^{2}','10^{2} < m')
legf=legend({'m \leq 10^{-2} ','10^{-2} < m \leq 10^{-1}','10^{-1} < m \leq 10^{1}','10^{1} < m \leq 10^{2}','10^{2} < m \leq 10^{3}','10^{3} < m'});
legend boxoff
legf.FontSize = 8;
title(legf,{'Legend panels B and D';'Copepods mass ranges'})
% set(gca,'yscale','log')
% set(gca,'Xticklabel',[])
xlabel('Days')

%%

cc_large=sum(cc_a(:,5:6),2)+sum(cc_p(:,5:6),2);

figure
subplot(1,2,1)
% plot(t(end-365*yr:end-yrend),NPPtot(end-365*yr:end-yrend),'k','linewidth',1.5)
plot(t(end-365*yr:end-yrend),(sum(cc_a,2)+sum(cc_p,2)).*zmld(end-365*yr:end-yrend),'k','linewidth',1.5)
xlabel('Month')
ylabel('mgC m^{-2} d^{-1}')
title('Total GPP in the mixed layer')
axis tight
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)


subplot(1,2,2)
% plot(t(end-365*yr:end-yrend),sum(pred_C_on_D_diag(end-365*yr:end-yrend,:).*D(end-365*yr:end-yrend,:),2).*zmld(end-365*yr:end-yrend),'k','linewidth',1.5)
plot(t(end-365*yr:end-yrend),sum(pred_C_on_D_diag(end-365*yr:end-yrend,:),2),'k','linewidth',1.5)
xlabel('Month')
axis tight
title('Specific consumption of FP')
ylabel('mgC m^{-2} d^{-1}')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)

fsize=10;
plot_settings(fsize)

%% integrated over the mixed layer
figure
subplot(3,2,1)
plot(t(end-365*yr:end-yrend),(sum(cc_a,2)+sum(cc_p,2)).*zmld(end-365*yr:end-yrend),'k','linewidth',1.5)
hold on
plot(t(end-365*yr:end-yrend),cc_large.*zmld(end-365*yr:end-yrend),'r','linewidth',1.5)
plot(t(end-365*yr:end-yrend),(1-(cc_large./(sum(cc_a,2)+sum(cc_p,2)))).*(sum(cc_a,2)+sum(cc_p,2)).*zmld(end-365*yr:end-yrend),'b','linewidth',1.5)
legend('Total','m>10\mugC','m<10\mugC')
legend boxoff
title('Integrated copepod biomass in mixed layer')
ylabel('mgC m^{-2}')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
xlim([t(end-365*yr) t(end)-yrend])

subplot(3,2,2)
plot(t(end-365*yr:end-yrend),cc_large./(sum(cc_a,2)+sum(cc_p,2)),'k','linewidth',1.5)
title('fraction of copepods larger than 10\mugC')
ylabel('[-]')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
xlim([t(end-365*yr) t(end)-yrend])

subplot(3,2,3)
plot(t(end-365*yr:end-yrend),sum(D(end-365*yr:end-yrend,:),2).*zmld(end-365*yr:end-yrend),'k','linewidth',1.5)
title('integrated biomass of fecal pellets in the mixed layer')
ylabel('mgC m^{-3}')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
xlim([t(end-365*yr) t(end)-yrend])

subplot(3,2,5)
plot(t(end-365*yr:end-yrend),fluxmld(end-365*yr:end-yrend)./NPPtot(end-365*yr:end-yrend),'k','linewidth',1.5)
title('flux_{ML}:FPP')
ylabel('[-]')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
xlim([t(end-365*yr) t(end)-yrend])

subplot(3,2,4)
plot(t(end-365*yr:end-yrend),fluxmld(end-365*yr:end-yrend),'k','linewidth',1.5)
hold on
plot(t(end-365*yr:end-yrend),flux1000(end-365*yr:end-yrend),'k--','linewidth',1.5)
legend('flux out of mixed layer','flux at 1000 m')
legend boxoff
title('total flux fecal pellets')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
xlim([t(end-365*yr) t(end)-yrend])

subplot(3,2,6)
plot(t(end-365*yr:end-yrend),flux1000(end-365*yr:end-yrend)./fluxmld(end-365*yr:end-yrend),'k','linewidth',1.5)
title('flux_{1000}:flux_{ML}')
ylabel('[-]')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
xlim([t(end-365*yr) t(end)-yrend])

plot_settings(10)
set(findall(gcf,'-property','FontSize'),'FontSize',10)
%%

figure
plot(t(end-365*yr:end-yrend),-zmld(end-365*yr:end-yrend),'k','linewidth',1.5)
title('zmld')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)

figure
plot(t(end-365*yr:end-yrend),temp(end-365*yr:end-yrend),'k','linewidth',1.5)
title('zmld')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)

figure
plot(t(end-365*yr:end-yrend),max(0,dzdt(end-365*yr:end-yrend)),'k','linewidth',1.5)
title('zmld')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)

figure
plot(t(end-365*yr:end-yrend),D(end-365*yr:end-yrend,:).*zmld(end-365*yr:end-yrend),'linewidth',1.5)
hold on
plot(t(end-365*yr:end-yrend),sum(D(end-365*yr:end-yrend,:),2).*zmld(end-365*yr:end-yrend),'r','linewidth',1.5)
legend
title('zmld')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
%%

figure

for i=1:20
    
    orderplots=[19 17 15 13 11 9 7 5 3 1 20 18 16 14 12 10 8 6 4 2];
   subplot(10,2,orderplots(i)) 
    
yyaxis right
semilogy(t(end-365*yr:end-yrend),reproduction_diag(end-365*yr:end-yrend,i),'o','color',cp,'markerfacecolor',cp,'markersize',1)
ylim([1e-3 1e-1])
yticks([1e-3 1e-2 1e-1])
yticklabels({'10^{-3}','10^{-2}','10^{-1}'})
hold on
yyaxis left
plot(t(end-365*yr:end-yrend),reproduction_diag(end-365*yr:end-yrend,i).*C(end-365*yr:end-yrend,param.ind_a(i)),'ko','markerfacecolor','k','markersize',2)
ylim([0 0.1])
xlim([t(end-365*yr) t(end)-yrend])
% yticks([0.05 0.1])
% yticklabels({'0.05','0.1'})
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cp;
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
if i==1
set(gca,'XtickLabels',months)
elseif i==11
    set(gca,'XtickLabels',months)
else
    set(gca,'XtickLabels',[])
end
alk=sprintf('m_a= %9.1E',param.Wa(i));
text(4.841e4,0.075,alk)
plot_settings(fsize)
 
    
    
end
ylabel({'Rep. rate ';'[mgC m^{-3} d^{-1}]'})
   yyaxis right
ylabel({'Specific Rep. rate ';'[mgC m^{-3} d^{-1}]'})
set(findall(gcf,'-property','FontSize'),'FontSize',10)