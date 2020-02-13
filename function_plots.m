function [plot_fncs]= function_plots()

        plot_fncs.plot_full_run       = @function_plot_full_run;
        plot_fncs.plot_annual         = @function_plot_annual;        
        plot_fncs.plot_annual_biomass = @function_plot_annual_biomass; 
        plot_fncs.plot_lat55          = @function_plot_lat55; 
        plot_fncs.plot_lat60          = @function_plot_lat60; 
        plot_fncs.plot_protists       = @function_plot_protists_stacked;
        plot_fncs.plot_chla           = @function_plot_chla;
        plot_fncs.plot_cohorts        = @function_plot_cohorts;
        plot_fncs.plot_diagnostics    = @function_plot_diagnostics;
        plot_fncs.plot_growth_rates   = @function_plot_growth_rates;
        plot_fncs.plot_mortalities    = @function_plot_mortalities;
        plot_fncs.plot_nu_cohorts     = @function_plot_nu_cohorts;
        plot_fncs.plot_reproduction   = @function_plot_reproduction;
        plot_fncs.plot_numberspec     = @function_plot_numberspec;
        plot_fncs.plot_all_stacked    = @function_plot_all_stacked;
        plot_fncs.plot_predations     = @function_plot_predations;
        plot_fncs.plot_all            = @function_plot_all;
        plot_fncs.plot_diagnostics_protists          = @function_diagnostics_protists;
        plot_fncs.plot_Flvlsources           = @function_Flvlsources;
end


function function_plot_full_run(t,N,P,C,param,fsize)

figure

subplot(4,1,1)
plot(t,N,'k','linewidth',1.5);
ylabel('N [\mugN L^{-1}]')

cmap=flip(viridis(param.nbr_P));
subplot(4,1,2)
h1=plot(t,P,'linewidth',1.5);
for i=1:length(cmap)
  set(h1(i) ,'Color', cmap(i,:));
end
legend('Location','northwest')
ylabel('P [\mugC L^{-1}]')

subplot(4,1,3)
cmap=flip(viridis(param.C_sp/2));
h1=plot(t,C(:,param.ind_a(1:param.C_sp_act)),'linewidth',1.5);
hold on
legend('Location','northwest')
% for i=1:length(cmap)
%   set(h1(i) ,'Color', cmap(i,:));
% end
ylabel('C [\mugC L^{-1}]')

subplot(4,1,4)
cmap=flip(viridis(param.C_sp/2));
h2=plot(t,C(:,param.ind_a(param.C_sp_act+1:end)),'linewidth',1.5);
legend('Location','northwest')
% for i=1:length(cmap)
%   set(h2(i) ,'Color', cmap(i,:));
% end
ylabel('C [\mugC L^{-1}]')

plot_settings(fsize)

end




function function_plot_annual(t,N,P,C,param,fsize,yr,yrend)



figure
subplot(4,1,1)
plot(t(end-365*yr:end-yrend),N(end-365*yr:end-yrend),'k','linewidth',1.5)
months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
set(gca,'XTick',t(end-365*yr):30.5:t(end)); set(gca,'XtickLabels',months)
title('N')
ylabel('[\mugN L^{-1}]')
plot_settings(fsize)

subplot(4,1,2)
cmap=flip(viridis(param.nbr_P));
for i=1:param.nbr_P
plot(t(end-365*yr:end-yrend:end),P(end-365*yr:end-yrend,i), 'color', cmap(i,:),'linewidth',1.5)
hold on
end
months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
set(gca,'XTick',t(end-365*yr)+365/4:365/4:t(end)-yrend); 
set(gca,'XtickLabels',months)
for j=1:param.nbr_P
    sizeG=param.V(j);
legendInfo2{j} = ['m = ' num2str(sizeG) '\mugC'];
end
title('Protists')
ylabel('[\mugC L^{-1}]')
legend(legendInfo2,'Location','northwest')
plot_settings(fsize)


subplot(4,1,3)
cmap=flip(viridis(param.C_sp/2));
h1=plot(t(end-365*yr:end-yrend),C(end-365*yr:end-yrend,param.ind_a(1:end/2)),'linewidth',1.5);
hold on
h2=plot(t(end-365*yr:end-yrend),C(end-365*yr:end-yrend,param.ind_a(end/2+1:end)),'--','linewidth',1.5);
months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
set(gca,'XTick',t(end-365*yr)+365/4:365/4:t(end)-yrend); 
set(gca,'XtickLabels',months)
for j=1:param.nbr_cops
    sizeC=param.Wa(j);
    legendInfo2{j} = ['m= ' num2str(sizeC) '\mugC'];
end
for i=1:length(cmap)
  set(h1(i) ,'Color', cmap(i,:));
  set(h2(i) ,'Color', cmap(i,:));
end
legend(legendInfo2,'Location','northwest')
title('Concentration adult copepods')
ylabel('[\mugC L^{-1}]')
plot_settings(fsize)

subplot(4,1,4)
h1=plot(t(end-365*yr:end-yrend),C(end-365*yr:end-yrend,param.ind_a(1:end/2)).*1000./param.Wvec(param.ind_a(1:end/2))','linewidth',1.5);
hold on
h2=plot(t(end-365*yr:end-yrend),C(end-365*yr:end-yrend,param.ind_a(end/2+1:end)).*1000./param.Wvec(param.ind_a(end/2+1:end))','--','linewidth',1.5);
months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
set(gca,'XTick',t(end-365*yr)+365/4:365/4:t(end)-yrend); set(gca,'XtickLabels',months)
for i=1:length(cmap)
  set(h1(i) ,'Color', cmap(i,:));
  set(h2(i) ,'Color', cmap(i,:));
end
title('Abundance Adult Copepod')
ylabel('[# m^{-3}]')
xlabel('Month')
plot_settings(fsize)        
end

function  function_plot_annual_biomass(t,N,P,C,param,fsize,yr,yrend)

figure
cmap=flip(viridis(length(param.ind_a(1:end/2))));
lwid=linspace(0.2,2.5,length(cmap));
lwid(:)=1.5;
subplot(2,1,1)
h1=plot(t(end-365*yr:end-yrend),C(end-365*yr:end-yrend,param.ind_a(1:end/2)),'linewidth',2);
hold on
h2=plot(t(end-365*yr:end-yrend),C(end-365*yr:end-yrend,param.ind_a(end/2+1:end)),'--','linewidth',2);
months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
set(gca,'XTick',t(end-365*yr)+365/4:365/4:t(end)-yrend); set(gca,'XtickLabels',months)
for j=1:param.nbr_cops
    sizeC=param.Wa(j);
legendInfo2{j} = ['m= ' num2str(sizeC) '\mugC'];
end
title('Biomass')
% xlim([t(end-365) t(end)])
ylabel('mgC m^{-3}')
for i=1:length(cmap)
  set(h1(i) ,'Color', cmap(i,:), 'linewidth', lwid(i));
  set(h2(i) ,'Color', cmap(i,:), 'linewidth', lwid(i));
end

legend(legendInfo2,'Location','northwest')
plot_settings(fsize)

subplot(2,1,2)
h1=semilogy(t(end-365*yr:end-yrend),C(end-365*yr:end-yrend,param.ind_a(1:end/2)).*1000./param.Wvec(param.ind_a(1:end/2))','linewidth',2);
hold on
h2=semilogy(t(end-365*yr:end-yrend),C(end-365*yr:end-yrend,param.ind_a(end/2+1:end)).*1000./param.Wvec(param.ind_a(end/2+1:end))','--','linewidth',2);
months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
set(gca,'XTick',t(end-365*yr)+365/4:365/4:t(end)-yrend); set(gca,'XtickLabels',months)
title('Abundance')
% xlim([t(end-365) t(end)])
ylabel('# m^{-3}')
ylim([1e-0 1e4])
xlabel('Month')
for i=1:length(cmap)
  set(h1(i) ,'Color', cmap(i,:), 'linewidth', lwid(i));
  set(h2(i) ,'Color', cmap(i,:), 'linewidth', lwid(i));
end
xlabel('Month')
plot_settings(fsize)

end


function  function_plot_lat55(t,N,P,C,param,yearsrun,fsize)

data_biomass = importdata('data_biomass_copepods.csv');
day_cop=data_biomass(:,1);
biom=data_biomass(:,2)./28;

data_oithona= importdata('data_oithona.csv');
day_oithona=data_oithona(:,1);
biom_oithona=data_oithona(:,2)./28;

yrs=5;
figure
% semilogx(t(end-365*yrs:end),sum(C(end-365*yrs:end,:),2).*zmld(end-365*yrs:end),'k:','linewidth',1.5)
hold on
plot(t(end-365*yrs:end),sum(C(end-365*yrs:end,:),2),'k','linewidth',2)
% plot(t(end-365*yrs:end),sum(C(end-365*yrs:end,param.ind_act),2)*30,'b','linewidth',1)
plot(t(end-365*yrs:end),sum(C(end-365*yrs:end,param.ind_pass),2),'r','linewidth',2)
scatter(day_cop+365*(yearsrun-1),biom,30,'markerfacecolor','k','markeredgecolor','k','MarkerFaceAlpha',0.5)
scatter(day_oithona+365*(yearsrun-1),biom_oithona,30,'markerfacecolor','r','markeredgecolor','r','MarkerFaceAlpha',0.5)
% ylim([100 max(sum(C(end-365*yrs:end,:),2)*30)])
% set(gca,'yscale','log')
months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
set(gca,'XTick',t(end-365*yrs):30.5:t(end)); set(gca,'XtickLabels',months)
% title(['mort ',num2str(param.m_cte),' '])
xlabel('Month')
ylabel('mgC m^{-3}')
legend({'Total biomass model','Passive biomass model','Total biomass data','Passive biomass data'})
set(gcf,'color','w');
xlim([t(end-365*yrs) t(end)])
plot_settings(fsize)


end

function  function_plot_lat60(t,N,P,C,param,fsize)

yrs=1;
figure
% semilogx(t(end-365*yrs:end),sum(C(end-365*yrs:end,:),2).*zmld(end-365*yrs:end),'k:','linewidth',1.5)
hold on
plot(t(end-365*yrs:end),sum(C(end-365*yrs:end,:),2),'k--','linewidth',1.5)
plot(t(end-365*yrs:end),sum(C(end-365*yrs:end,param.ind_act),2),'b','linewidth',1)
plot(t(end-365*yrs:end),sum(C(end-365*yrs:end,param.ind_pass),2),'r','linewidth',1)

months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
set(gca,'XTick',t(end-365*yrs):30.5:t(end)); set(gca,'XtickLabels',months)
% title(['mort ',num2str(param.m_cte),' '])
plot_settings(fsize)

end

function  function_plot_protists_stacked(t,N,P,C,param,fsize)

yr=10;%;
yrend=0;%yr-1;
tinterval=t(end-365*yr):30.5:t(end)-yrend;
Gyear=P(end-365*yr:end-yrend,:);
gg1=sum(Gyear(:,param.V<=1e-5),2);
gg2=sum(Gyear(:,param.V>1e-5 & param.V<=1e-4),2);
gg3=sum(Gyear(:,param.V>1e-4 & param.V<=1e-3),2);
gg4=sum(Gyear(:,param.V>1e-3 & param.V<=1e-2),2);
gg5=sum(Gyear(:,param.V>1e-2),2);

gg=cat(2,gg1,gg2,gg3,gg4,gg5);
months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
% bluemap=brewermap(5,'Blues');
bluemap=flip(viridis(5));
figure
h =bar(t(end-365*yr:end-yrend),gg,1,'stacked','EdgeColor','none');
hold on
for i=12:12:length(tinterval)-1
plot([tinterval(i), tinterval(i)], [1e-5 1e4],'k--')
end
title('Protists')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
ylabel('Biomass [mg m^{-3}]')
% set(gca,'yscale','log')
% set(h,{'FaceColor'},{bluemap(1,:);bluemap(2,:);'r';'r';bluemap(1,:)});
set(h,{'FaceColor'},{bluemap(1,:);bluemap(2,:);bluemap(3,:);bluemap(4,:);bluemap(5,:)});
% set(bar_handle(1),'FaceColor',bluemap(1,:)) ;
% set(bar_handle(:,2),'FaceColor','b',bluemap(2,:)) ;
% set(bar_handle(3),'FaceColor','y',bluemap(3,:)) ;
% set(bar_handle(4),'FaceColor','g',bluemap(4,:)) ;
grid on
xlabel('Month')
legend({'           m \leq 10^{-5}','10^{-5} < m \leq 10^{-4}','10^{-4} < m \leq 10^{-3}','10^{-3} < m \leq 10^{-2}','10^{-2} < m'})
plot_settings(fsize)

end

function function_plot_chla(t,N,P,C,param,fsize)
yr=1;
    figure
    plot(t(end-365*yr:end), sum(P(end-365*yr:end,:),2)./75,'k','linewidth',2)
    set(gca,'XTick',t(end-365*yr):30.5:t(end)); 
    months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
    set(gca,'XtickLabels',months)
    title('Chl-a ')
    ylabel('mg m^{-3}')
%     ylim([0 3])
    xlim([t(end-365*yr) t(end)])
    plot_settings(fsize)
    
end

function function_plot_cohorts(t,N,P,C,param,fsize,yr,yrend)

months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
tinterval=t(end-365*yr):30.5:t(end)-yrend;

viridis_cmap=viridis(10);
figure
st1=1;
st2=param.nbr_stages;
deltaC=param.deltaC;
% deltaC(end,:)=1e-1;

for i=1:param.nbr_cops
       orderplots=flip([2:2:param.C_sp_act*2,1:2:param.C_sp_act*2]);
%     subplot(2,param.nbr_cops/2,i)
    subplot(param.nbr_cops/2,2,orderplots(i))
% surface(t(end-365*yr:end-yrend),param.Wvec(st1:st2),log10(C(end-365*yr:end-yrend,st1:st2)'.*1000./(param.W(:,i).*deltaC(:,i))))
surface(t(end-365*yr:end-yrend),param.Wvec(st1:st2),log10(C(end-365*yr:end-yrend,st1:st2)'./(deltaC(:,i))))
shading interp
% imagesc(C(end-365:end,st1:st2)'./param.Wvec(st1:st2))
% set(gca,'Ydir','normal')
% set(gca,'zscale','log')
set(gca,'yscale','log')
set(gca,'XTick',tinterval); 
set(gca,'XtickLabels',months)
% colorbar
colormap(viridis_cmap)
st1=st2+1;
st2=st2+param.nbr_stages;
axis tight
caxis([-6 5])
end
subplot(param.nbr_cops/2,2,param.nbr_cops/2*2)
chb=colorbar('southoutside');
ylabel('Body mass [\mugC]')
set(gcf,'color','w');
chb.Label.String = 'Log_{10} Number spectrum [# m^{-3} mgC^{-1}]';
plot_settings(fsize)

end


function function_plot_diagnostics(t,N,P,C,param,fsize,fg_diag,Flvl_diag,pred_C_P_diag,pred_C_C_diag,pred_P_P_diag,dd_mortC_diag)

months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
viridis_cmap=viridis(40);

figure

subplot(3,3,1)
surface(t(end-365:end),param.V,fg_diag(end-365:end,:)')
shading interp
colorbar
set(gca,'yscale','log')
title('feeding lvl protists')
caxis([0 1])
set(gca,'XTick',t(end-365):30.5:t(end)); 
set(gca,'XtickLabels',months)
axis tight

subplot(3,3,2)
surface(t(end-365:end),param.Wvec(param.ind_act),Flvl_diag(end-365:end,param.ind_act)')
shading interp
colorbar
set(gca,'yscale','log')
title('feeding lvl active')
caxis([param.fc(1) 1])
set(gca,'XTick',t(end-365):30.5:t(end)); 
set(gca,'XtickLabels',months)
ylim([min(param.Wvec) max(param.Wvec)])
xlim([t(end-365) t(end)])

subplot(3,3,3)
surface(t(end-365:end),param.Wvec(param.ind_pass),Flvl_diag(end-365:end,param.ind_pass)')
shading interp
colorbar
title('fedding lvl passive')
set(gca,'yscale','log')
caxis([param.fc(1) 1])
set(gca,'XTick',t(end-365):30.5:t(end)); 
set(gca,'XtickLabels',months)
ylim([min(param.Wvec) max(param.Wvec)])
xlim([t(end-365) t(end)])


subplot(3,3,4)
surface(t(end-365:end),param.V,pred_C_P_diag(end-365:end,:)')
shading interp
colorbar
title('predat C on P')
set(gca,'yscale','log')
set(gca,'XTick',t(end-365):30.5:t(end)); 
set(gca,'XtickLabels',months)
axis tight

subplot(3,3,5)
surface(t(end-365:end),param.Wvec(param.ind_act),pred_C_C_diag(end-365:end,param.ind_act)')
shading interp
colorbar
set(gca,'yscale','log')
title('Predation from C on C_act')
set(gca,'XTick',t(end-365):30.5:t(end)); 
set(gca,'XtickLabels',months)
ylim([min(param.Wvec) max(param.Wvec)])
xlim([t(end-365) t(end)])

subplot(3,3,6)
surface(t(end-365:end),param.Wvec(param.ind_pass),pred_C_C_diag(end-365:end,param.ind_pass)')
shading interp
colorbar
set(gca,'yscale','log')
set(gca,'XTick',t(end-365):30.5:t(end)); 
set(gca,'XtickLabels',months)
title('Predation from C on C_pass')
ylim([min(param.Wvec) max(param.Wvec)])
xlim([t(end-365) t(end)])

subplot(3,3,7)
surface(t(end-365:end),param.V,pred_P_P_diag(end-365:end,:)')
shading interp
colorbar
set(gca,'yscale','log')
set(gca,'XTick',t(end-365):30.5:t(end)); 
set(gca,'XtickLabels',months)
axis tight

subplot(3,3,8)
surface(t(end-365:end),param.Wvec(param.ind_act),dd_mortC_diag(end-365:end,param.ind_act)')
shading interp
colorbar
set(gca,'yscale','log')
set(gca,'XTick',t(end-365):30.5:t(end)); 
set(gca,'XtickLabels',months)
title('dd_mort Act')
ylim([min(param.Wvec) max(param.Wvec)])
xlim([t(end-365) t(end)])


subplot(3,3,9)
surface(t(end-365:end),param.Wvec(param.ind_pass),dd_mortC_diag(end-365:end,param.ind_pass)')
shading interp
colorbar
set(gca,'yscale','log')
set(gca,'XTick',t(end-365):30.5:t(end)); 
title('dd_mort Pass')
set(gca,'XtickLabels',months)
title('density dependent mortality')
ylim([min(param.Wvec) max(param.Wvec)])
xlim([t(end-365) t(end)])

colormap(viridis_cmap)

plot_settings(fsize)

end


function function_plot_growth_rates(t,N,P,C,param,fsize,mu_diag,dPdt_diag,nu_diag,dCdt_diag)

months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
numap=brewermap(100,'PiYG');
figure
subplot(2,3,1)
surface(t(end-365:end),param.V,mu_diag(end-365:end,:)')
shading interp
colorbar
colormap(numap)
set(gca,'yscale','log')
title('mu')
caxis([-0.6 0.6])
set(gca,'XTick',t(end-365):30.5:t(end)); 
set(gca,'XtickLabels',months)
axis tight

subplot(2,3,4)
surface(t(end-365:end),param.V,dPdt_diag(end-365:end,:)'./P(end-365:end,:)')
shading interp
colorbar
colormap(numap)
set(gca,'yscale','log')
title('dPdt')
caxis([-0.6 0.6])
set(gca,'XTick',t(end-365):30.5:t(end)); 
set(gca,'XtickLabels',months)
axis tight

numinmax=0.15;

subplot(2,3,2)
surface(t(end-365:end),param.Wvec(param.ind_act),nu_diag(end-365:end,param.ind_act)')
shading interp
colorbar
colormap(numap)
set(gca,'yscale','log')
title('nu act')
caxis([-numinmax numinmax])
set(gca,'XTick',t(end-365):30.5:t(end)); 
set(gca,'XtickLabels',months)
ylim([min(param.Wvec) max(param.Wvec)])
xlim([t(end-365) t(end)])

subplot(2,3,3)
surface(t(end-365:end),param.Wvec(param.ind_pass),nu_diag(end-365:end,param.ind_pass)')
shading interp
colorbar
colormap(numap)
set(gca,'yscale','log')
title('nu pass')
caxis([-0.05 0.05])
set(gca,'XTick',t(end-365):30.5:t(end)); 
set(gca,'XtickLabels',months)
ylim([min(param.Wvec) max(param.Wvec)])
xlim([t(end-365) t(end)])

subplot(2,3,5)
surface(t(end-365:end),param.Wvec(param.ind_act),dCdt_diag(end-365:end,param.ind_act)'./C(end-365:end,param.ind_act)')
shading interp
colorbar
set(gca,'yscale','log')
title('dCdt act')
caxis([-numinmax numinmax])
set(gca,'XTick',t(end-365):30.5:t(end)); 
set(gca,'XtickLabels',months)
ylim([min(param.Wvec) max(param.Wvec)])
xlim([t(end-365) t(end)])

subplot(2,3,6)
surface(t(end-365:end),param.Wvec(param.ind_pass),(dCdt_diag(end-365:end,param.ind_pass)'./C(end-365:end,param.ind_pass)'))
shading interp
colorbar
title('dCdt pass')
set(gca,'yscale','log')
caxis([-0.01 0.01])
set(gca,'XTick',t(end-365):30.5:t(end)); 
set(gca,'XtickLabels',months)
ylim([min(param.Wvec) max(param.Wvec)])
xlim([t(end-365) t(end)])

plot_settings(fsize)


end


function function_plot_mortalities(t,N,P,C,param,fsize,starvation_diag,pred_C_C_diag, dd_mortC_diag,yr,yrend)

months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
tinterval=t(end-365*yr):30.5:t(end)-yrend;

main_mort=zeros(size(pred_C_C_diag));
for j=1:param.nbr_cops*param.nbr_stages
    for i=t(end)-365*yr:t(end-yrend)
        morts_vec=[starvation_diag(i,j),pred_C_C_diag(i,j), param.d_c(j), dd_mortC_diag(i,j)];
      [~,  main_mort(i,j)]=max(morts_vec);
    end        
end

bluemap2=brewermap(4,'Blues');
viridismap=viridis(4);
figure
st1=1;
st2=param.nbr_stages;

for i=1:param.nbr_cops
orderplots=flip([2:2:param.C_sp_act*2,1:2:param.C_sp_act*2]);
%     subplot(2,param.nbr_cops/2,i)
    subplot(param.nbr_cops/2,2,orderplots(i))
surface(t(end-365*yr:end-yrend),param.Wvec(st1:st2),main_mort(end-365*yr:end-yrend,st1:st2)')
shading interp
% imagesc(C(end-365:end,st1:st2)'./param.Wvec(st1:st2))
% set(gca,'Ydir','normal')
% set(gca,'zscale','log')
set(gca,'yscale','log')
set(gca,'XTick',tinterval); 
set(gca,'XtickLabels',months)
st1=st2+1;
st2=st2+param.nbr_stages;
axis tight
colormap(bluemap2)
caxis([1 4])
end
subplot(param.nbr_cops/2,2,param.nbr_cops/2*2)
cbh=colorbar('southoutside');
cbh.Ticks = [1 2 3 4] ; %Create 8 ticks from zero to 1
cbh.TickLabels = {'Starvation','Predation','Linear','Density dependen'} ;

plot_settings(fsize)

end


function function_plot_nu_cohorts(t,N,P,C,param,fsize,nu_diag)

viridis_cmap=viridis(10);
months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
figure
st1=1;
st2=param.nbr_stages;
deltaC=param.deltaC;
deltaC(end,:)=1;

for i=1:param.nbr_cops
       orderplots=flip([2:2:param.C_sp_act*2,1:2:param.C_sp_act*2]);
%     subplot(2,param.nbr_cops/2,i)
    subplot(param.nbr_cops/2,2,orderplots(i))
surface(t(end-365:end),param.Wvec(st1:st2),(nu_diag(end-365:end,st1:st2)'))
shading interp
% imagesc(C(end-365:end,st1:st2)'./param.Wvec(st1:st2))
% set(gca,'Ydir','normal')
% set(gca,'zscale','log')
set(gca,'yscale','log')
set(gca,'XTick',t(end-365):30.5:t(end)); 
set(gca,'XtickLabels',months)
colorbar
colormap(viridis_cmap)
st1=st2+1;
st2=st2+param.nbr_stages;
axis tight
caxis([-0.1 0.1])
end
% colormap(numap)
title('nu')
plot_settings(fsize)

end


function function_plot_reproduction(t,N,P,C,param,fsize,reproduction_diag,yr,yrend,mu_diag)

tinterval=t(end-365*yr):30.5:t(end)-yrend;
viridis_cmap=viridis(30);

% months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
months = ['A';'J';'O';'J'];
months = ['1/4';'1/2';'3/4';''];
months = ['0.25';'0.50';'0.75';''];
months = ['   ';'1/2';'   ';'   '];
figure
numap=brewermap(100,'PiYG');
ax1=subplot(3,1,1)
contourf(t(end-365*yr:end-yrend),param.V,mu_diag(end-365*yr:end-yrend,:)');
hold on
for i=13:12:length(tinterval)-1
plot([tinterval(i), tinterval(i)], [1e-8 1e8],'k--')
end
% shading interp
colorbar
colormap(ax1,numap)
set(gca,'yscale','log')
title('Division rate protists [d^{-1}]')
caxis([-0.6 0.6])
% set(gca,'XTick',t(end-365):30.5:t(end)); 
% set(gca,'XtickLabels',months)
set(gca,'XTick',t(end-365*yr)+365/4:365/4:t(end)-yrend); 
set(gca,'XtickLabels',months)
ylim([min(param.V) max(param.V)])
xlim([t(end-365*yr) t(end)-yrend])
ylabel({'Body mass'; '[\mugC]'})
set(gca,'yTick',[1e-5, 1e-4, 1e-3]); 
set(gca,'ytickLabels',{'10^{-5}','10^{-4}','10^{-3}'})

subplot(3,1,2)
contourf(t(end-365*yr:end-yrend),param.Wvec(param.ind_a(1:param.C_sp_act)),log10(reproduction_diag(end-365*yr:end-yrend,1:param.C_sp_act)'))
% shading interp
hold on
colorbar
set(gca,'yscale','log')
% set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
% set(gca,'XtickLabels',months)
set(gca,'XTick',t(end-365*yr)+365/4:365/4:t(end)-yrend); 
set(gca,'XtickLabels',months)
xlim([t(end-365*yr) t(end)-yrend])
title('Log_{10} Reproduction active copepods [d^{-1}]')
for i=13:12:length(tinterval)-1
plot([tinterval(i), tinterval(i)], [1e-5 1e8],'k--')
end
ylim([min(param.Wvec) max(param.Wvec)])
plot_settings(fsize)
% caxis([-4 -1])
ylabel({'Body mass'; '[\mugC]'})
set(gca,'yTick',[1e-3, 1e-0, 1e3]); 
set(gca,'ytickLabels',{'10^{-3}','10^{0}','10^{3}'})

subplot(3,1,3)
surface(t(end-365*yr:end-yrend),param.Wvec(param.ind_a(param.C_sp_act+1:end)),log10(reproduction_diag(end-365*yr:end-yrend,param.C_sp_act+1:end)'))
shading interp
hold on
colorbar
set(gca,'yscale','log')
% set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
% set(gca,'XtickLabels',months)
set(gca,'XTick',t(end-365*yr)+365/4:365/4:t(end)-yrend); 
set(gca,'XtickLabels',months)
xlim([t(end-365*yr) t(end)-yrend])
title('Log_{10} Reproduction passive copepods [d^{-1}]')
for i=13:12:length(tinterval)-1
plot([tinterval(i), tinterval(i)], [1e-5 1e8],'k--')
end
colormap(viridis_cmap)
ylim([min(param.Wvec) max(param.Wvec)])
caxis([-4 -1])
plot_settings(fsize)
set(gca,'yTick',[1e-3, 1e-0, 1e3]); 
set(gca,'ytickLabels',{'10^{-3}','10^{0}','10^{3}'})
ylabel({'Body mass'; '[\mugC]'})
xlabel('Fraction of year')


% figure
% % plot(t(end-365:end),reproduction_diag(end-365:end,1:end/2)./C(end-365:end,param.ind_a(1:end/2)),'linewidth',2)
% plot(t(end-365:end),reproduction_diag(end-365:end,1:end/2),'o','linewidth',2)
% hold on
% % plot(t(end-365:end),reproduction_diag(end-365:end,end/2+1:end)./C(end-365:end,param.ind_a(end/2+1:end)),'--','linewidth',2)
% plot(t(end-365:end),reproduction_diag(end-365:end,end/2+1:end),'*','linewidth',2)
% set(gca,'XTick',t(end-365):30.5:t(end)); 
% set(gca,'XtickLabels',months)
% for j=1:param.nbr_cops
%     sizeC=param.Wa(j);
% legendInfo2{j} = ['m= ' num2str(sizeC) '\mugC'];
% end
% xlabel('Month')
% ylabel('Reproduction [mgC m^{-3} d^{-1}]')
% set(gca,'FontSize',18)
% set(gcf,'color','w');
% xlim([t(end-365) t(end)])
% 
% legend(legendInfo2,'Location','northwest')

end


function function_plot_numberspec(t,N,P,C,param,fsize,fg_diag,Flvl_diag,total_mortality_diag,pred_C_C_diag,dd_mortC_diag,...
                                    pred_P_P_diag,pred_C_P_diag,ontogeny_switch,seasonal_switch,mback)
                      
                                
[Wsort_act, idxWsort_act]=sort(param.Wvec(param.ind_act));
[Wsort_pass, idxWsort_pass]=sort(param.Wvec(param.ind_pass));
                                
if seasonal_switch==0

Pmean=mean(P(ceil(t(end)*3/4):end,:),1).*1000./param.V;
Cmean=mean(C(ceil(t(end)*3/4):end,:),1).*1000./param.Wvec';

fg_mean=mean(fg_diag(ceil(t(end)*3/4):end,:),1);
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
mback_mean=mean(mback(ceil(t(end)*3/4):end,:),1);

Flvl_max_act=max(Flvl_diag(ceil(t(end)*3/4):end,param.ind_act));
Flvl_max_pass=max(Flvl_diag(ceil(t(end)*3/4):end,param.ind_pass));

Flvl_min_act=min(Flvl_diag(ceil(t(end)*3/4):end,param.ind_act));
Flvl_min_pass=min(Flvl_diag(ceil(t(end)*3/4):end,param.ind_pass));

else

Pmean=P.*1000./param.V;
Cmean=C.*1000./param.Wvec';

fg_mean=fg_diag;
Flvl_mean_act=Flvl_diag(param.ind_act);
Flvl_mean_pass=Flvl_diag(param.ind_pass);
total_mortality_mean_act=total_mortality_diag(param.ind_act);
total_mortality_mean_pass=total_mortality_diag(param.ind_pass);
pred_C_C_mean_act=pred_C_C_diag(param.ind_act);
pred_C_C_mean_pass=pred_C_C_diag(param.ind_pass);
ddmort_mean_act=dd_mortC_diag(param.ind_act);
ddmort_mean_pass=dd_mortC_diag(param.ind_pass);
pred_P_P_mean=pred_P_P_diag;
pred_C_P_mean=pred_C_P_diag;
    
    
end


deltaC=param.deltaC;
deltaC(end,:)=param.deltaC(param.ind_a-1);

lowxlim=1e-8;

figure
subplot(3,1,1)
loglog(param.V,Pmean./param.delta_V,'g','linewidth', 1.5)
hold on
if ontogeny_switch==1
st1=1;
st2=param.nbr_stages;
for i=1:param.C_sp_act
loglog(param.Wvec(st1:st2),Cmean(st1:st2)./deltaC(st1:st2),'k','linewidth', 1.5)
st1=st1+param.nbr_stages;
st2=st2+param.nbr_stages;
end
st1=param.nbr_stages.*param.C_sp_act+1;
st2=param.nbr_stages.*param.C_sp_act+param.nbr_stages;
for i=1:param.C_sp_act
loglog(param.Wvec(st1:st2),Cmean(st1:st2)./deltaC(st1:st2),'r','linewidth', 1.5)
st1=st1+param.nbr_stages;
st2=st2+param.nbr_stages;
end
else
    loglog(param.Wvec(param.ind_a(1:param.C_sp_act)),Cmean(param.ind_a(1:param.C_sp_act))./deltaC(param.ind_a(1:param.C_sp_act)),'k','linewidth', 1.5)
    loglog(param.Wvec(param.ind_a(param.C_sp_act+1:end)),Cmean(param.ind_a(param.C_sp_act+1:end))./deltaC(param.ind_a(param.C_sp_act+1:end)),'r','linewidth', 1.5)
    
end

% loglog(param.Wvec,1e6*param.D*param.Wvec.^(-2),'k--','linewidth',1.5)
loglog(logspace(-7,4),1e6*param.D*(logspace(-7,4)).^(-2),'k--','linewidth',1)
ylim([1e-5 1e20])
xlim([lowxlim 1e4])
grid on
title(['D= ',num2str(param.D),' '])
grid on

subplot(3,1,2)
plot(param.V, fg_mean,'g')
hold on
plot(Wsort_act, Flvl_mean_act(idxWsort_act),'k')
plot(Wsort_pass, Flvl_mean_pass(idxWsort_pass),'r')
fc_pass=param.fc(param.ind_pass);
plot(Wsort_pass, fc_pass(idxWsort_pass),'r--')
plot(Wsort_act, param.fc(param.ind_act),'k--')
if seasonal_switch==0
funtion_plots_shadings(Wsort_act,Flvl_min_act(idxWsort_act),Flvl_max_act(idxWsort_act),'k')    
funtion_plots_shadings(Wsort_pass,Flvl_min_pass(idxWsort_pass),Flvl_max_pass(idxWsort_pass),'r')   
end
set(gca,'xscale','log')
ylim([0 1])
xlim([lowxlim 1e4])

subplot(3,1,3)
loglog(param.V, pred_P_P_mean,'g')
hold on
loglog(param.V, mback_mean,'m')
loglog(param.V, pred_C_P_mean,'b')
loglog(Wsort_act, total_mortality_mean_act(idxWsort_act),'k')
loglog(Wsort_pass, total_mortality_mean_pass(idxWsort_pass),'r')
loglog(Wsort_act, pred_C_C_mean_act(idxWsort_act),'k--')
loglog(Wsort_act, ddmort_mean_act(idxWsort_act),'k:','linewidth', 1.5)
loglog(Wsort_pass, ddmort_mean_pass(idxWsort_pass),'r:','linewidth', 1.5)
loglog(Wsort_pass, pred_C_C_mean_pass(idxWsort_pass),'r--')
ylim([0 1])
xlim([lowxlim 1e4])
ylim([1e-5 1e1])
title(['mort ',num2str(param.m_coef),' '])
grid on

plot_settings(fsize)

end

function function_plot_all_stacked(t,N,P,C,param,fsize,yr,yrend,zmld)

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
% months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
months = ['M';'J';'S';'D'];
months = ['Ap';'Jl';'Oc';'Ja'];
% bluemap=brewermap(5,'Blues');
bluemap=flip(viridis(6));
figure

subplot(3,rownbs,1)
h =bar(t(end-365*yr:end-yrend),gg,1,'stacked','EdgeColor','none');
hold on
for i=13:12:length(tinterval)-1
plot([tinterval(i), tinterval(i)], [1e-5 1e8],'k--')
end
ylim([0 max(sum(gg,2))])
title('Protists')
% set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XTick',t(end-365*yr)+30.5*3:30.5*3:t(end)-yrend); 
set(gca,'XtickLabels',months)
ylabel('[mg m^{-2}]')
set(h,{'FaceColor'},{bluemap(1,:);bluemap(2,:);bluemap(3,:);bluemap(4,:);bluemap(5,:);bluemap(6,:)});
xlabel('Month')
% legend({'           m \leq 10^{-6}','10^{-6} < m \leq 10^{-5}','10^{-5} < m \leq 10^{-4}','10^{-4} < m \leq 10^{-3}','10^{-3} < m \leq 10^{-2}','10^{-2} < m'})
plots=get(gca, 'Children');
legf=legend(plots([1 2 3 4 5 6]),flip({'           m < 10^{-6}','10^{-6} \leq m < 10^{-5}','10^{-5} \leq m < 10^{-4}','10^{-4} \leq m < 10^{-3}','10^{-3} \leq m < 10^{-2}','10^{-2} \leq m'}))
title(legf,{'Legend panel a';'Protiosts mass ranges'})
legf.FontSize = 8;
legend boxoff
plot_settings(fsize)


%%%%%%%%%
%copepods

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

bluemap=brewermap(6,'OrRd');
set(0,'DefaultAxesColorOrder',bluemap);
% subplot(3,2,1)
% axes(ha(2))
subplot(3,rownbs,2)
bar(t(end-365*yr:end-yrend),cc_a,1,'stacked')
hold on
title('Active copepods')
% set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XTick',t(end-365*yr)+30.5*3:30.5*3:t(end)-yrend); 
for i=13:12:length(tinterval)-1
plot([tinterval(i), tinterval(i)], [1e-5 1e8],'k--')
end
set(gca,'XtickLabels',months)

ylabel('[mg m^{-2}]')
% set(gca,'yscale','log')
ylim([0 max(sum(cc_a,2))])
% grid on
plot_settings(fsize)


bluemap=brewermap(6,'OrRd');
set(0,'DefaultAxesColorOrder',bluemap);
% subplot(3,2,2)
subplot(3,rownbs,3)
% axes(ha(3))
bar(t(end-365*yr:end-yrend),cc_p,1,'stacked')
hold on
title('Passive copepods')
% set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XTick',t(end-365*yr)+30.5*3:30.5*3:t(end)-yrend); 
for i=13:12:length(tinterval)-1
plot([tinterval(i), tinterval(i)], [1e-5 1e8],'k--')
end
set(gca,'XtickLabels',months)
ylabel('[mg m^{-2}]')
% legend('m \leq 10^{-1} ','10^{-1} < m \leq 10^{0}','10^{0} < m \leq 10^{1}','10^{1} < m \leq 10^{2}','10^{2} < m')
% legend('m \leq 10^{-2} ','10^{-2} < m \leq 10^{-1}','10^{-1} < m \leq 10^{1}','10^{1} < m \leq 10^{2}','10^{2} < m \leq 10^{3}','10^{3} < m')
plots=get(gca, 'Children');
legf=legend(plots([1 2 3 4 5 6]),flip({'m < 10^{-2} ','10^{-2} \leq m < 10^{-1}','10^{-1} \leq m < 10^{0}','10^{0} \leq m < 10^{1}','10^{1} \leq m < 10^{2}','10^{2} \leq m'}));

legend boxoff
% set(gca,'yscale','log')

ylim([0 max(sum(cc_p,2))])
plot_settings(fsize)
% grid on

% Cyear=C(end-365*yr:end-yrend,param.ind_act)*1000./param.Wvec(param.ind_act)';
% cc1_a=sum(Cyear(:,param.W(param.ind_act)<=1e0),2);
% cc2_a=sum(Cyear(:,param.W(param.ind_act)>1e0 & param.W(param.ind_act)<=1e1),2);
% cc3_a=sum(Cyear(:,param.W(param.ind_act)>1e1 & param.W(param.ind_act)<=1e2),2);
% cc4_a=sum(Cyear(:,param.W(param.ind_act)>1e2 & param.W(param.ind_act)<=1e3),2);
% cc5_a=sum(Cyear(:,param.W(param.ind_act)>1e3),2);
% 
% Cyear_p=C(end-365*yr:end-yrend,param.ind_pass)*1000./param.W(param.ind_pass);
% cc1_p=sum(Cyear_p(:,param.W(param.ind_pass)<=1e0),2);
% cc2_p=sum(Cyear_p(:,param.W(param.ind_pass)>1e0 & param.W(param.ind_pass)<=1e1),2);
% cc3_p=sum(Cyear_p(:,param.W(param.ind_pass)>1e1 & param.W(param.ind_pass)<=1e2),2);
% cc4_p=sum(Cyear_p(:,param.W(param.ind_pass)>1e2 & param.W(param.ind_pass)<=1e3),2);
% cc5_p=sum(Cyear_p(:,param.W(param.ind_pass)>1e3),2);
% 
% cc_a=cat(2,cc1_a,cc2_a,cc3_a,cc4_a,cc5_a);
% cc_p=cat(2,cc1_p,cc2_p,cc3_p,cc4_p,cc5_p);
% 
% bluemap=brewermap(5,'OrRd');
% set(0,'DefaultAxesColorOrder',bluemap);
% % subplot(3,2,1)
% % axes(ha(2))
% subplot(3,rownbs,4)
% bar(t(end-365*yr:end-yrend),cc_a,1,'stacked')
% hold on
% title('Active copepods')
% set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
% for i=13:12:length(tinterval)-1
% plot([tinterval(i), tinterval(i)], [1e-5 1e8],'k--')
% end
% set(gca,'XtickLabels',months)
% ylabel('Abundance [# L^{-1}]')
% % set(gca,'yscale','log')
% ylim([0 max(sum(cc_a,2))])
% % set(gca,'yscale','log')
% 
% 
% bluemap=brewermap(5,'OrRd');
% set(0,'DefaultAxesColorOrder',bluemap);
% % subplot(3,2,2)
% subplot(3,rownbs,6)
% % axes(ha(3))
% bar(t(end-365*yr:end-yrend),cc_p,1,'stacked')
% hold on
% title('Ambush copepods')
% set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
% for i=13:12:length(tinterval)-1
% plot([tinterval(i), tinterval(i)], [1e-5 1e8],'k--')
% end
% set(gca,'XtickLabels',months)
% ylabel('Abundance [# L^{-1}]')
% % legend('m \leq 10^{-1} ','10^{-1} < m \leq 10^{0}','10^{0} < m \leq 10^{1}','10^{1} < m \leq 10^{2}','10^{2} < m')
% legend('m \leq 10^{0} ','10^{0} < m \leq 10^{1}','10^{1} < m \leq 10^{2}','10^{2} < m \leq 10^{3}','10^{3} < m')
% % ylim([0 5])
% ylim([0 max(sum(cc_p,2))])
% % set(gca,'yscale','log')

end


function function_plot_predations(t,N,P,C,param,fsize,pred_C_P_diag,pred_C_C_diag,pred_P_P_diag,dd_mortC_diag,yr,yrend)

months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
% months = ['   ';'1/2';'   ';'   '];
viridis_cmap=viridis(30);
tinterval=t(end-365*yr):30.5:t(end)-yrend;
figure
subplot(4,1,2)
% surface(t(end-365*yr:end-yrend),param.V,pred_C_P_diag(end-365*yr:end-yrend,:)')
% shading flat
contourf(t(end-365*yr:end-yrend),param.V,pred_C_P_diag(end-365*yr:end-yrend,:)',10,'edgecolor','none')

hold on
for i=13:12:length(tinterval)-1
plot([tinterval(i), tinterval(i)], [1e-8 1e8],'w--')
end
title('predat from copepods on Protists')
set(gca,'yscale','log')
set(gca,'XTick',tinterval); 
set(gca,'XtickLabels',months)
set(gca,'yTick',[1e-5, 1e-4, 1e-3]); 
set(gca,'ytickLabels',{'10^{-5}','10^{-4}','10^{-3}'})
ylim([min(param.V) max(param.V)])
xlim([t(end-365*yr) t(end)-yrend])
ylabel({'Body mass'; '[\mugC]'})
colorbar

pred_C_A=pred_C_C_diag(end-365*yr:end-yrend,param.ind_act)';
[Wsort_a, idxsort_a]=sort(param.Wvec(param.ind_act));
pred_C_P=pred_C_C_diag(end-365*yr:end-yrend,param.ind_pass)';
[Wsort_p, idxsort_p]=sort(param.Wvec(param.ind_pass));

subplot(4,1,3)
% surface(t(end-365*yr:end-yrend),param.Wvec(param.ind_act),log10(pred_C_C_diag(end-365*yr:end-yrend,param.ind_act)'))
% shading interp
contourf(t(end-365*yr:end-yrend),Wsort_a,log10(pred_C_A(idxsort_a,:)),30,'edgecolor','none');

hold on
for i=13:12:length(tinterval)-1
plot([tinterval(i), tinterval(i)], [1e-8 1e8],'w--')
end
colorbar
set(gca,'yscale','log')
title('Log_{10} Predation from Copepods on active cops.')
set(gca,'XTick',tinterval); 
set(gca,'XtickLabels',months)
set(gca,'yTick',[1e-3, 1e-0, 1e3]); 
set(gca,'ytickLabels',{'10^{-3}','10^{0}','10^{3}'})
xlim([t(end-365*yr) t(end)-yrend])
ylim([min(param.Wvec) max(param.Wvec)])
caxis([-8 -1.5])
cbh = colorbar ; %Create Colorbar
cbh.Ticks = -8:2:-2; %Create 8 ticks from zero to 1
cbh.TickLabels = [{'10^{-8}'}, {'10^{-6}'}, {'10^{-4}'}, {'10^{-2}'}] ;
ylabel({'Body mass'; '[\mugC]'})

subplot(4,1,4)
% surface(t(end-365*yr:end-yrend),param.Wvec(param.ind_pass),log10(pred_C_C_diag(end-365*yr:end-yrend,param.ind_pass)'))
% shading interp
contourf(t(end-365*yr:end-yrend),Wsort_p,log10(pred_C_P(idxsort_p,:)),30,'edgecolor','none');
hold on
for i=13:12:length(tinterval)-1
plot([tinterval(i), tinterval(i)], [1e-8 1e8],'w--')
end
colorbar
set(gca,'yscale','log')
set(gca,'XTick',tinterval); 
set(gca,'XtickLabels',months)
title('Log_{10} Predation from copepods on passive cops.')
ylim([min(param.Wvec) max(param.Wvec)])
xlim([t(end-365*yr) t(end)-yrend])
set(gca,'yTick',[1e-3, 1e-0, 1e3]); 
set(gca,'ytickLabels',{'10^{-3}','10^{0}','10^{3}'})
caxis([-8 -1.5])
cbh = colorbar ; %Create Colorbar
cbh.Ticks = -8:2:-2; %Create 8 ticks from zero to 1
cbh.TickLabels = [{'10^{-8}'}, {'10^{-6}'}, {'10^{-4}'}, {'10^{-2}'}] ;

ylabel({'Body mass'; '[\mugC]'})
xlabel('Fraction of year')

subplot(4,1,1)
% surface(t(end-365*yr:end-yrend),param.V,pred_P_P_diag(end-365*yr:end-yrend,:)')
% shading interp
contourf(t(end-365*yr:end-yrend),param.V,pred_P_P_diag(end-365*yr:end-yrend,:)',10,'edgecolor','none')
hold on
colormap(viridis_cmap)
for i=13:12:length(tinterval)-1
plot([tinterval(i), tinterval(i)], [1e-8 1e8],'w--')
end
colorbar
set(gca,'yscale','log')
set(gca,'XTick',tinterval); 
set(gca,'XtickLabels',months)
xlim([t(end-365*yr) t(end)-yrend])
ylim([min(param.V) max(param.V)])
set(gca,'yTick',[1e-5, 1e-4, 1e-3]); 
set(gca,'ytickLabels',{'10^{-5}','10^{-4}','10^{-3}'})
title('Predation from protists on protists.')
ylabel({'Body mass'; '[\mugC]'})

plot_settings(fsize)

end

function function_plot_all(t,N,P,C,param,fsize,yr,yrend)


tinterval=t(end-365*yr):30.5:t(end)-yrend;

Gyear=P(end-365*yr:end-yrend,:)./param.V;
gg1=sum(Gyear(:,param.V<=1e-5),2);
gg2=sum(Gyear(:,param.V>1e-5 & param.V<=1e-4),2);
gg3=sum(Gyear(:,param.V>1e-4 & param.V<=1e-3),2);
gg4=sum(Gyear(:,param.V>1e-3 & param.V<=1e-2),2);
gg5=sum(Gyear(:,param.V>1e-2),2);

gg=cat(2,gg1,gg2,gg3,gg4,gg5);
months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
% bluemap=brewermap(5,'Blues');
bluemap=flip(viridis(5));
figure

subplot(3,1,1)
h =plot(t(end-365*yr:end-yrend),gg,'linewidth',1.5);
hold on
for i=13:12:length(tinterval)-1
plot([tinterval(i), tinterval(i)], [1e-5 1e8],'k--')
end

ylim([0 max(sum(gg,2))])
xlim([t(end-365*yr) t(end)-yrend])
title('Protists')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)
ylabel('Biomass [# L^{-1}]')
set(h,{'Color'},{bluemap(1,:);bluemap(2,:);bluemap(3,:);bluemap(4,:);bluemap(5,:)});
xlabel('Month')
legend({'           m \leq 10^{-5}','10^{-5} < m \leq 10^{-4}','10^{-4} < m \leq 10^{-3}','10^{-3} < m \leq 10^{-2}','10^{-2} < m'})
plot_settings(fsize)
grid on
set(gca,'yscale','log')

Cyear=C(end-365*yr:end-yrend,param.ind_act)*1000./param.Wvec(param.ind_act)';
cc1_a=sum(Cyear(:,param.W(param.ind_act)<=1e0),2);
cc2_a=sum(Cyear(:,param.W(param.ind_act)>1e0 & param.W(param.ind_act)<=1e1),2);
cc3_a=sum(Cyear(:,param.W(param.ind_act)>1e1 & param.W(param.ind_act)<=1e2),2);
cc4_a=sum(Cyear(:,param.W(param.ind_act)>1e2 & param.W(param.ind_act)<=1e3),2);
cc5_a=sum(Cyear(:,param.W(param.ind_act)>1e3),2);

Cyear_p=C(end-365*yr:end-yrend,param.ind_pass)*1000./param.W(param.ind_pass);
cc1_p=sum(Cyear_p(:,param.W(param.ind_pass)<=1e0),2);
cc2_p=sum(Cyear_p(:,param.W(param.ind_pass)>1e0 & param.W(param.ind_pass)<=1e1),2);
cc3_p=sum(Cyear_p(:,param.W(param.ind_pass)>1e1 & param.W(param.ind_pass)<=1e2),2);
cc4_p=sum(Cyear_p(:,param.W(param.ind_pass)>1e2 & param.W(param.ind_pass)<=1e3),2);
cc5_p=sum(Cyear_p(:,param.W(param.ind_pass)>1e3),2);

% Cyear=C(end-365*yr:end-yrend,param.ind_a(1:param.C_sp_act));%*1000./param.Wvec(param.ind_act)';
% cc1_a=sum(Cyear(:,param.W(param.ind_a(1:param.C_sp_act))<=1e0),2);
% cc2_a=sum(Cyear(:,param.W(param.ind_a(1:param.C_sp_act))>1e0 & param.W(param.ind_a(1:param.C_sp_act))<=1e1),2);
% cc3_a=sum(Cyear(:,param.W(param.ind_a(1:param.C_sp_act))>1e1 & param.W(param.ind_a(1:param.C_sp_act))<=1e2),2);
% cc4_a=sum(Cyear(:,param.W(param.ind_a(1:param.C_sp_act))>1e2 & param.W(param.ind_a(1:param.C_sp_act))<=1e3),2);
% cc5_a=sum(Cyear(:,param.W(param.ind_a(1:param.C_sp_act))>1e3),2);
% 
% Cyear_p=C(end-365*yr:end-yrend,param.ind_a(param.C_sp_act+1:end));%*1000./param.W(param.ind_pass);
% cc1_p=sum(Cyear_p(:,param.W(param.ind_a(param.C_sp_act+1:end))<=1e0),2);
% cc2_p=sum(Cyear_p(:,param.W(param.ind_a(param.C_sp_act+1:end))>1e0 & param.W(param.ind_a(param.C_sp_act+1:end))<=1e1),2);
% cc3_p=sum(Cyear_p(:,param.W(param.ind_a(param.C_sp_act+1:end))>1e1 & param.W(param.ind_a(param.C_sp_act+1:end))<=1e2),2);
% cc4_p=sum(Cyear_p(:,param.W(param.ind_a(param.C_sp_act+1:end))>1e2 & param.W(param.ind_a(param.C_sp_act+1:end))<=1e3),2);
% cc5_p=sum(Cyear_p(:,param.W(param.ind_a(param.C_sp_act+1:end))>1e3),2);

cc_a=cat(2,cc1_a,cc2_a,cc3_a,cc4_a,cc5_a);%./sum(Cyear,2);
cc_p=cat(2,cc1_p,cc2_p,cc3_p,cc4_p,cc5_p);%./sum(Cyear_p,2);

bluemap=flip(viridis(5));
set(0,'DefaultAxesColorOrder',bluemap);
% subplot(3,2,1)
% axes(ha(2))
subplot(3,1,2)
plot(t(end-365*yr:end-yrend),cc_a,'linewidth',1.5)
hold on
title('Active copepods')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
for i=13:12:length(tinterval)-1
plot([tinterval(i), tinterval(i)], [1e-5 1e8],'k--')
end
set(gca,'XtickLabels',months)
ylabel('Abundance [# m^{-3}]')
% set(gca,'yscale','log')
ylim([0 max(sum(cc_a,2))])
xlim([t(end-365*yr) t(end)-yrend])
grid on
set(gca,'yscale','log')

bluemap=flip(viridis(5));
set(0,'DefaultAxesColorOrder',bluemap);
% subplot(3,2,2)
subplot(3,1,3)
% axes(ha(3))
plot(t(end-365*yr:end-yrend),cc_p,'linewidth',1.5)
hold on
title('Ambush copepods')
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
for i=13:12:length(tinterval)-1
plot([tinterval(i), tinterval(i)], [1e-5 1e8],'k--')
end
set(gca,'XtickLabels',months)
ylabel('Abundance [# m^{-3}]')
% legend('m \leq 10^{-1} ','10^{-1} < m \leq 10^{0}','10^{0} < m \leq 10^{1}','10^{1} < m \leq 10^{2}','10^{2} < m')
legend('m \leq 10^{0} ','10^{0} < m \leq 10^{1}','10^{1} < m \leq 10^{2}','10^{2} < m \leq 10^{3}','10^{3} < m')
% ylim([0 5])
ylim([0 max(sum(cc_p,2))])
xlim([t(end-365*yr) t(end)-yrend])
set(gca,'yscale','log')
grid on

plot_settings(fsize)

end

function function_diagnostics_protists(param,fsize,JL_diag,JN_diag,JF_diag,mu_diag,pred_P_P_diag,pred_C_P_diag,mback_diag)

figure
subplot(1,2,1)
loglog(param.V,mean(JL_diag(ceil(end*3/4):end,:)),'g','linewidth',1.5)
hold on
loglog(param.V,mean(JN_diag(ceil(end*3/4):end,:)),'b','linewidth',1.5)
loglog(param.V,mean(JF_diag(ceil(end*3/4):end,:)),'m','linewidth',1.5)
loglog(param.V,mean(mu_diag(ceil(end*3/4):end,:)),'k--','linewidth',1.5)
loglog(param.V,param.mu_max,'k:','linewidth',1)
loglog(param.V,param.R,'r:','linewidth',1)
xlim([min(param.V) max(param.V)])
ylim([1e-3 1e1])
legend({'J_L','J_N','J_F','division rate','J_{max}','J_R'},'Location','southwest')
xlabel('Cell mass [\mugC]')
ylabel('[d^{-1}]')
title('Gains protists')
plot_settings(fsize)

subplot(1,2,2)
loglog(param.V,mean(pred_P_P_diag(ceil(end*3/4):end,:)),'g','linewidth',1.5)
hold on
loglog(param.V,mean(pred_C_P_diag(ceil(end*3/4):end,:)),'b','linewidth',1.5)
loglog(param.V,mean(mback_diag(ceil(end*3/4):end,:)),'r:','linewidth',1.5)
loglog(param.V,mean(mu_diag(ceil(end*3/4):end,:)),'k--','linewidth',1.5)
loglog(param.V,param.mu_max,'b:','linewidth',1)
loglog(param.V,mean(pred_P_P_diag(ceil(end*3/4):end,:)) + mean(mback_diag(ceil(end*3/4):end,:)) + mean(pred_C_P_diag(ceil(end*3/4):end,:)),'m--','linewidth',1)
legend('predation from protists','predation from cops','background mort','J_{max}','division rate')
xlim([min(param.V) max(param.V)])
ylim([1e-3 1e1])
xlabel('Cell mass [\mugC]')
ylabel('[d^{-1}]')
title('Losses protists')

plot_settings(fsize)



end


function function_Flvlsources(param,fsize,t,Pfrac_diag,Cfrac_diag,Dfrac_diag,Flvl_diag,yr,yrend)


months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];

Pfrac_a=Pfrac_diag(end-365*yr:end-yrend,param.ind_act)';
Cfrac_a=Cfrac_diag(end-365*yr:end-yrend,param.ind_act)';
Dfrac_a=Dfrac_diag(end-365*yr:end-yrend,param.ind_act)';
Flvl_a=Flvl_diag(end-365*yr:end-yrend,param.ind_act)';
[Wsort_a, idxsort_a]=sort(param.Wvec(param.ind_act));
Pfrac_p=Pfrac_diag(end-365*yr:end-yrend,param.ind_pass)';
Cfrac_p=Cfrac_diag(end-365*yr:end-yrend,param.ind_pass)';
Dfrac_p=Dfrac_diag(end-365*yr:end-yrend,param.ind_pass)';
Flvl_p=Flvl_diag(end-365*yr:end-yrend,param.ind_pass)';
[Wsort_p, idxsort_p]=sort(param.Wvec(param.ind_pass));

figure 
subplot(4,2,3)
% surface(t(end-365*yr:end-yrend),param.Wvec(param.ind_act),Pfrac_diag(end-365*yr:end-yrend,param.ind_act)'); 
% shading flat; 
contourf(t(end-365*yr:end-yrend),Wsort_a,Pfrac_a(idxsort_a,:),10,'edgecolor','none'); 
set(gca,'yscale','log'), colorbar
title('Fractio F_{lvl} from eating protists')
ylabel('Body mass')
xlabel('Month')
axis tight
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)

subplot(4,2,4)
% surface(t(end-365*yr:end-yrend),param.Wvec(param.ind_pass),Pfrac_diag(end-365*yr:end-yrend,param.ind_pass)'); 
contourf(t(end-365*yr:end-yrend),Wsort_p,Pfrac_p(idxsort_p,:),10,'edgecolor','none'); 
% shading flat; 
set(gca,'yscale','log'), colorbar
title('Fractio F_{lvl} from eating protists')
ylabel('Body mass')
xlabel('Month')
axis tight
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)

subplot(4,2,5)
% surface(t(end-365*yr:end-yrend),param.Wvec(param.ind_act),Cfrac_diag(end-365*yr:end-yrend,param.ind_act)'); 
% shading flat;
contourf(t(end-365*yr:end-yrend),Wsort_a,Cfrac_a(idxsort_a,:),10,'edgecolor','none'); 
set(gca,'yscale','log'), colorbar
title('Fractio F_{lvl} from eating copepods')
ylabel('Body mass')
xlabel('Month')
axis tight
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)

subplot(4,2,6)
% surface(t(end-365*yr:end-yrend),param.Wvec(param.ind_pass),Cfrac_diag(end-365*yr:end-yrend,param.ind_pass)'); 
% shading flat;
contourf(t(end-365*yr:end-yrend),Wsort_p,Cfrac_p(idxsort_p,:),10,'edgecolor','none'); 
set(gca,'yscale','log'), colorbar
title('Fractio F_{lvl} from eating copepods')
ylabel('Body mass')
xlabel('Month')
axis tight
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)

subplot(4,2,7)
% surface(t(end-365*yr:end-yrend),param.Wvec(param.ind_act),Dfrac_diag(end-365*yr:end-yrend,param.ind_act)'); 
% shading flat; 
contourf(t(end-365*yr:end-yrend),Wsort_a,Dfrac_a(idxsort_a,:),10,'edgecolor','none'); 
set(gca,'yscale','log'), colorbar
title('Fractio F_{lvl} from eating fecal pellets')
ylabel('Body mass')
xlabel('Month')
axis tight
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)

subplot(4,2,8)
% surface(t(end-365*yr:end-yrend),param.Wvec(param.ind_pass),Dfrac_diag(end-365*yr:end-yrend,param.ind_pass)'); 
% shading flat; 
contourf(t(end-365*yr:end-yrend),Wsort_p,Dfrac_p(idxsort_p,:),10,'edgecolor','none'); 
set(gca,'yscale','log'), colorbar
title('Fractio F_{lvl} from eating pellets')
ylabel('Body mass')
xlabel('Month')
axis tight
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)

subplot(4,2,1)
% surface(t(end-365*yr:end-yrend),param.Wvec(param.ind_act),Flvl_diag(end-365*yr:end-yrend,param.ind_act)'); 
% shading flat; 
contourf(t(end-365*yr:end-yrend),Wsort_a,Flvl_a(idxsort_a,:),10,'edgecolor','none'); 
set(gca,'yscale','log'), colorbar
title('Feeding level active feeders')
ylabel('Body mass')
xlabel('Month')
axis tight
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)

subplot(4,2,2)
% surface(t(end-365*yr:end-yrend),param.Wvec(param.ind_pass),Flvl_diag(end-365*yr:end-yrend,param.ind_pass)'); 
% shading flat; 
contourf(t(end-365*yr:end-yrend),Wsort_p,Flvl_p(idxsort_p,:),10,'edgecolor','none'); 
set(gca,'yscale','log'), colorbar
title('Feeding level passive feeders')
ylabel('Body mass')
xlabel('Month')
axis tight
set(gca,'XTick',t(end-365*yr):30.5:t(end)-yrend); 
set(gca,'XtickLabels',months)

colormap(viridis(10))
fsize=10;
plot_settings(fsize)
end