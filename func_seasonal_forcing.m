function [temp,Diff,zmld,Ls,dzdt]=func_seasonal_forcing(param, years, seasonal_switch,bif_switch,Diff_bif)




% zmld1=load('mldfile_lat55.txt')./6;
fash=load('mld_fasham_mat.mat');
zmld1=fash.alk3;
zmld1(end)=zmld1(1);
% zmld1(:)=30;
% diff1=load('difffile_lat55.txt');
% diff=load('MLD_vec_L60.mat');
temp1=load('tempfile_lat55.txt')+3;
% temp1=load('Temp_vec_L55.mat');
% temp1=temp1.vq'+273;
clouds=1-(5/8);%2.5



zmld=[zmld1; repmat(zmld1(2:end),years-1,1)]; %mixed layer depth
dzdt=zmld(2:end)-zmld(1:end-1); %change in mixed layer depth
dzdt=[dzdt(1); dzdt]; %correct for vector size

% Diff=[diff; repmat(diff(2:end),years-1,1)] + max(0,dzdt)./zmld; %mixed layer depth

Diff=(param.Diff_min+max(0,dzdt))./zmld;
% Diff=diff1(1:end1);
% Diff=(param.Diff_min)./zmld;
% Diff=(0.1+max(0,dzdt))./zmld;

kyear=0;
day_type=1;
tsol=0:1:365*years;
Lss=daily_insolation(kyear,param.lat,tsol,day_type); %surface light
Lss=daily_insolation(kyear,55,tsol,day_type); %surface light
% Ls=Lss;
Ls=Lss.*clouds.*0.4.*4.6; %correct for cloud cover, convert to PAR (PAR/I=0.4) and convert to microeinsteins (1Wm^-2=4.6microeinteins)
%temperature witihin the mixed layer
temp=[temp1(1); repmat(temp1(2:end),years,1)];
temp=temp-273;


if seasonal_switch==0
    temp(:)=15;
%     Diff(:)=0.001;
    zmld(:)=40;
    Ls(:)=50;
    dzdt(:)=0;
%     
% elseif bif_switch==1
%     temp(:)=15;
%     Diff(:)=Diff_bif;
%     zmld(:)=25;
%     Ls(:)=150;
%     dzdt(:)=0;
%     
else

% months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
% 
% figure; 
% 
% subplot(4,1,1)
% plot(Ls(1:365), 'k','linewidth',1.5)
% set(gca,'XTick',1:30.5:365); set(gca,'XtickLabels',months)
% xlim([1 365])
% title(sprintf('Surface PAR Lat %d', param.lat))
% ylabel('[\muE s^{-1} m^{-2}]')
% 
% subplot(4,1,2)
% plot(-zmld(1:365), 'k','linewidth',1.5)
% set(gca,'XTick',1:30.5:365); set(gca,'XtickLabels',months)
% xlim([1 365])
% ylim([-250 0])
% title('Mixed layer depth')
% ylabel('[m]')
% 
% cp= [0.5882    0.5882    0.5882];
% subplot(4,1,3)
% yyaxis left
% plot(max(0,dzdt), 'k--','linewidth',1.5)
% hold on
% plot(param.Diff_min.*ones(size(1:365)), 'k:','linewidth',1.5)
% plot(param.Diff_min+max(0,dzdt), 'k','linewidth',1.5)
% ylabel('[m d^{-1}]')
% ylim([0 4])
% yyaxis right
% plot(Diff(1:365),'color',cp,'linewidth',1)
% set(gca,'XTick',1:30.5:365); set(gca,'XtickLabels',months)
% xlim([1 365])
% title('Input of nutrients')
% ylabel('[d^{-1}]')
% set(gca,'Yscale','log')
% ylim([1e-3 1e0])
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.YAxis(2).Color = cp;
% columnlegend(2,{'(max(0,dzdt))','\omega','(max(0,dzdt))+\omega','(max(0,dzdt))+\omega/z_{mld}'})
% % legend('(max(0,dzdt))','\omega','(max(0,dzdt))+\omega','(max(0,dzdt))+\omega/z_{mld}')
% legend boxoff
% 
% subplot(4,1,4)
% plot(temp(1:365), 'k','linewidth',1.5)
% set(gca,'XTick',1:30.5:365); set(gca,'XtickLabels',months)
% xlim([1 365])
% title('Temperature')
% ylabel('C')
% 
% fsize=10;
% plot_settings(fsize)

end

end