clc; clear;tic; 
add_all;
level='sky';
inputpath_Rheating = 'E:\Repeat_the_experiment\testdata\data_Tsradtrend\';
inputpath_component = 'E:\Repeat_the_experiment\testdata\data_radcomponent\';
load([inputpath_Rheating,'trend_CERESrad_HadTs.mat']);% 加载变量: 'lon', 'lat', 'time',cons_rhs, cons_Ts, lat, lon, p_rhs, p_Ts, time, trend_rhs, trend_Ts
load([inputpath_component,'trend_radset_cerescld_allsky',level,'.mat']);% 加载变量:mon_trend_dR_q,mon_trend_dR_alb,mon_trend_dR_ts,mon_trend_dR_t,mon_trend_dR_cloud,sea_trend_dR_q,sea_trend_dR_alb,sea_trend_dR_ts,sea_trend_dR_t,sea_trend_dR_cloud,yr_trend_dR_q,yr_trend_dR_alb,yr_trend_dR_ts ,yr_trend_dR_t,yr_trend_dR_cloud
% load([inputpath_component,'trendradset_erai_allskysfc.mat']);% 加载变量:mon_trend_dR_q,mon_trend_dR_alb,mon_trend_dR_ts,mon_trend_dR_t,mon_trend_dR_cloud,sea_trend_dR_q,sea_trend_dR_alb,sea_trend_dR_ts,sea_trend_dR_t,sea_trend_dR_cloud,yr_trend_dR_q,yr_trend_dR_alb,yr_trend_dR_ts ,yr_trend_dR_t,yr_trend_dR_cloud
% load([inputpath_component,'trendradset_erai_allskysfc.mat']);% 加载变量:mon_trend_dR_q,mon_trend_dR_alb,mon_trend_dR_ts,mon_trend_dR_t,mon_trend_dR_cloud,sea_trend_dR_q,sea_trend_dR_alb,sea_trend_dR_ts,sea_trend_dR_t,sea_trend_dR_cloud,yr_trend_dR_q,yr_trend_dR_alb,yr_trend_dR_ts ,yr_trend_dR_t,yr_trend_dR_cloud
% load([inputpath_component,'trendradset_erai_allskysfc.mat']);% 加载变量:mon_trend_dR_q,mon_trend_dR_alb,mon_trend_dR_ts,mon_trend_dR_t,mon_trend_dR_cloud,sea_trend_dR_q,sea_trend_dR_alb,sea_trend_dR_ts,sea_trend_dR_t,sea_trend_dR_cloud,yr_trend_dR_q,yr_trend_dR_alb,yr_trend_dR_ts ,yr_trend_dR_t,yr_trend_dR_cloud

load('D:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\tools\map\a_world_map\mask\mask_cp144.mat') % load word land mask
load('D:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\tools\map\a_world_map\mask\mask_ce72.mat') % load word land mask
load('D:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\tools\map\a_world_map\correct_worldmap.mat') % 订正后的世界地图文件word_mapx(:),word_mapy(:)
load('D:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\tools\map\a_china_map\mask14472.mat')


month = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
season = {'MAM','JJA','SON','DJF'};
unite = {' trend(W·m-2/10a)'};
dataname = {'ERAi','CERES'};
component={'Ts radeffect','RHeating radiation','Cloud radeffect','Ta radeffect','q radeffect','albedo radeffect'};
component_label={'dR_ts','CErhs','dR_cloud','dR_t','dR_q','dR_alb'};

% multiply 10a 
% month 
trend_var_m = zeros(144,72,12,6);
for i=1:6
  temp=eval(['mon_trend_',component_label{i}]);
  trend_var_m(:,:,:,i) = squeeze(temp(:,:,:,1))*365*10;
end
% mask world land or china land(pmask=1 or 2)
pmask=1;
if pmask==1
    mask1=maskworld_cp;
elseif pmask==2
    mask1=maskchina_cp;
end
% mask_temp=repmat(mask1,[1 1 12 6]);

% trend_var_m(~mask_temp)=NaN;
% trend_var_m(:,late>60|late<-60,:,:)=NaN;%划定纬度剔除两个极地


%season(trend of season(im==1 is MAM, im==2 is JJA, im==3 is SON, im==4 is DJF))
trend_var_s = zeros(144,72,4,6);
for i=1:6
    temp=eval(['sea_trend_',component_label{i}]);
    trend_var_s(:,:,:,i) = squeeze(temp(:,:,:))*365*10;
end
% mask world land
% mask_temp=repmat(mask1,[1 1 4 6]);
% trend_var_s(~mask_temp)=NaN;
% trend_var_s(:,late>60|late<-60,:,:)=NaN;%划定纬度剔除两个极地

%year
trend_var_y = zeros(144,72,6);
for i=1:6
    temp=eval(['yr_trend_',component_label{i}]);
    trend_var_y(:,:,i) = squeeze(temp(:,:))*365*10;
end
% mask world land 
% mask_temp=repmat(mask1,[1 1 6]);
% trend_var_y(~mask_temp)=NaN;
% trend_var_y(:,late>60|late<-60,:,:)=NaN;%划定纬度剔除两个极地

%----------------------section1:year图----------------------------------
% cal the relative correlation component_label={'dR_ts','rhs','dR_cloud','dR_t','dR_q','dR_alb'};
cc = zeros(6,1); pp = zeros(6,1);
for index = 1:6
    [cc0,pp0] = corrcoef(trend_var_y(:,:,1),trend_var_y(:,:,index),'Rows','complete');
    cc(index,1)=cc0(1,2);
    pp(index,1) = pp0(1,2);
end
% pos = zeros(4,2);
%----------------------plot figure (only 60S-60N land)----------------------------------
% MAX = squeeze(max(max(max(abs(trend_var_y)))));
% mmin = -roundn(MAX,0); mmax = roundn(MAX,0);
% mmin=mmin';mmax=mmax';
mmin=-5;
mmax=-mmin;
% parameter
p1=[1 2 2 1 1 1];% dataname 

%画图前设置
set(0,'defaultfigurecolor','w')%设置画布底色为白色
% set(gcf,'outerposition',get(0,'screensize'));%设置figure全屏
% lon1=[2.5 357.5];lat1=[-75 75];%画布经纬度范围
if pmask==1
    lon1=[2.5 357.5];lat1=[-60 60];% world area
elseif pmask==2
    lon1=[70,140];lat1=[0,60];%china area
end
figure
hold on
[Xlon,Ylat] = meshgrid(lone,late);
for ii = 1:6
    trendm=squeeze(trend_var_y(:,:,ii));
    if ii==1||ii==3||ii==5
       i=(ii-1)/2+1;j=1;
    else
       i=ii/2;j=2;
    end
    subplot_yc(3,2,i,j);% 画子图
    hold on
    m_proj('Mercator','lon',lon1, 'lat',lat1);%投影坐标等预设参数Mercator,Equidistant cylindrical,Miller Cylindrical
    m_contourf(Xlon,Ylat,trendm',30,'LineStyle','none')
    colormap(mycolor(16));       %翻转颜色条colormap(flipud(mycolor(13)));%colormap(jet(4))%快速看增加和减小趋势
    caxis([mmin mmax]);
    hold on
    m_line(world_mapx(:),world_mapy(:),'color',[0 0 0],'LineWidth',0.5);
    m_grid('linestyle','none','tickdir','out','yaxislocation','left','fontsize',10,'color','k');%边框坐标轴选项

    % title({[dataname{p1},' ',component{ii}];season{i}});% cc=',num2str(corr))
%     if ii==1
%         title({[dataname{p1},' ',component{ii},unite];'year mean '});% cc=',num2str(corr))
%     else
%         title({[dataname{p1},' ',component{ii},unite];['year mean (cc = ',num2str(cc(ii)),')']});% cc=',num2str(corr))
%     end
title({[dataname{p1(ii)},' ',component{ii},unite{1}];['year mean (cc = ',num2str(cc(ii)),')']});% cc=',num2str(corr))
    hold on
    pos = get(gca,'Position');
end
colorbar('Position',[pos(1)+pos(3)+0.005 pos(2) 0.01 pos(4)*1.5])
caxis([mmin mmax])
sgtitle(upper(level),'Fontsize',14,'position',[-3.5,1.4])





