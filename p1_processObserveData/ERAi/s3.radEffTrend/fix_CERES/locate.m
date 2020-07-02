clc; clear;tic; 
add_all;
inputpath = 'D:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\project\repeat_experiment\04.plot_radio_component\fix_CERES\';
load([inputpath,'trend_cldrad_ceres_year.mat'])
load('D:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\tools\map\a_world_map\mask\mask_cp144.mat') % load word land mask
load('D:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\tools\map\a_world_map\mask\mask_ce72.mat') % load word land mask
load('D:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\tools\map\a_world_map\correct_worldmap.mat') % 订正后的世界地图文件word_mapx(:),word_mapy(:)
load('D:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\tools\map\a_china_map\mask14472.mat')
inputpath_component = 'E:\Repeat_the_experiment\testdata\data_radcomponent\';
load([inputpath_component,'trend_radset_cerescld_allskysfc.mat']);% 加载变量:mon_trend_dR_q,mon_trend_dR_alb,mon_trend_dR_ts,mon_trend_dR_t,mon_trend_dR_cloud,sea_trend_dR_q,sea_trend_dR_alb,sea_trend_dR_ts,sea_trend_dR_t,sea_trend_dR_cloud,yr_trend_dR_q,yr_trend_dR_alb,yr_trend_dR_ts ,yr_trend_dR_t,yr_trend_dR_cloud

pmask=1;
month = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
season = {'MAM','JJA','SON','DJF'};
unite = {' trend(W·m-2/10a)'};
dataname = {'ERAi','CERES'};
component={'Ts radeffect','RHeating radiation','Cloud radeffect','Ta radeffect','q radeffect','albedo radeffect'};
component_label={'dR_ts','CErhs','dR_cloud','dR_t','dR_q','dR_alb'};


mmin=-4;
mmax=-mmin;
% parameter
p1=[1 2 2 1 1 1];% dataname 

%画图前设置
set(0,'defaultfigurecolor','w')%设置画布底色为白色
% set(gcf,'outerposition',get(0,'screensize'));%设置figure全屏
% lon1=[2.5 357.5];lat1=[-75 75];%画布经纬度范围
if pmask==1
    lon1=[2.5 357.5];lat1=[-75 75];% world area
elseif pmask==2
    lon1=[70,140];lat1=[0,60];%china area
end
figure
hold on
[Xlon,Ylat] = meshgrid(lone,late);
for ii = 3:3
    if ii==1||ii==3||ii==5
       i=(ii-1)/2+1;j=1;
    else
       i=ii/2;j=2;
    end
    hold on
    m_proj('Mercator','lon',lon1, 'lat',lat1);%投影坐标等预设参数Mercator,Equidistant cylindrical,Miller Cylindrical
    m_pcolor(lone,late,trendm')
%     m_contourf(Xlon,Ylat,trendm',30,'LineStyle','none')
    colormap(mycolor(16));       %翻转颜色条colormap(flipud(mycolor(13)));%colormap(jet(4))%快速看增加和减小趋势
    caxis([mmin mmax]);
    hold on
    m_line(world_mapx(:),world_mapy(:),'color',[0 0 0],'LineWidth',0.5);

    xnum=[45,182.5,320];%[42.5,112.5,182.5,255,327.5];%[];%
    ynum=[];
    m_grid('linestyle','-','xtick',xnum,'ytick',ynum,'tickdir','in','yaxislocation','left','fontsize',10,'color','k');%边框坐标轴选项

    h1=title({[dataname{p1(ii)},' ',component{ii},unite{1}];['year mean']});% cc=',num2str(corr))
    hold on
    pos = get(gca,'Position');
end
colorbar
% colorbar('Position',[pos(1)+pos(3)+0.005 pos(2) 0.01 pos(4)*1.5])
caxis([mmin mmax])






