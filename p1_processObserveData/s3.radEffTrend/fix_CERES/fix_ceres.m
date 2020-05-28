clc; clear;tic;
inputpath = 'D:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\project\repeat_experiment\04.plot_radio_component\fix_CERES\';
load([inputpath,'trend_cldrad_ceres_year.mat'])%trendm(144,72)
inputpath_component = 'E:\Repeat_the_experiment\testdata\data_radcomponent\';
load([inputpath_component,'trend_radset_cerescld_allskysfc.mat']);% 加载变量:mon_trend_dR_q,mon_trend_dR_alb,mon_trend_dR_ts,mon_trend_dR_t,mon_trend_dR_cloud,sea_trend_dR_q,sea_trend_dR_alb,sea_trend_dR_ts,sea_trend_dR_t,sea_trend_dR_cloud,yr_trend_dR_q,yr_trend_dR_alb,yr_trend_dR_ts ,yr_trend_dR_t,yr_trend_dR_cloud

trendm_666=trendm;
input_p=[42.5 182.5 320];
n_section=length(input_p);

%计算对应的格点标签
%第一个永远是拼接而成
input_lab=input_p;
for index = 1:n_section
    input_lab(index)=find(lone==input_p(index));
end

%计算格点数
input_n=input_p;
for index = 2:n_section
    input_n(index)=(input_lab(index)-input_lab(index-1));
end
input_n(1)=144-sum(input_n(2:n_section));

%预先赋予各分布部矩阵的size

trendm=trendm';
%赋值
section{1}=[trendm(:,input_lab(3):144) trendm(:,1:input_lab(1)-1)];
section{2}=trendm(:,input_lab(1):input_lab(2)-1);
section{3}=trendm(:,input_lab(2):input_lab(3)-1);

%比较并迭代
k=0;
for jj = 1:111122
    k=k+1;
    nt1=k;
    temp1=section{nt1};
    if k==n_section
        k=0;
    end
    nt2=k+1;
    temp2=section{nt2};
    min_t1=min(size(temp1));
    min_t2=min(size(temp2));
    
    for ii = 1:72
        no_gradient=(temp1(ii,end)+temp2(ii,1)*(min_t2/min_t1))/(1+(min_t2/min_t1));
        pluse_t1=no_gradient-temp1(ii,end);
        pluse_t2=no_gradient-temp2(ii,1);
        temp1(ii,:)=temp1(ii,:)+pluse_t1;
        temp2(ii,:)=temp2(ii,:)+pluse_t2;
    end
    section{nt1}=temp1;
    section{nt2}=temp2;
    
end

%拼接回原矩阵
temp_1=section{1};
num_1=length(1:input_lab(1)-1);
trendm666(:,1:input_lab(1)-1)=temp_1(:,end-num_1+1:end);
trendm666(:,input_lab(3):144)=temp_1(:,1:end-num_1);
trendm666(:,input_lab(1):input_lab(2)-1)=section{2};
trendm666(:,input_lab(2):input_lab(3)-1)=section{3};

inputpath = 'D:\onedrive\OneDrive - smail.nju.edu.cn\a_research\P02.Ts_change_research\code\p1_processObserveData\s4.RadeffTrend\fix_CERES\';
load([inputpath,'trend_cldrad_ceres_year.mat'])
load([inputpath_Rheating,'trend_tsrad.mat']);% 加载变量: 'lon', 'lat', 'time',cons_rhs, cons_Ts, lat, lon, p_rhs, p_Ts, time, trend_rhs, trend_Ts

load('D:\onedrive\OneDrive - smail.nju.edu.cn\a_research study\xx.code set\matlab\tools\map\02.world map\mat_file\mask\mask_cp144.mat') % load word land mask
load('D:\onedrive\OneDrive - smail.nju.edu.cn\a_research study\xx.code set\matlab\tools\map\02.world map\mat_file\mask\mask_ce72.mat') % load word land mask
load('D:\onedrive\OneDrive - smail.nju.edu.cn\a_research study\xx.code set\matlab\tools\map\02.world map\mat_file\correct_worldmap.mat') % 订正后的世界地图文件word_mapx(:),word_mapy(:)
load('D:\onedrive\OneDrive - smail.nju.edu.cn\a_research study\xx.code set\matlab\tools\map\01.china map\mat_file\mask14472.mat')
% inputpath_component = 'E:\Repeat_the_experiment\testdata\data_radcomponent\';
% load([inputpath_component,'trend_radset_cerescld_allskysfc.mat']);% 加载变量:mon_trend_dR_q,mon_trend_dR_alb,mon_trend_dR_ts,mon_trend_dR_t,mon_trend_dR_cloud,sea_trend_dR_q,sea_trend_dR_alb,sea_trend_dR_ts,sea_trend_dR_t,sea_trend_dR_cloud,yr_trend_dR_q,yr_trend_dR_alb,yr_trend_dR_ts ,yr_trend_dR_t,yr_trend_dR_cloud

pmask=1;
month = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
season = {'MAM','JJA','SON','DJF'};
unite = {' trend(0.2K/10a)',' trend(W\cdotm-2/10a)',...
' trend(W\cdotm-2/10a)',' trend(W\cdotm-2/10a)',...
' trend(W\cdotm-2/10a)',' trend(W\cdotm-2/10a)'};
dataname = {'ERAi','CERES'};
component={'Ts','RHeating radiation','Ta radeffect','q radeffect','Cloud radeffect','albedo radeffect'};
component_label={'Ts','CErhs','dR_t','dR_q','dR_cloud','dR_alb'};


mmin=-5;
mmax=-mmin;
% parameter
p1=[1 2 1 1 2 1];% dataname

%画图前设置
set(0,'defaultfigurecolor','w')%设置画布底色为白色
ss = get(0,'ScreenSize');                                   % ???????????
h=figure('Position',[ss(4)/10-2000 ss(3)/20 ss(3)*3/7 ss(4)*4/6.8]);   % ???????????????

% set(gcf,'outerposition',get(0,'screensize'));%设置figure全屏
% lon1=[2.5 357.5];lat1=[-75 75];%画布经纬度范围
if pmask==1
    lon1=[2.5 357.5];lat1=[-75 75];% world area
elseif pmask==2
    lon1=[70,140];lat1=[0,60];%china area
end

hold on
[Xlon,Ylat] = meshgrid(lone,late);
for ii = 1:6
    if ii==1||ii==3||ii==5
        i=(ii-1)/2+1;j=1;
    else
        i=ii/2;j=2;
    end
    hold on
    m_proj('Mercator','lon',lon1, 'lat',lat1);%投影坐标等预设参数Mercator,Equidistant cylindrical,Miller Cylindrical
    m_pcolor(lone,late,trendm666)
    %     m_contourf(Xlon,Ylat,trendm',30,'LineStyle','none')
    colormap(mycolor(18));       %翻转颜色条colormap(flipud(mycolor(13)));%colormap(jet(4))%快速看增加和减小趋势
    caxis([mmin mmax]);
    hold on
    m_line(world_mapx(:),world_mapy(:),'color',[0 0 0],'LineWidth',0.5);

    m_grid('linestyle','none','tickdir','out','yaxislocation','left','fontsize',10,'color','k');%边框坐标轴选项
% %     xnum=[45,182.5,320];%[42.5,112.5,182.5,255,327.5];%[];%
%     ynum=[];
%     m_grid('linestyle','-','xtick',input_p,'ytick',ynum,'tickdir','in','yaxislocation','left','fontsize',10,'color','k');%边框坐标轴选项
    
    h1=title({[dataname{p1(ii)},' ',component{ii},unite{ii}];['(cc = ',num2str(cc(ii)),')']});% cc=',num2str(corr))
    hold on
    pos = get(gca,'Position');
end
colorbar('Position',[pos(1)+pos(3)+0.01 pos(2) 0.015 pos(4)*3.7])

% colorbar('Position',[pos(1)+pos(3)+0.005 pos(2) 0.01 pos(4)*1.5])
caxis([mmin mmax])
