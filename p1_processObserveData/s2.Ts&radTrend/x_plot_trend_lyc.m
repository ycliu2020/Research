% x_plot_month_trend_lyc.m
% plot Ts/rhs trend;trend1:(12,2) k and b, from 1-12 month
% note: units:/10a eg: trend*365*10
% Ts: K/10 a , rhs:W/m^(-2)/10 a
% color map: 16
clc; clear;tic; %计时
add_all;
pre_load;
inputpath = 'E:\Repeat_the_experiment\testdata\ERAi\';
load([inputpath,'trend_tsrad.mat']);% 加载变量: 'lon', 'lat', 'time',cons_rhs, cons_Ts, lat, lon, p_rhs, p_Ts, time, trend_rhs, trend_Ts
load('D:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\tools\绘图\map\a_world_map\mask144.mat') % load word land mask
load('D:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\tools\绘图\map\a_world_map\correct_worldmap.mat') % 订正后的世界地图文件word_mapx(:),word_mapy(:)
[Xlon,Ylat] = meshgrid(lon,lat);
month = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
season = {'MAM','JJA','SON','DJF'};
action = {'Ts trend(K/10a)','LW↓+SW trend(W·m-2/10a)'};

% multiply 10a
trend_Ts_m = squeeze(mon_trend_Ts(:,:,:,1))*365*10;
trend_rhs_m = squeeze(mon_trend_rhs(:,:,:,1))*365*10;
%season
trend_Ts_s = sea_trend_Ts(:,:,:)*365*10;
trend_rhs_s = sea_trend_rhs(:,:,:)*365*10;
%year
trend_Ts_y = yr_trend_Ts(:,:)*365*10;
trend_rhs_y = yr_trend_rhs(:,:)*365*10;

%----------------------this is every season and both equation sides(only land)----------------------------------
% trend = trend_Ts_m; % decide varible to plot
trend1 = trend_Ts_s; % decide varible to plot
trend2 = trend_rhs_s; % decide varible to plot
for ii=1:4
temp1=squeeze(trend1(:,:,ii));
temp2=squeeze(trend2(:,:,ii));
temp1(~maskworld)=NaN;%mask land
temp2(~maskworld)=NaN;%mask land
trend1(:,:,ii)=temp1;
trend2(:,:,ii)=temp2;
end
trend1(:,lat>60&lat<-60,:)=NaN;%划定纬度剔除两个极地
trend2(:,lat>60&lat<-60,:)=NaN;
cc = zeros(4,1); pp = zeros(4,1);
% Compute correlation(seasonly)
for ii = 1:4
    [cc0,pp0] = corrcoef(squeeze(trend1(:,:,ii)),squeeze(trend2(:,:,ii)),'Rows','complete');
    cc(ii,1) = -cc0(1,2);
    pp(ii,1) = -pp0(1,2);
end

lon1=[2.5 357.5];lat1=[-60 60];%画布经纬度范围


%----------------------plot figure about every season and both equation sides(only land)----------------------------------
% MAX = squeeze(max(max(max(abs(trend)))));
% mmin = -roundn(MAX,0); mmax = roundn(MAX,0);
% mmin=mmin';mmax=mmax';
mmin=-1;
mmax=1;
%画图前设置
% set(0,'defaultfigurecolor','w')%设置画布底色为白色
% set(gcf,'outerposition',get(0,'screensize'));%设置figure全屏
figure(1)
hold on
% 画4*2的season分布子图, 行为不同seasons,列为方程两项, i为行, j为列
for ii = 1:8
    if ii<5
        i=ii;j=1;kk=1;
        trendm = squeeze(trend1(:,:,ii));   
        corr = roundn(cc(ii,1),-2);% 相关系数   
    elseif ii>=5
        i=ii-4;j=2;kk=2;
        trendm = squeeze(trend2(:,:,ii-4));
    end

    subplot_yih(4,2,i,j);% 画子图
    hold on
    
    m_proj('Mercator','lon',lon1, 'lat',lat1);%投影坐标等预设参数Mercator,Equidistant cylindrical,Miller Cylindrical
    m_contourf(Xlon,Ylat,trendm',30,'LineStyle','none')
    colormap(mycolor(16));       %翻转颜色条colormap(flipud(mycolor(13)));%colormap(jet(4))%快速看增加和减小趋势
    caxis([mmin mmax]);
    hold on
    m_line(world_mapx(:),world_mapy(:),'color',[0 0 0],'LineWidth',0.5);
    
    
    m_grid('linestyle','none','tickdir','out','yaxislocation','left','fontsize',10,'color','k');%边框坐标轴选项
    if ii<5
        title(strcat(season(ii),'  cc=',num2str(corr),action(kk)));% cc=',num2str(corr))
    elseif ii>=5
        title(strcat(season(ii-4),action(kk)));% cc=',num2str(corr))
    end
    hold on
    pos = get(gca,'Position');
    
end
colorbar('Position',[pos(1)+pos(3)+0.005 pos(2) 0.01 pos(4)*3.73])
caxis([mmin mmax])


% %----------------------plot every month figure----------------------------------
% % MAX = squeeze(max(max(max(abs(trend)))));
% % mmin = -roundn(MAX,0); mmax = roundn(MAX,0);
% % mmin=mmin';mmax=mmax';
% mmin=-10;
% mmax=10;
% %画图前设置
% set(0,'defaultfigurecolor','w')%设置画布底色为白色
% set(gcf,'outerposition',get(0,'screensize'));%设置figure全屏
% figure(1)
% hold on
% % 1月开始, 计算每个月所在子图位置, 画4*3的月分布子图, 行为连续月,列为不同季节, i为行, j为列
% for ii = 1:12
%     if ii<11
%         mm = ii+2;
%     else
%         mm = ii-10;
%     end
%     if ii<4
%        i = 1;
%        j = ii;
%     elseif ii<7
%        i = 2;
%        j = ii-3;
%     elseif ii<10
%        i = 3;
%        j = ii - 6;
%     else
%        i = 4;
%        j = ii - 9;
%     end
%     trendm = squeeze(trend(:,:,mm));
%     subplot_yih(4,3,i,j);% 画子图
%     hold on
%     m_proj('Mercator','lon',lon1, 'lat',lat1);%投影坐标等预设参数Mercator,Equidistant cylindrical,Miller Cylindrical
%     m_contourf(Xlon,Ylat,trendm',30,'LineStyle','none')
%     colormap(mycolor(16));       %翻转颜色条colormap(flipud(mycolor(13)));%colormap(jet(4))%快速看增加和减小趋势
%     caxis([mmin mmax]);
%     hold on
%     m_line(word_mapx(:),word_mapy(:),'color',[0 0 0],'LineWidth',0.5);


%     m_grid('linestyle','none','tickdir','out','yaxislocation','left','fontsize',10,'color','k');%边框坐标轴选项
%     title(strcat(month(mm)));% cc=',num2str(corr))
%     hold on
%     pos = get(gca,'Position');

% end
% colorbar('Position',[pos(1)+pos(3)+0.005 pos(2) 0.01 pos(4)*3.73])
% caxis([mmin mmax])



% ------output------
% name=['Ta_',t_name{season+1},'_',plv_name{level},'_',num2str(year(1)),'-',num2str(year(2))];
% fig_name='C:\Users\LYC\Desktop\当前任务接口\222';
% export_fig( gcf , '-png','-r500'  , fig_name );

%----------------------plot figure about every season and both equation sides ----------------------------------

t=toc;
disp(t)