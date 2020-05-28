clc
clear
filepath = 'C:\Yourfilepath\'; 
filepath1 = 'C:\Yourfilepath\';
filename = strcat(filepath,'Trendrad_Chinasupp.nc');
filename1 = strcat(filepath1,'TrendChina_ERAits.nc');
trend0 = ncread(filename,'trend_Rxs');
spectrum = 2;
trend = squeeze(trend0(5,1,spectrum,2,:,:,:,1));                          % x,sky,spectrum,budget,lon,lat,month.coef
trend_tot0 = ncread(filename1,'trend_ts'); 
trend_tot = squeeze(trend_tot0(:,:,:,1));  
cc = zeros(12,1); pp = zeros(12,1);
for month = 1:12
    [cc0,pp0] = corrcoef(squeeze(trend(:,:,month)),squeeze(trend_tot(:,:,month)),'Rows','complete');
    cc(month,1) = -cc0(1,2);
    pp(month,1) = -pp0(1,2);
end
% trend0 = ncread(filename,'trend_ts'); 
trend = trend*365*10;
lon = ncread(filename,'longitude'); nlon = length(lon);
lat = ncread(filename,'latitude'); nlat = length(lat);
[Xlon,Ylat] = meshgrid(lon,lat);
map3 = 'C:\Yourfilepath\maps\bou2_4l.shp';
readmap3 = shaperead(map3);
bou1_4lx=[readmap3(:).X];%提取经度信息
bou1_4ly=[readmap3(:).Y];%提取纬度信息

MAX = max(max(max(abs(trend)))); 
% if MAX > roundn(MAX,1)
%    mmin = -roundn(MAX,1)-1; mmax = roundn(MAX,1)+1;
% else
%    mmin = -roundn(MAX,1); mmax = roundn(MAX,1);
% end
mmin = -roundn(MAX,0); mmax = roundn(MAX,0);
month = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
%% Plot
figure(1)
hold all
load ('C:\Yourfilepath\mylib\color_bwr','mycmap')
for ii = 1:12
   if ii<11
       mm = ii+2;
   else
       mm = ii-10;
   end
   if ii<4
      i = ii;
      j = 1;
   elseif ii<7
      i = ii-3;
      j = 2;
   elseif ii<10
      i = ii - 6;
      j = 3;
   else
      i = ii - 9;
      j = 4;
   end
   trendm = squeeze(trend(:,:,mm));
   corr = roundn(cc(mm,1),-2);
   subplot_yih(3,4,i,j);
   hold all
   m_proj('miller','lon',[72,137],'lat',[16,55])%选择投影方式
   m_contourf(Xlon,Ylat,trendm',30,'LineStyle','none')
   m_plot(bou1_4lx,bou1_4ly)%绘国界
   colormap(gca,mycmap)
   caxis([mmin mmax]);
   m_grid('linestyle','none','tickdir','out','fontsize',8)%添加坐标
   title(strcat(month(mm),'  cc=',num2str(corr)))
%    colorbar('southoutside')   
   pos = get(gca,'Position');
   scale = 4.5;
   h1=axes('Position',[pos(1) pos(2)+0.01 pos(3)/scale pos(4)/scale]);%创建坐标系时返回它的句柄
   axes(h1);%将坐标系h1置为当前坐标系
   m_proj('miller','lon',[106,122],'lat',[2,26]) %设置南海区域
   m_contourf(Xlon,Ylat,trendm','LineStyle','none')
   colormap(gca,mycmap)
   caxis([mmin mmax]);
   hold all
   m_plot(bou1_4lx,bou1_4ly)%绘国界
   m_grid('XTick',[],'YTick',[])%添加坐标
   hold off
end
colorbar('Position',[pos(1)+pos(3)+0.005 pos(2) 0.01 pos(4)*3.73])
caxis([mmin mmax])
% saveas(gcf,'try2.tiff')
% hold off