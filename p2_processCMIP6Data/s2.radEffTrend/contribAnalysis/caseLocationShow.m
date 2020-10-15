%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-10-15 18:42:44
% LastEditTime : 2020-10-15 19:37:27
% LastEditors  : LYC
% Description  : plot the exact area of east china 
% FilePath     : /code/p2_processCMIP6Data/s2.radEffTrend/contribAnalysis/caseLocationShow.m
%  
%%---------------------------------------------------------
clc; clear;
nowpath = pwd;


%% Part 1 mask east china
% load mask map
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/中国干湿以及青藏高原分区地图/main_use/mask720.mat')
% 国境线
bou_china=shaperead('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/map_data/中国及各省shp/国界线/bou1_4p.shp'); % load china  boundary
bou_chinaX=[bou_china(:).X];bou_chinaY=[bou_china(:).Y];
% 国境线(带南海线)
bou_china_line=shaperead('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/map_data/中国及各省shp/国界线/bou1_4l.shp'); % load china  boundary
bou_china_lineX=[bou_china_line(:).X];bou_china_lineY=[bou_china_line(:).Y];
% 省界线
bou_chinaProvince=shaperead('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/map_data/中国及各省shp/省界线/bou2_4p.shp'); % load china  boundary
bou_chinaProvinceX=[bou_chinaProvince(:).X];bou_chinaProvinceY=[bou_chinaProvince(:).Y];

%%  plot test
lon1 = [65 144]; lat1 = [12 55]; % world area
% lon1 = [60 150]; lat1 = [0 60]; % world area
set(0, 'defaultfigurecolor', 'w')
set(0, 'DefaultFigureVisible', 'off')
ss = get(0, 'ScreenSize');
h = figure('Position', [1128 238 715 600]);  %[1128 238 715 600], [-1051 -733 983 600]
% clf reset
set(h, 'Color', [1 1 1]);

m_proj('Mercator', 'lon', lon1, 'lat', lat1); %Mercator,Equidistant cylindrical,lambert,Miller Cylindrical

% m_line(world_mapx(:), world_mapy(:), 'color', [0 0 0], 'LineWidth', 0.5);  % world boudary
maskchina_cpNum=double(mask1);

lonMap=lon720;latMap=lat720;
maskchina_cpNum(lonMap<112,:)=0;
maskchina_cpNum(:,latMap>38)=0;
m_pcolor(lonMap+(lonMap(1)-lonMap(2))/2, latMap+(latMap(1)-latMap(2))/2, maskchina_cpNum');
colorDIY=[1,1,1;106/180, 61/180, 154/180];
colormap(colorDIY)
hold on

% 国界线和省界线
m_line(bou_chinaX,bou_chinaY,'color','k','LineWidth',0.5);
hold on
m_line(bou_chinaProvinceX,bou_chinaProvinceY,'color','k','LineWidth',0.5);
hold on
m_line(bou_china_lineX,bou_china_lineY,'color','k','LineWidth',0.5);
hold on
% 外框线 
m_plot([112,112],[20,38],'linestyle','-','color','r','LineWidth',1.5);
hold on
m_plot([125,125],[20,38],'linestyle','-','color','r','LineWidth',1.5);
hold on
m_plot([112, 125],[38,38],'linestyle','-','color','r','LineWidth',1.5);
hold on
m_plot([112, 125],[20,20],'linestyle','-','color','r','LineWidth',1.5);
hold on
m_grid('linestyle', 'none', 'tickdir', 'out', 'yaxislocation', 'left', 'fontsize', 10, 'color', 'k'); %


h1=axes('Position',[0.8 0.2 0.09 0.2]);%创建坐标系时返回它的句柄[left bottom width height]
axes(h1);%将坐标系h1置为当前坐标系
set(gcf,'PaperPositionMode','auto')
m_proj('miller','lon',[106,122],'lat',[2,26]) %设置南海区域
m_plot(bou_china_lineX,bou_china_lineY,'color','k','LineWidth',0.5)%绘国界
m_grid('XTick',[],'YTick',[])%添加坐标

