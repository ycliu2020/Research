%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-10-15 18:42:44
% LastEditTime : 2020-11-18 21:38:09
% LastEditors  : LYC
% Description  : plot the exact area of east china 
% FilePath     : /code/p2_processCMIP6Data/s2.radEffTrend/contribAnalysis/caseLocationShow_All.m
%  
%%---------------------------------------------------------
clc; clear;
nowpath = pwd;


%% Part 1 mask east china
% load mask map
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load world land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load world land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')
% world
latRange = 90;
lon1 = [2.5-180 357.5-180]; lat1 = [-latRange latRange]; 
% % china east
% lon1 = [65 144]; lat1 = [12 55]; 
% the Northern Hemisphere
% lon1 = [60 150]; lat1 = [0 60]; % world area
% lon1 = [60 150]; lat1 = [0 60]; % world area
set(0, 'defaultfigurecolor', 'w')
set(0, 'DefaultFigureVisible', 'on')
ss = get(0, 'ScreenSize');
h = figure('Position', [1128 238 715 600]);  %[1128 238 715 600], [-1051 -733 983 600]
% clf reset
set(h, 'Color', [1 1 1]);

m_proj('Equidistant Cylindrical', 'lon', lon1, 'lat', lat1); %Mercator,Equidistant cylindrical,lambert,Miller Cylindrical
m_coast('patch',[218, 220, 224]/360,'edgecolor','k'); 
m_grid('linestyle', 'none', 'tickdir', 'out', 'yaxislocation', 'left',  'fontsize', 16, 'fontweight','bold', 'linewi',1.5, 'color', 'k','backcolor','w'); %

% m_line(world_mapx(:), world_mapy(:), 'color', [0 0 0], 'LineWidth', 0.5);
hold on
% 行政区
shp_World=shaperead('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/map_data/中国、世界多种shp文件/世界行政区_国家_NE版.shp');
% bndryX_World=shp_World.X;
% bndryY_World=shp_World.Y;
% m_line(shp_World(:).X, shp_World(:).Y, 'color', [0 0 0], 'LineWidth', 0.5);
% hold on 
size_shp_World=size(shp_World);
for areaNum = 1:size_shp_World(1)
    m_line(shp_World(areaNum).X, shp_World(areaNum).Y, 'color', [0 0 0], 'LineWidth', 0.5);
    hold on
end
% m_line(shp_World.X(:), shp_World.Y(:), 'color', [0 0 0], 'LineWidth', 0.5);
% %填色
% % m_line(world_mapx(:), world_mapy(:), 'color', [0 0 0], 'LineWidth', 0.5);  % world boudary
% maskchina_cpNum=double(mask1);

% lonMap=lon720;latMap=lat720;
% maskchina_cpNum(lonMap<112,:)=0;
% maskchina_cpNum(:,latMap>38)=0;
% m_pcolor(lonMap+(lonMap(1)-lonMap(2))/2, latMap+(latMap(1)-latMap(2))/2, maskchina_cpNum');
% colorDIY=[1,1,1;106/180, 61/180, 154/180];
% colormap(colorDIY)
% hold on

% 外框线 
linWi=1.5;
% china east
cnEast_bndry_lon=[112 112 125 125 112];
cnEast_bndry_lat=[20 38 38 20 20];
m_line(cnEast_bndry_lon, cnEast_bndry_lat, 'linewi',linWi,'color','r')
hold on

% USA east
USAEast_bndry_lon=[-90 -90 -70 -70 -90];
USAEast_bndry_lat=[30 45 45 30 30];
m_line(USAEast_bndry_lon, USAEast_bndry_lat, 'linewi',linWi,'color','r')
hold on

% EUR west
EURwest_bndry_lon=[-9 -9 30 30 -9];
EURwest_bndry_lat=[35 60 60 35 35];
m_line(EURwest_bndry_lon, EURwest_bndry_lat, 'linewi',linWi,'color','r')
hold on

save_png('/home/liuyc/Research/P02.Ts_change_research/figure/proj2_cmip6Result/caseShow.png')%high resolution





