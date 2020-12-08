%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-10-15 18:42:44
% LastEditTime : 2020-12-08 15:00:35
% LastEditors  : Please set LastEditors
% Description  : plot the exact area of east china 
% FilePath     : /code/p3_paperFigIntegrate/Fig8_radContribAnalysis_caseLocShow_All.m
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
lon1 = [-150 150]; lat1 = [-75 75]; 
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
linWi=1.25;
% china east
lonCN=[112 122];
latCN=[22 37];
cnEast_bndry_lon=[lonCN(1) lonCN(1) lonCN(2) lonCN(2) lonCN(1)];
cnEast_bndry_lat=[latCN(1) latCN(2) latCN(2) latCN(1) latCN(1)];
m_line(cnEast_bndry_lon, cnEast_bndry_lat, 'linewi',linWi,'color','r')
hold on

% USA east
lonUSA=[-90 -80];
latUSA=[30 45];
USAEast_bndry_lon=[lonUSA(1) lonUSA(1) lonUSA(2) lonUSA(2) lonUSA(1)];
USAEast_bndry_lat=[latUSA(1) latUSA(2) latUSA(2) latUSA(1) latUSA(1)];
m_line(USAEast_bndry_lon, USAEast_bndry_lat, 'linewi',linWi,'color','r')
hold on

% EUR west
lonEUR=[-1 14];
latEUR=[44 54];
EURwest_bndry_lon=[lonEUR(1) lonEUR(1) lonEUR(2) lonEUR(2) lonEUR(1)];
EURwest_bndry_lat=[latEUR(1) latEUR(2) latEUR(2) latEUR(1) latEUR(1)];
m_line(EURwest_bndry_lon, EURwest_bndry_lat, 'linewi',linWi,'color','r')
hold on

save_png('/home/liuyc/Research/P02.Ts_change_research/figure/proj3_PaperFig/v0.0/Fig8_caseShow.png')%high resolution





