load('/data1/liuyincheng/cmip6-process/amip_2000-2014/BCC-CSM2-MR/rawdata_regrid/alb.mat')
load('/data1/liuyincheng/cmip6-process/amip_2000-2014/BCC-CSM2-MR/rawdata_regrid/global_vars.mat')

plotVar=squeeze(alb(:,:,1)); % 2000 year
plotLat=latf;
plotLon=lonf;

load('/home/lyc/lib/tools/matlab/map/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/lyc/lib/tools/matlab/map/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/lyc/lib/tools/matlab/map/02.world_map/mat_file/correct_worldmap.mat')% ????????????????word_mapx(:),word_mapy(:)
load('/home/lyc/lib/tools/matlab/map/01.china_map/mat_file/mask14472.mat')


% figure set
set(0,  'defaultfigurecolor',  'w')
areaNum = 1; % world land
p_3 = 90; % Latitude range
lon1 = [2.5 357.5]; lat1 = [-p_3 + 1 p_3 - 1]; % world area
% lon1=[70,140];lat1=[0,60]; % China area

colorBarMin = min(min(plotVar)); 
colorBarMax = max(max(plotVar));

ss = get(0, 'ScreenSize'); % ???????????
h = figure('Position', [ss(4)/2 ss(3)/35 ss(3)*3/9.5 ss(4)*4/5]); % ???????????????
set(h, 'Color', [1 1 1]);
f_matrix = reshape(1:1, [1, 1])';

% figure
subplot_yc(1, 1, 1, 1); % ?????
hold on
m_proj('miller',  'lon', lon1,  'lat', lat1); %?????????????Mercator,Equidistant cylindrical,lambert,Miller Cylindrical
m_pcolor(plotLon, plotLat, plotVar');
colormap(mycolor(18)); %????????colormap(flipud(mycolor(13)));%colormap(jet(4))%????????????????????
caxis([colorBarMin colorBarMax]);
hold on
m_line(world_mapx(:), world_mapy(:),  'color', [0 0 0],  'LineWidth', 0.5);
m_grid('linestyle',  'none',  'tickdir',  'out',  'yaxislocation',  'left',  'fontsize', 10,  'color',  'k'); %????????????
hold on
pos = get(gca,  'Position');
c = colorbar('southoutside');
c.TickLength = 0.02;
% c.LineWidth = 2;
c.Limits = [colorBarMin colorBarMax];






