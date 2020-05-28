% ECMs homework
% 形势场分析
clear; clc; tic;
% map file
load('/home/lyc/lib/tools/matlab/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/lyc/lib/tools/matlab/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/lyc/lib/tools/matlab/myMap/02.world_map/mat_file/correct_worldmap.mat')% ????????????????word_mapx(:),word_mapy(:)
load('/home/lyc/lib/tools/matlab/myMap/01.china_map/mat_file/mask14472.mat')

nowpath = pwd;
inputPath = '/data1/liuyincheng/z_homework/201907.nc';
title_level={'500','850','1000'};
title_time={'14:00','20:00'};

% read(480,241,3,124)
lat = ncread(inputPath, 'latitude');
lon = ncread(inputPath, 'longitude');
level = ncread(inputPath, 'level');
U = ncread(inputPath, 'u');
V = ncread(inputPath, 'v');
H = ncread(inputPath, 'z');
T = ncread(inputPath, 't');

set(0, 'DefaultFigureVisible', 'on')
set(0, 'defaultfigurecolor', 'w')
lon1 = [100 140]; lat1 = [20 50]; % world area
% lon1 = [0 180]; lat1 = [20 50]; % world area

%regrid 144x72(unite grids)
latw = 88.75:-2.5:-88.75; nlatw = length(lat);
lonw = 2.5:2.5:357.5; nlonw = length(lon);
nlon = length(lon); nlat = length(lat);

T_range{1}=[-14 -12 -10 -8 -6 -4 -2 0 2 4 6 8 10 12 14 16];
T_range{2}=[12 15 18 21 24 27];
T_range{3}=[3 9 12 15 18 24 27 30 33 36 40];
for varLevel = 1:1% 1.500, 2.850, 3.1000

    for varTime = 1:1% 1.14:00, 2.20:00
        timeNum = varTime + 9;
        U_temp = squeeze(squeeze(U(:, :, varLevel, timeNum)));
        V_temp = squeeze(squeeze(V(:, :, varLevel, timeNum)));
        H_temp = squeeze(squeeze(H(:, :, varLevel, timeNum)));
        T_temp = squeeze(squeeze(T(:, :, varLevel, timeNum)));
        H_temp=H_temp./10;
        T_temp=T_temp-273;
        [Xlon, Ylat] = meshgrid(lat, lon);
        [Xlonw, Ylatw] = meshgrid(latw, lonw);
        U_temp = interp2(Xlon, Ylat, U_temp, Xlonw, Ylatw);
        V_temp = interp2(Xlon, Ylat, V_temp, Xlonw, Ylatw);
        
        
        figure
        m_proj('Mercator', 'lon', lon1, 'lat', lat1); %?????????????Mercator,Equidistant cylindrical,lambert,Miller Cylindrical
        m_coast('patch',[.7 .7 .7],'edgecolor','none');
        hold on;
        H_range=5000:40:6400 ;

        [LN,LT]=meshgrid(lon,lat);
        [CS, ch] = m_contour(LN, LT, T_temp',T_range{varLevel},'color','r'); %设定数据的范围,'edgecolor','none'
        clabel(CS,ch,'fontsize',8,'color','r');
        hold on;

        [CS1, ch1] = m_contour(LN, LT, H_temp','color','b'); %设定数据的范围
        clabel(CS1,ch1,'fontsize',8,'color','b');
        hold on;
        U_temp=U_temp(lonw<=140&lonw>=100,latw<=50&latw>=20);
        V_temp=V_temp(lonw<=140&lonw>=100,latw<=50&latw>=20);
        latw=latw(latw<=50&latw>=20);
        lonw=lonw(lonw<=140&lonw>=100);
        [LNw,LTw]=meshgrid(lonw,latw);
        ch2=m_windbarb(LNw,LTw,U_temp',V_temp',1,'units','m/s','linewi',1,'color','k');
        hold off;

        % m_line(world_mapx(:), world_mapy(:), 'color', [0 0 0], 'LineWidth', 0.5);

        m_grid('linestyle', 'none', 'tickdir', 'out', 'yaxislocation', 'left', 'fontsize', 10, 'color', 'k'); %????????????
        
        % ax=m_contfbar([0 1],-0.02,CS,ch,'endpiece','yes');
        % set(ax,'fontsize',12)
        % xlabel(ax,'Mean Daily Precipitation Rate/(kg/m^2/s)');
        % colormap(mycolor(22))
        legend([ch,ch1],'Temperature','Geopotential height','Location','northwest') %带方位
        title({['ERAi ',title_level{varLevel},'hPa ',title_time{varTime},'BJT']}); % cc=',num2str(corr))\
        [h1,h2,h3,h4]=legend_wb([0.85,.85,.12,.12],'4','k');
        % [h11,h21,h31]=m_legendvc( 'northeast',LNw,LTw,U_temp',V_temp',1,2);


        % outputPath=fullfile(nowpath,'figure');
        % f_tt = [title_level{varLevel},'hPa_','20190703_',title_time{varTime},'BJT'];
        % figurename=[outputPath,'/',f_tt,'.png'];
        % saveas(gcf,figurename)
 

        
        % ax=m_contfbar(0.05,[0 1],CS,ch, 'axfrac',.04,'endpiece','no','fontsize',12);
        % %([左 上下宽度],归一化高度（横色棒两者倒一下就变成竖坐标）CS,ch,'axfrac'左右宽度,'endpiece','no'去掉色块外小三角；'edgecolor','none/b'去掉框线，'levels':‘set’/'match'精准显示色阶还是模糊处理使数据和颜色对应更加规整

        % set(ax, 'xtick', a); %这个是颜色棒中具体的值
    end

end

% 20190703 14:00 and 20:00
% U_14=squeeze(U(:,:,:,10));
% V_14=squeeze(V(:,:,:,10));
% H_14=squeeze(H(:,:,:,10));
% T_14=squeeze(T(:,:,:,10));

% U_20=squeeze(U(:,:,:,11));
% V_20=squeeze(V(:,:,:,11));
% H_20=squeeze(H(:,:,:,11));
% T_20=squeeze(T(:,:,:,11));
