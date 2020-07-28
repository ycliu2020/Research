%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% mothly data
% plot the Ts trend and Rheating, cloud, Ta, q, albedo trend
% raw data: ERAi
% add 200207_201706 cloud Effect trend by MODIS data using kernel method
% time series: 200207_201706
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

clc; clear; tic;
inputpath_Rheating = '/data1/liuyincheng/Observe-process/200207-201706/ERAi/anomaly_trend/trend_dtsdradERAi.mat'; %cons_rhs, cons_ts, lat, lon, p_rhs, p_ts, time, trendm_drhs, trendm_dts, trendyr_drhs, trendyr_dts
inputpath_component = '/data1/liuyincheng/Observe-process/200207-201706/ERAi/radEffect_trend/trend_dradEfect_all_ERAi.mat'; %latPlot, lonPlot, time, trendyr_dRsfc_alb, trendyr_dRsfc_cloud, trendyr_dRsfc_q, trendyr_dRsfc_realRad, trendyr_dRsfc_ta, trendyr_dRsfc_total, trendyr_dRsfc_total_c, trendyr_dRsfc_ts, trendyr_dRtoa_alb, trendyr_dRtoa_cloud, trendyr_dRtoa_q, trendyr_dRtoa_realRad, trendyr_dRtoa_ta, trendyr_dRtoa_total, trendyr_dRtoa_total_c, trendyr_dRtoa_ts
inputpath_cloudMOD = '/data1/liuyincheng/Observe-process/200207-201706/MODIS/radEffect_trend/trend_dcloud_MOD.mat'; %latPlot, lonPlot, readme, time, trendm_dRsfc_cloudMOD, trends_dRsfc_cloudMOD, trendyr_dRsfc_cloudMOD,trendyr_dRsfc_residualMOD',trendyr_dRsfcHeat_varsCal_Resi',trendyr_dRsfcHeat_varsCal_noResi',trendyr_testVar_equal_Resi','trendyr_dRsfcHeat_varsCal_ResiERAi'
load(inputpath_Rheating)
load(inputpath_component)
load(inputpath_cloudMOD)

outputPath = '/home/lyc/research/P02.Ts_change_research/figure/Observe/';
auto_mkdir(outputPath)

load('/home/lyc/lib/tools/matlab/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/lyc/lib/tools/matlab/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/lyc/lib/tools/matlab/myMap/02.world_map/mat_file/correct_worldmap.mat')%
load('/home/lyc/lib/tools/matlab/myMap/01.china_map/mat_file/mask14472.mat')

month = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
season = {'MAM', 'JJA', 'SON', 'DJF'};
unite = {' trend(W\cdotm-2/10a)',' trend(W\cdotm-2/10a)', ...
    ' trend(W\cdotm-2/10a)',' trend(W\cdotm-2/10a)', ...
    ' trend(W\cdotm-2/10a)',' trend(W\cdotm-2/10a)'};

component = {'Ts', 'RHeating radiation', 'Cloud radeffect', 'Ta radeffect', 'q radeffect', 'albedo radeffect'};
component_label = {'dRsfc_ts', 'drhs','dRsfc_cloudMOD', 'dRsfc_ta', 'dRsfc_q',  'dRsfc_alb'};
vars_label = strcat('trendyr_', component_label);
dataname = {'ERAi', 'CERES', 'MODIS'};
level_label = {'SFC', 'TOA', 'ATM(TOA-SFC)'};

% figure set
set(0, 'defaultfigurecolor', 'w')
areaNum = 0; % world land
latRange = 60; % Latitude range
lon1 = [2.5 357.5]; lat1 = [-latRange + 1 latRange - 1]; % world area
% lon1=[70,140];lat1=[0,60]; % China area
% colorRange={[-5 -5 -4 -4 -2 -2];[-5 -5 -4 -4 -1 -1];[-3 -4 -1 -2 -.5 -1];[-5 -7 -1.5 -3 -1 -1]};
mmin = -4; %colorRange{p_1};
mmax = -mmin;

% read plot data(use one var to plot)
nlon=length(lonPlot);nlat=length(latPlot);
trendyr = zeros(nlon, nlat, 6);

for jj = 1:6
    eval(['trendyr(:,:,jj)=', vars_label{jj}, ';'])
end

trendyr = trendyr * 365 * 10;
% trendyr(:, :, 1) = trendyr(:, :, 1) .* 5;
% mask and cal the cc
[trendyr, yr_cc, yr_pp] = maskArea(trendyr, lat, latRange, -latRange, areaNum);

set(0, 'DefaultFigureVisible', 'on')
ss = get(0, 'ScreenSize'); % ???????????
h = figure('Position', [ss(4) / 2 ss(3) / 35 ss(3) * 3/9.5 ss(4) * 4/5]); % ???????????????

set(h, 'Color', [1 1 1]);
f_matrix = reshape(1:6, [2, 3])';
p_4=[1 1 3 1 1 1];
% figure
for ii = 1:6
    trendz = squeeze(trendyr(:, :, ii));
    [i, j] = find(f_matrix == ii);

    subplot_yc(3, 2, i, j); % ?????
    hold on
    m_proj('Mercator', 'lon', lon1, 'lat', lat1); %?????????????Mercator,Equidistant cylindrical,lambert,Miller Cylindrical
    m_pcolor(lonPlot, latPlot, trendz');
    %contourf
    % [Xlon,Ylat] = meshgrid(lon,lat);
    % m_contourf(Xlon,Ylat,trendz',30,'LineStyle','none')
    colormap(mycolor(18)); %????????colormap(flipud(mycolor(13)));%colormap(jet(4))%????????????????????
    % caxis([mmin(ii) mmax(ii)]);
    caxis([mmin mmax]);
    hold on
    m_line(world_mapx(:), world_mapy(:), 'color', [0 0 0], 'LineWidth', 0.5);
    m_grid('linestyle', 'none', 'tickdir', 'out', 'yaxislocation', 'left', 'fontsize', 10, 'color', 'k'); %????????????
    title({[dataname{p_4(ii)},' ',component{ii}, unite{ii}]; ['year mean (cc = ', num2str(yr_cc{ii}), ')']}); % cc=',num2str(corr))
    % c=colorbar;
    % % c.Limits=[mmin(ii) mmax(ii)];
    % c.Box='off';
    hold on
    pos = get(gca, 'Position');

end

c = colorbar('southoutside', 'Position', [pos(1) - 0.35 pos(2) - 0.06 0.6 pos(4) / 9]);
c.TickLength = 0.026;
% c.LineWidth = 2;
c.Limits = [mmin mmax];
if mmin == -5
    c.Ticks = [-4, -3, -2, -1, 0, 1, 2, 3, 4];
elseif mmin == -4
    c.Ticks = [-3.2, -2.4, -1.6, -.8, 0, .8, 1.6, 2.4, 3.2];
elseif mmin == -3
    c.Ticks = [-2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5];
elseif mmin == -2.5
    c.Ticks = [-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2];
end
% c.LineWidth = 'white';
% c.Box='off';
tt = ['Level:', level_label{1}, ', Era: 200207-201706', ];
sgtt = sgtitle(headLineTxt, 'Fontsize', 14, 'Interpreter', 'none');
% save file 
% figTitle = [level.time1{p_1}(1:end - 1), '_', level.model2{level1}];
figurename = '/home/lyc/research/P02.Ts_change_research/figure/Observe/1.2/sfc_cloudMODIS_landsea.png';
saveas(gcf, figurename)
% save_png(figurename)%high resolution
close gcf
