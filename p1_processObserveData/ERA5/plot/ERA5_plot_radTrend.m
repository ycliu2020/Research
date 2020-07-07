%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-07-06 15:05:35
% LastEditTime : 2020-07-07 11:20:56
% LastEditors  : LYC
% Description  :
% FilePath     : /code/p1_processObserveData/ERA5/plot/ERA5_plot_radTrend.m
%
%%---------------------------------------------------------
clc; clear;
nowpath = pwd;
% load mask map
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')
[mlabels, areaNum] = obsPlotParameters('toa', 'land', 'ERA5-radEffect-tsRad');
[readme, level, tLin, vars] = obsParameters();
% Latitude range
p_3 = 60;
lon1 = [2.5 357.5]; lat1 = [-p_3 + 1 p_3 - 1]; % world area
set(0, 'defaultfigurecolor', 'w')

%% different time series, 1mean 2000-03 to 2018-02(18*12). 2 mean 200207-201706(15*12)
for p_1 = 1:2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %read data
    varsPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{1}); %rawdata
    dvarsPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{2}); %anomaly
    dvarsTrendPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{3}); %anomaly_trend
    kernelCalPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{4}); % kernelCal
    radEfectPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{5}); %radEffect
    dradTrendPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{6}); %/data1/liuyincheng/cmip6-proces/aimp_2000-2014/MRI-ESM2-0/ensemble/radEffect_trend/
    outPutPath = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/01.observe/1.3/', tLin.time{p_1}, 'ERA5', 'radEffect');
    auto_mkdir(outPutPath)

    load([dvarsTrendPath, 'global_vars.mat'])%% 'lon_f', 'lat_f', 'lon_k', 'lat_k', 'plev_k', 'time'
    load([dradTrendPath, 'trend_dradEfect_toa_cld.mat'])% 10 vars:'trendyr_dRtoa_ta','trendyr_dRtoa_taOnly2', 'trendyr_dRtoa_tas2., 'trendyr_dRtoa_tsAtom', 'trendyr_dRtoa_mainEffect', 'trendyr_dRtoa_residual', 'trendyr_dRtoa_cloud', 'trendyr_dRtoa_q', 'trendyr_dRtoa_alb', 'trendyr_dRtoa_ts'
    load([dradTrendPath, 'trend_dradEfect_sfc_cld.mat'])% 10 vars:'trendyr_dRsfc_ta','trendyr_dRsfc_taOnly2', 'trendyr_dRsfc_tas2., 'trendyr_dRsfc_tsAtom', 'trendyr_dRsfc_mainEffect', 'trendyr_dRsfc_residual', 'trendyr_dRsfc_cloud', 'trendyr_dRsfc_q', 'trendyr_dRsfc_alb', 'trendyr_dRsfc_ts'
    load([dvarsTrendPath, 'trend_dnetTOA.mat'])% trendyr_dnetTOA
    load([dvarsTrendPath, 'trend_drhs.mat'])% trendyr_drhs
    load([dvarsTrendPath, 'trend_dhFlux.mat'])% trendyr_dhFlux
    load([dvarsTrendPath, 'trend_drhsPlus.mat'])% trendyr_dhFlux
    load([dvarsTrendPath, 'trend_dts.mat'])% trendyr_dts
    nlonf = length(lon_f);
    nlatf = length(lat_f);
    % use one var to plot
    trendyr = zeros(nlonf, nlatf, 9);

    for jj = 1:9
        eval(['trendyr(:,:,jj)=', mlabels.vars{jj}, ';'])
    end

    trendyr = trendyr * 365 * 10;
    trendyr(:, :, 1) = trendyr(:, :, 1);
    % mask and cal the cc
    [trendyr, yr_cc, yr_pp] = maskArea(trendyr, lat_f, p_3, -p_3, areaNum);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot
    f_row = 3; f_col = 3; % 设置画图的行列
    set(0, 'DefaultFigureVisible', 'on')

    ss = get(0, 'ScreenSize');
    h = figure('Position', [ss(4)/2-100 ss(3) / 35 ss(3)/5*f_col (ss(4)-80)/5*f_row]);
    % h = figure('Position', [ss(4)/2-100 ss(3) / 35 ss(3)/2+100 ss(4) * 4/5]);
    % set(h,'visible','off');
    % clf reset;
    set(h, 'Color', [1 1 1]);
    f_matrix = reshape(1:f_row * f_col, [f_col, f_row])';
    % figure
    for varNum = 1:f_row * f_col
        trendz = squeeze(trendyr(:, :, varNum));
        [plotRow, plotCol] = find(f_matrix == varNum);
        subplot_yc(3, 3, plotRow, plotCol);
        hold on
        m_proj('Mercator', 'lon', lon1, 'lat', lat1); %Mercator,Equidistant cylindrical,lambert,Miller Cylindrical
        m_pcolor(lon_f, lat_f, trendz');
        % [x, y]=find(trendz(:,:)<=0);
        % m_plot(lon(x),lat(y),'g.','markersize',5,'color','k');
        % [Xlon,Ylat] = meshgrid(lon,lat);
        % m_contourf(Xlon,Ylat,trendz',30,'LineStyle','none')
        colormap(mycolor(18)); %mycolor(100)is soden color????????colormap(flipud(mycolor(13)));%colormap(jet(4))
        col_SeriesNum = 10;
        [colorbar_Series] = findSuit_colorInt(trendz, col_SeriesNum);
        % caxis([mmin(varNum) mmax(varNum)]);
        caxis([min(colorbar_Series) max(colorbar_Series)]);
        hold on
        m_line(world_mapx(:), world_mapy(:), 'color', [0 0 0], 'LineWidth', 0.5);

        if plotCol == 1 && plotRow == f_row
            m_grid('linestyle', 'none', 'tickdir', 'out', 'yaxislocation', 'left', 'fontsize', 8, 'color', 'k');
        elseif plotCol == 1 && plotRow ~= f_row
            m_grid('linestyle', 'none', 'tickdir', 'out', 'xticklabels', [], 'yaxislocation', 'left', 'fontsize', 8, 'color', 'k');
        elseif plotCol ~= 1 && plotRow == f_row
            m_grid('linestyle', 'none', 'tickdir', 'out', 'yticklabels', [], 'fontsize', 8, 'color', 'k');
        elseif plotCol ~= 1 && plotRow ~= f_row
            m_grid('linestyle', 'none', 'tickdir', 'out', 'xticklabels', [], 'yticklabels', [], 'fontsize', 8, 'color', 'k');
        end

        title({[mlabels.component{varNum}, mlabels.unite{varNum}]; ['cc = ', num2str(yr_cc{varNum})]}, 'Interpreter', 'latex', 'fontsize', 10); % cc=',num2str(corr))
        % c=colorbar;
        % % c.Limits=[mmin(varNum) mmax(varNum)];%
        % c.Box='off';
        hold on
        c = colorbar;
        % c.TickLength = 0.0245;
        c.Ticks = colorbar_Series(2:end - 1);
        c.Limits = [min(colorbar_Series) max(colorbar_Series)];
    end

    tt = {['Level:', mlabels.level, ', Era: ', tLin.time{p_1}], ['Data:', mlabels.dataName{1}, ', Trend(year mean)']};
    sgtt = sgtitle(tt, 'Fontsize', 14, 'Interpreter', 'none');
    figureName = [mlabels.dataName{1}, '_', tLin.time{p_1}, '_radEffect', '_', mlabels.area];
    saveFileName = [outPutPath, '/', figureName, '.png'];
    saveas(gcf, saveFileName)
    % save_png(saveFileName)%high resolution
end
