%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-07-06 15:05:35
% LastEditTime : 2020-09-15 21:13:56
% LastEditors  : LYC
% Description  :
% FilePath     : /code/p1_processObserveData/ERA5/plot/ERA5_plot_radTrend_fig2.m
%
%%---------------------------------------------------------
clc; clear;
nowpath = pwd;
% load mask map
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')
[mlabels, areaNum] = obsPlotParameters('sfc', 'land', 'ERA5-radEffect-tsRad');
[readme, level, tLin, vars] = obsParameters('ERA5');
% Latitude range
latRange = 60;
lon1 = [2.5 357.5]; lat1 = [-latRange + 1 latRange - 1]; % world area
set(0, 'defaultfigurecolor', 'w')

%% different time series, 1mean 2000-03 to 2018-02(18*12). 2 mean 200207-201706(15*12)
for exmNum = 1:2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %read data
    varsPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{1}); %rawdata
    dvarsPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{2}); %anomaly
    dvarsTrendPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{3}); %anomaly_trend
    kernelCalPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{4}); % kernelCal
    radEfectPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{5}); %radEffect
    dradTrendPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{6}); %/data1/liuyincheng/cmip6-proces/aimp_2000-2014/MRI-ESM2-0/ensemble/radEffect_trend/
    outPutPath = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/proj1_Observe/Radiative_effect/v0.1/',  ['ERA5_',tLin.time{exmNum}]);
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

    for varNum = 1:9
        eval(['trendyr(:,:,varNum)=', mlabels.vars{varNum}, ';'])
    end

    trendyr = trendyr * 365 * 10;
    trendyr(:, :, 1) = trendyr(:, :, 1);
    % mask and cal the cc
    [trendyr, yr_cc, yr_pp] = maskArea(trendyr, lat_f, latRange, -latRange, areaNum);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot
    f_row = 1; f_col = 2; % 设置画图的行列
    set(0, 'DefaultFigureVisible', 'on')

    ss = get(0, 'ScreenSize');
    h = figure('Position', [96 433 1361 455]);
    % h = figure('Position', [ss(4)/2-100 ss(3) / 35 ss(3)/2+100 ss(4) * 4/5]);
    % set(h,'visible','off');
    % clf reset;
    set(h, 'Color', [1 1 1]);
    f_matrix = reshape(1:f_row * f_col, [f_col, f_row])';
    % figure
    for varNum = 1:f_row * f_col

        [plotRow, plotCol] = find(f_matrix == varNum);
        subplot_yc(f_row, f_col, plotRow, plotCol);
        hold on
        if varNum==2
            varNum=3;
        end
        trendz = squeeze(trendyr(:, :, varNum));

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
        cax_limt=5;
        caxis([-cax_limt cax_limt]);
        hold on
        m_line(world_mapx(:), world_mapy(:), 'color', [0 0 0], 'LineWidth', 0.5);
        % plot ticklabels only on the margin
        % if plotCol == 1 && plotRow == f_row
        %     m_grid('linestyle', 'none', 'tickdir', 'out', 'yaxislocation', 'left', 'fontsize', 8, 'color', 'k');
        % elseif plotCol == 1 && plotRow ~= f_row
        %     m_grid('linestyle', 'none', 'tickdir', 'out', 'xticklabels', [], 'yaxislocation', 'left', 'fontsize', 8, 'color', 'k');
        % elseif plotCol ~= 1 && plotRow == f_row
        %     m_grid('linestyle', 'none', 'tickdir', 'out', 'yticklabels', [], 'fontsize', 8, 'color', 'k');
        % elseif plotCol ~= 1 && plotRow ~= f_row
        %     m_grid('linestyle', 'none', 'tickdir', 'out', 'xticklabels', [], 'yticklabels', [], 'fontsize', 8, 'color', 'k');
        % end
        m_grid('linestyle', 'none', 'tickdir', 'out', 'yaxislocation', 'left', 'fontsize', 8, 'color', 'k');

        title({[mlabels.component{varNum}, mlabels.unite{varNum}]; ['spatial cc = ', num2str(yr_cc{varNum})]}, 'Interpreter', 'latex', 'fontsize', 10);
        % sequence
        figSeq = char(97 + varNum - 1);
        yLoc = -0.3;
        if plotRow == f_row
            yLoc = -0.3;
        end
        txt = text(0.5, yLoc, ['(', figSeq, ')'], 'Units', 'normalized', 'FontSize', 10);
        % c=colorbar;
        % % c.Limits=[mmin(varNum) mmax(varNum)];%
        % c.Box='off';
        hold on
        c = colorbar;
        % c.TickLength = 0.0245;
        c.Ticks = [-5 -2.5 0 2.5 5];
        c.Limits = [-cax_limt cax_limt];
    end

    headLineTxt = {['Level:', mlabels.level, ', Era: ', tLin.time{exmNum}], ['Data:', mlabels.dataName{1}, ', Trend(year mean)']};
    sgtt = sgtitle(headLineTxt, 'Fontsize', 14, 'Interpreter', 'none');
    figureName = [mlabels.dataName{1}, '_', tLin.time{exmNum}, '_', mlabels.figType{1}, '_', mlabels.area, '_', mlabels.level,'_fig2'];
    saveFileName = [outPutPath, '/', figureName, '.png'];
    saveas(gcf, saveFileName)
    % save_png(saveFileName)%high resolution
    close gcf

end
