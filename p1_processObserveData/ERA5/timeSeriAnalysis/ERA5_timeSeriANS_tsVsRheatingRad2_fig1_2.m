%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-09-13 19:06:17
% LastEditTime : 2020-11-16 22:16:50
% LastEditors  : LYC
% Description  : paper用图: ERA5 dataset timeSeries analyse(dR_ts and rhs)
% FilePath     : /code/p1_processObserveData/ERA5/timeSeriAnalysis/ERA5_timeSeriANS_tsVsRheatingRad2_fig1_2.m
%
%%---------------------------------------------------------
clc; clear;
nowpath = pwd;
% load mask map
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')
[mlabels, areaNum] = obsPlotParameters('sfc', 'land', 'ERA5-radEffect-ts');
[readme, level, tLin, vars] = obsParameters('ERA5');

latRange = 90; % Latitude range
lon1 = [2.5 357.5]; lat1 = [-latRange + 1 latRange - 1]; % world area
toaSfc = {'toa', 'sfc'};
lon_k = 0:2.5:357.5; nlonk = length(lon_k); % kernel lat lon
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);

for exmNum = 1:2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %read data
    varsPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{1}); %rawdata
    dvarsPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{2}); %anomaly
    dvarsTrendPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{3}); %anomaly_trend
    kernelCalPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{4}); % kernelCal
    radEfectPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{5}); %radEffect
    dradTrendPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{6}); %/data1/liuyincheng/cmip6-proces/aimp_2000-2014/MRI-ESM2-0/ensemble/radEffect_trend/
    outPutPath = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/proj1_Observe/TimeSeries_analysis/', ['ERA5_', tLin.time{exmNum}]);
    auto_mkdir(outPutPath)
    load([dvarsTrendPath, 'global_vars.mat'])%% 'lon_f', 'lat_f', 'lon_k', 'lat_k', 'plev_k', 'time'
    load([radEfectPath, 'dradEfect_sfc_cld.mat'])% load albEffect, husEffect, mainEffect, taEffect, taOnlyEffect, taOnlyEffect2, tasEffect, tasEffect2, totalEffect, tsEffect, wvlwEffect, wvswEffect
    load([dvarsPath, 'drhs.mat'])% load drhs
    nlonf = length(lon_f);
    nlatf = length(lat_f);
    ntime = length(time);
    % regrid
    drhs = autoRegrid3(lon_k, lat_k, time, drhs, lon_f, lat_f, time);

    varUsed = zeros(nlonf, nlatf, ntime, 2);
    varUsed(:, :, :, 1) = tsEffect;
    varUsed(:, :, :, 2) = drhs;
    varNames = {'dR_{Ts}', 'dRHeating', 'ta RadEfect', 'ts RadEfect', 'q RadEfect', 'alb RadEfect'};
    yLabel = {'W/m2', 'W/m2', 'W/m2', 'W/m2', 'W/m2', 'W/m2'};

    sizeVarUsed = size(varUsed);
    sizeLon = sizeVarUsed(1); sizeLat = sizeVarUsed(2); sizeTime = sizeVarUsed(3); sizeVar = sizeVarUsed(4);
    varUsedYearly = zeros(sizeLon, sizeLat, sizeTime / 12, sizeVar);

    for varNum = 1:sizeVar
        varUsedYearly(:, :, :, varNum) = monthlyToYearly(varUsed(:, :, :, varNum));
    end

    % mask only -60-60 land.
    for varNum = 1:sizeVar
        [varUsedYearly(:, :, :, varNum), ~, ~] = maskArea(squeeze(varUsedYearly(:, :, :, varNum)), lat_f, latRange, -latRange, 'world');
    end

    % cal areaMeanLatWeight
    sizeVarUsedYearly = size(varUsedYearly);
    varUsedYearly_weightMean = zeros(sizeVarUsedYearly(3), sizeVar);

    for varNum = 1:sizeVar

        for timeNum = 1:sizeVarUsedYearly(3)
            varUsedYearly_weightMean(timeNum, varNum) = areaMeanLatWeight(varUsedYearly(:, :, timeNum, varNum), lat_f);
        end

    end

    ear5_time = time;
    clear time
    time.vec = datevec(ear5_time);
    timeYearly = unique(time.vec(:, 1));
    timeYearly = timeYearly(1:end - 1);
    % cal cc of areaMeanLatWeight
    sizeWeightMean = size(varUsedYearly_weightMean);
    cc_weightMean = cell(1, sizeWeightMean(2));
    pp_weightMean = cell(1, sizeWeightMean(2));

    for varNum = 1:sizeWeightMean(2)
        [cc0, pp0] = corrcoef(varUsedYearly_weightMean(:, 1), varUsedYearly_weightMean(:, varNum), 'Rows', 'complete');
        cc_weightMean{varNum} = roundn(cc0(1, 2), -2); %保留俩位小数
        pp_weightMean{varNum} = pp0(1, 2); % confidence interval
    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % plot time series and CC
    % % time series
    % timeSer = [1985 1990 1995 2000 2005 2010 2015 2020 2025 2035 2045 2055 2065 2075 2085 2095 2105];
    % char_timeSer = cellstr(string(timeSer));
    % % fig set
    % f_row = 1; f_col = 1; % 设置画图的行列
    % set(0, 'DefaultFigureVisible', 'on')
    % ss = get(0, 'ScreenSize');
    % coef_amplify = 1.5;
    % % h = figure('Position', [ss(4) / 2 - 100 ss(3) / 35 ss(3) / 5 * f_col * coef_amplify (ss(4) - 80) / 5 * f_row * coef_amplify]);
    % h = figure('Position', [317 203 672 384]);
    % % clf reset;
    % set(h, 'Color', [1 1 1]);
    % % zero line
    % zeroLine = zeros(length(timeYearly), 1);

    % varUsedTemp1 = squeeze(varUsedYearly_weightMean(:, 1));
    % varUsedTemp2 = squeeze(varUsedYearly_weightMean(:, 2));

    % plot(timeYearly, varUsedTemp1', 'b', 'LineWidth', 1.5)
    % hold on
    % plot(timeYearly, varUsedTemp2', 'r', 'LineWidth', 1.5)
    % hold on
    % plot(timeYearly, zeroLine', 'k--')
    % hold on;
    % % set x axes
    % xticks(timeSer)
    % datetick('x', 'yyyy', 'keepticks'); % 用日期的格式显示横坐标
    % xlim([timeYearly(1) timeYearly(end)])
    % xlabel('time line')
    % xticklabels(char_timeSer)

    % ylabel(yLabel{2})
    % % set y axes
    % ymax = 5;

    % if strcmp(varNames{2}, 'dTs') == 1
    %     ymax = 0.5;
    % end

    % ylim([-ymax ymax])

    % ax = gca;
    % ax.XMinorTick = 'on'; ax.YMinorTick = 'on'; % 开启次刻度线
    % ax.TickLength = [0.02 0.01]; %刻度线长度      set(gca,'ticklength', [0.02 0.01]);
    % ax.XColor = 'k'; ax.YColor = 'k'; % 设置刻度线颜色
    % title({['Data: ', mlabels.dataName{1}, ', Level: ', mlabels.level, ', Era: ', tLin.time{exmNum}], '60N-60S land, global mean', [varNames{2}, ' & ', varNames{1}, ' cc=', num2str(cc_weightMean{2})]}, 'FontWeight', 'normal')
    % % title(['60N-60S land, global mean: ', varNames{2}, '&', varNames{1}, ' cc=', num2str(cc_weightMean{2})], 'FontWeight', 'normal')

    % % ylabel({Level{i}, ' Rad Anomaly (Wm^{-2})'}, 'Fontsize', 10)
    % lgd = legend('dR_{Ts}', 'dRHeating', 'Fontsize', 8, 'Location', 'northwest');
    % legend('boxoff')%删除图例背景和轮廓
    % lgd_inf = get(lgd);
    % text(lgd_inf.Position(1) - 0.12, lgd_inf.Position(2) + 0.095, ['cc= ', num2str(cc_weightMean{2})], 'FontWeight', 'normal', 'Interpreter', 'none', 'Units', 'normalized')
    % hold on

    % % save figures
    % figName = [mlabels.dataName{1}, '_', tLin.time{exmNum}, '_', mlabels.area, '_', mlabels.level, '_globalLandMeanTimeSeries'];
    % figPath = [outPutPath, '/', figName, '.png'];
    % saveas(gcf, figPath)
    % save_png(figPath)%high resolution
    % close gcf

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % cal cc of every point of map
    cc_global = zeros(nlonf, nlatf, sizeVar);
    pp_global = zeros(nlonf, nlatf, sizeVar);

    for varNum = 1:sizeVar

        for latNum = 1:nlatf

            for lonNum = 1:nlonf
                varUsedTemp = squeeze(squeeze(varUsedYearly(lonNum, latNum, :, :)));
                [cc0, pp0] = corrcoef(varUsedYearly(lonNum, latNum, :, 1), varUsedYearly(lonNum, latNum, :, varNum), 'Rows', 'complete');
                cc_global(lonNum, latNum, varNum) = roundn(cc0(1, 2), -2); %保留俩位小数
                pp_global(lonNum, latNum, varNum) = pp0(1, 2); % confidence interval
            end

        end

    end

    % plot cc global dustribution
    f_row = 1; f_col = 1; % 设置画图的行列
    set(0, 'DefaultFigureVisible', 'on')
    ss = get(0, 'ScreenSize');
    coef_amplify = 3;
    h = figure('Position', [ss(4) / 2 - 100 ss(3) / 35 ss(3) / 5 * f_col * coef_amplify (ss(4) - 80) / 5 * f_row * coef_amplify]);
    % clf reset;
    set(h, 'Color', [1 1 1]);
    f_matrix = reshape(1:f_row * f_col, [f_col, f_row])';

    % figure
    for varNum = 1:f_row * f_col
        trendz = squeeze(cc_global(:, :, 2));
        [plotRow, plotCol] = find(f_matrix == varNum);
        subplot_yc(f_row, f_col, plotRow, plotCol);
        hold on
        m_proj('Mercator', 'lon_f', lon1, 'lat_f', lat1); %Mercator,Equidistant cylindrical,lambert,Miller Cylindrical
        m_pcolor(lon_f, lat_f, trendz');
        colormap(mycolor(18)); %mycolor(100)is soden color????????colormap(flipud(mycolor(13)));%colormap(jet(4))
        col_SeriesNum = 10;
        % [colorbar_Series] = findSuit_colorInt(trendz, col_SeriesNum);
        max_color = 1;
        min_color = -max_color;
        caxis([min_color max_color]);
        hold on
        m_line(world_mapx(:), world_mapy(:), 'color', [0 0 0], 'LineWidth', 0.5);
        m_grid('linestyle', 'none', 'tickdir', 'out', 'yaxislocation', 'left', 'fontsize', 8, 'color', 'k');

        headLineTxt = {'The correlation coefficient of dRTs and dRHeating', ['Data: ', mlabels.dataName{1}, ', Level: ', mlabels.level, ', Era: ', tLin.time{exmNum}], ['global mean CC= ', num2str(cc_weightMean{2})]};
        title(headLineTxt, 'Interpreter', 'none', 'fontsize', 14); % cc=',num2str(corr))
        hold on
        c = colorbar;
        c.Location='southoutside';

        % c.TickLength = 0.0245;
        % c.Limits = [min_color min_color];
    end

    % save figures
    figName = [mlabels.dataName{1}, '_', tLin.time{exmNum}, '_', mlabels.area, '_', mlabels.level, '_globalDistribution'];
    figPath = [outPutPath, '/', figName, '.png'];
    saveas(gcf, figPath)
    % save_png(figPath)%high resolution
    % close gcf

end
