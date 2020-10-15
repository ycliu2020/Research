%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-09-13 19:06:17
% LastEditTime : 2020-09-23 09:53:59
% LastEditors  : LYC
% Description  : paper用图: ERA5 dataset timeSeries analyse(dR_ts and rhs)
% FilePath     : /code/p1_processObserveData/ERA5/timeSeriAnalysis/ERA5_timeSeriANS_tsVsVarsRad.m
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

latRange = 60; % Latitude range
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
    load([radEfectPath, 'dR_residual_cld_sfc'])% load dR_residual_cld_sfc
    load([radEfectPath, 'dR_cloud_sfc'])% load dR_cloud_sfc
    load([dvarsPath, 'drhs.mat'])% load drhs
    nlonf = length(lon_f);
    nlatf = length(lat_f);
    ntime = length(time);
    % regrid
    drhs = autoRegrid3(lon_k, lat_k, time, drhs, lon_f, lat_f, time);

    varUsed = zeros(nlonf, nlatf, ntime, 2);
    % varUsed(:, :, :, 1) = tsEffect;
    % varUsed(:, :, :, 2) = drhs;

    % fig1
    dRheating = dR_cloud_sfc + husEffect + taEffect + albEffect + dR_residual_cld_sfc;

    varUsed(:, :, :, 1) = tsEffect;
    varUsed(:, :, :, 2) = dR_cloud_sfc;
    varUsed(:, :, :, 3) = husEffect;
    varUsed(:, :, :, 4) = taEffect;
    varUsed(:, :, :, 5) = albEffect;
    varUsed(:, :, :, 6) = dR_residual_cld_sfc;
    % varUsed(:, :, :, 7) = dR_netSfc; %sum(varUsed(:,:,:,1:6),4);
    % varUsed(:, :, :, 8) = tsEffect+dRheating;
    % varUsed(:, :, :, 9) = dRheating;

    varNames = {'dR_{Ts}', 'dR_{cloud}', 'dR_{q}', 'dR_{ta}', 'dR_{albedo}', 'dR_{residual}'};
    yLabel = {'Wm-2', 'Wm-2', 'Wm-2', 'Wm-2', 'Wm-2', 'Wm-2', 'Wm-2', 'Wm-2', 'Wm-2'};
    varColor = {'#000000', '#4450A1', '#F08212', '#C2162E', '#90C64D', '#67C1EE'}; %黑色 深蓝 橘色 红色 浅绿 浅蓝

    % fig2
    % dRheating=dR_cloud_sfc+husEffect+taEffect+albEffect+dR_residual_cld_sfc;
    % varUsed(:, :, :, 1) = tsEffect;
    % varUsed(:, :, :, 2) = dRheating;
    % varUsed(:, :, :, 3) = -dR_netSfc;
    % varUsed(:, :, :, 4) = dhFlux;
    % varUsed(:, :, :, 5) = tsEffect+dRheating-dR_netSfc;
    % varNames = {'dR_{Ts}', 'dR_{heating}', 'dR_{netSFC}', 'dR_{LH+SH}','dR_{ts+Rheating+dRnetSFC}'};
    % yLabel = {'Wm-2', 'Wm-2', 'Wm-2', 'Wm-2', 'Wm-2', 'Wm-2', 'Wm-2', 'Wm-2', 'Wm-2', 'Wm-2', 'Wm-2'};
    % varColor = {'#000000', '#4450A1', '#F08212', '#C2162E', '#90C64D'}; %黑色 深蓝 橘色 红色 浅绿 浅蓝

    % transfor to yearly Data
    sizeVarUsed = size(varUsed);
    sizeLon = sizeVarUsed(1); sizeLat = sizeVarUsed(2); sizeTime = sizeVarUsed(3); sizeVar = sizeVarUsed(4);
    varUsedYearly = zeros(sizeLon, sizeLat, sizeTime / 12, sizeVar);

    for varNum = 1:sizeVar
        varUsedYearly(:, :, :, varNum) = monthlyToYearly(squeeze(varUsed(:, :, :, varNum)));
    end

    % mask only -60-60 land.
    for varNum = 1:sizeVar
        [varUsedYearly(:, :, :, varNum), ~, ~] = maskArea(squeeze(varUsedYearly(:, :, :, varNum)), lat_f, latRange, -latRange, 'world');
    end

    % % method 1 cal areaMeanLatWeight
    sizeVarUsedYearly = size(varUsedYearly);
    varUsedYearly_weightMean = zeros(sizeVarUsedYearly(3), sizeVar);

    for varNum = 1:sizeVar

        for timeNum = 1:sizeVarUsedYearly(3)
            varUsedYearly_weightMean(timeNum, varNum) = areaMeanLatWeight(squeeze(squeeze(varUsedYearly(:, :, timeNum, varNum))), lat_f);
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot time series and CC
    % time series
    timeSer = [1985 1990 1995 2000 2005 2010 2015 2020 2025 2035 2045 2055 2065 2075 2085 2095 2105];
    char_timeSer = cellstr(string(timeSer));
    % fig set
    set(0, 'DefaultFigureVisible', 'on')
    ss = get(0, 'ScreenSize');
    coef_amplify = 1.5;
    % h = figure('Position', [ss(4) / 2 - 100 ss(3) / 35 ss(3) / 5 * f_col * coef_amplify (ss(4) - 80) / 5 * f_row * coef_amplify]);
    h = figure('Position', [84 246 889 514]);
    % clf reset;
    set(h, 'Color', [1 1 1]);

    for varNum = 1:sizeWeightMean(2)
        varUsedTemp = squeeze(varUsedYearly_weightMean(:, varNum));
        plot(timeYearly, varUsedTemp', 'color', varColor{varNum}, 'LineWidth', 2)
        hold on

    end

    % zero line
    zeroLine = zeros(length(timeYearly), 1);
    plot(timeYearly, zeroLine', 'k--')
    hold on;

    % set x axes
    xticks(timeSer)
    datetick('x', 'yyyy', 'keepticks'); % 用日期的格式显示横坐标
    xlim([timeYearly(1) timeYearly(end)])
    xlabel('time line')
    xticklabels(char_timeSer)

    ylabel(yLabel{2})
    % set y axes
    ymax = 4;

    if exmNum == 1
        ymax = 0.8;
    elseif exmNum == 2
        ymax = 0.8;
    end

    ylim([-ymax ymax])

    ax = gca;
    ax.XMinorTick = 'on'; ax.YMinorTick = 'on'; % 开启次刻度线
    ax.TickLength = [0.02 0.01]; %刻度线长度      set(gca,'ticklength', [0.02 0.01]);
    ax.XColor = 'k'; ax.YColor = 'k'; % 设置刻度线颜色
    title_txt = {['Data: ', mlabels.dataName{1}, ', Level: ', mlabels.level, ', Era: ', tLin.time{exmNum}], '60N-60S land, global mean '};
    title(title_txt, 'FontWeight', 'normal')

    lgdTxtTmp(1:6)={',cc='};
    cc_weightMeanTxt=cellfun(@num2str,cc_weightMean,'UniformOutput',false);
    lgdTxt=cellfun(@strcat,varNames,lgdTxtTmp,cc_weightMeanTxt,'UniformOutput',false);
    lgd = legend(lgdTxt, 'Fontsize', 8, 'Location', 'northwest', 'NumColumns', 3);
    % lgd = legend(varNames, 'Fontsize', 8, 'Location', 'northwest', 'NumColumns', 3);
    legend('boxoff')%删除图例背景和轮廓
    lgd_inf = get(lgd);
    % text(lgd_inf.Position(1) - 0.12, lgd_inf.Position(2) + 0.085, ['cc= ', num2str(cc_weightMean{2})], 'FontWeight', 'normal', 'Interpreter', 'none', 'Units', 'normalized')
    % hold on

    % save figures
    figName = [mlabels.dataName{1}, '_', tLin.time{exmNum}, '_', mlabels.area, '_', mlabels.level, '_globalLandMean_ts&VarsRad'];
    figPath = [outPutPath, '/', figName, '.png'];
    saveas(gcf, figPath)
    % save_png(figPath)%high resolution
    % close gcf

end
