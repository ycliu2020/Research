%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-08-31 17:00:15
% LastEditTime : 2021-02-09 00:10:20
% LastEditors  : LYC
% Description  : 同时画时间序列和相关性分布图
% FilePath     : /P02.Ts_change_research/code/p3_paperFigIntegrate/Fig3_timeSeries_tsVsRheatingRad/CMIP6_plot_fig3.m
%
%%---------------------------------------------------------
clc; clear; tic;
nowpath = pwd;
% load mask map
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')

toaSfc = {'toa', 'sfc'};
lon_k = 0:2.5:357.5; nlonk = length(lon_k); % kernel lat lon
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment
for exmNum = 1:2%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    %CAMS-CSM1-0 didn't have sfc clear sky radiation, delete it
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % exmPath
    exmPath = ['/data1/liuyincheng/cmip6-process/', level.time1{exmNum}]; %/data1/liuyincheng/cmip6-process/2000-2014/

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model
    for mdlNum = 3:3%1:length(level.model2)
        % model path
        mdlName = level.model2{mdlNum};
        mdlPath = fullfile(exmPath, level.model2{mdlNum});
        eval(['cd ', mdlPath]);
        disp(' ')
        disp([level.model2{mdlNum}, ' model start!'])

        % ensemble member path
        esmName = getPath_fileName(mdlPath, '.');
        eval(['cd ', nowpath]);

        %% 暂时只看esm实验
        esm = 'r1i1p1f1';

        if sum(strcmp(esmName, esm)) == 0
            disp(['the ', esm, ' ensemble of ', mdlName, ' didnt exist']);
            continue
        end

        specificNum = find(strcmp(esmName, 'r1i1p1f1') == 1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensemble member
        for esmNum = specificNum:specificNum%1:length(esmName)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% load and read
            esmPath = fullfile(mdlPath, esmName{esmNum, 1});
            % data path
            varsPath = fullfile(esmPath, level.process3{1}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/rawdata_regrid
            dvarsPath = fullfile(esmPath, level.process3{2}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
            dvarsTrendPath = fullfile(esmPath, level.process3{3}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly_trend
            kernelPath = fullfile(esmPath, level.process3{5}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/kernelsCal
            dradEffectPath = fullfile(esmPath, level.process3{6}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
            dnonLocalCldPath = fullfile(esmPath, level.process3{8}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/non_localCld/
            load([dradEffectPath, 'global_vars.mat'])% lat_f lon_f time plev_k readme
            ntime = length(time.date);

            % cal surface temperature
            load([varsPath, 'ts.mat'])% ts
            load([dvarsPath, 'dts.mat'])%dts
            dts = autoRegrid3(lon_k, lat_k, time.date, dts, lon_f, lat_f, time.date);
            ts = autoRegrid3(lon_k, lat_k, time.date, ts, lon_f, lat_f, time.date);

            % cal rhs(RHeating)
            load([dvarsPath, 'drlds.mat'])% surface_downwelling_longwave_flux_in_air
            load([dvarsPath, 'drsds.mat'])% surface_downwelling_shortwave_flux_in_air
            load([dvarsPath, 'drsus.mat'])% surface_upwelling_shortwave_flux_in_air
            dR_swnet = drsds - drsus; % sfc net shortwave flux
            drhs = drlds + dR_swnet; % equilibrium equation's RHS, nearly equal to sfc upward rad
            drhs = autoRegrid3(lon_k, lat_k, time.date, drhs, lon_f, lat_f, time.date);

            % cloud radRffect
            load([dradEffectPath, 'dR_cloud_toa.mat'])
            load([dradEffectPath, 'dR_cloud_sfc.mat'])

            % TotalEffect
            load([dradEffectPath, 'dradEfect_sfc_cld.mat'])% 'totalEffect', 'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', 'taEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect', 'mainEffect'
            load([dradEffectPath, 'real_dradEfect.mat'])% 'dR_allsky', 'l_rad', 's_rad', 'dR_clr', 'readme_realradEfect'
            load([dradEffectPath, 'dR_residual_cld_sfc.mat'])% dR_residual_cld_sfc
            load([dradEffectPath, 'dR_residual_cld_toa.mat'])% dR_residual_cld_toa
            dR_netSfc = squeeze(dR_allsky(:, :, :, 1));
            dR_netTOA = squeeze(dR_allsky(:, :, :, 2));
            drhsKern = albEffect + husEffect + taEffect + dR_cloud_sfc;

            if exmNum == 1
                % cut to 2000.03-2014-02
                timeStr = string(datestr(datenum(time.date), 'yyyy-mm'));
                cutStart = find(timeStr == '2000-03');
                cutEnd = find(timeStr == '2014-02');
                time = time.date(cutStart:cutEnd);
                ntime = length(time);
                drhs = drhs(:, :, cutStart:cutEnd);
                drhsKern = drhsKern(:, :, cutStart:cutEnd);
                tsEffect = tsEffect(:, :, cutStart:cutEnd);
                startMonth=3;
            else
                time = time.date;
                ntime = length(time);
                startMonth=1;

            end

            varUsed = zeros(nlonf, nlatf, ntime, 2);
            varUsed(:, :, :, 1) = -tsEffect;
            varUsed(:, :, :, 2) = drhs;
            varUsed(:, :, :, 3) = drhsKern;
            [outPutPath] = plotFig(varUsed,time,lat_f,startMonth);
            if exmNum == 1
                title_txt = {['Level:', num2str(toaSfc{2}), ', Era: 200003-201402'], ['Model:', level.model2{mdlNum} ', Ensemble: ', esmName{esmNum}], ['90N-90S land, global mean'],[]};
            else
                title_txt = {['Level:', num2str(toaSfc{2}), ', Era: ', level.time1{exmNum}(1:end - 1)], ['Model:', level.model2{mdlNum} ', Ensemble: ', esmName{esmNum}], ['90N-90S land, global mean'],[]};

            end
            title(title_txt, 'FontWeight', 'normal', 'Interpreter', 'none')
            
            % save figures
            outPutPath=fullfile(outPutPath,'Fig3_CMIP6_land_SFC_globalLandMeanTimeSeries');
            auto_mkdir(outPutPath)
            figName = ['Fig3_', level.model2{mdlNum}, '_', esmName{esmNum},'_',level.time1{exmNum}(1:end - 1),];
            figPath = [outPutPath, '/', figName, '.eps'];
            export_fig(gcf,figPath,'-r600','-cmyk')
            % saveas(gcf, figurePath)
            % save_png(figurePath)%high resolution
            % close gcf

        end

    end

end

eval(['cd ', nowpath]);
t = toc; disp(t)

function [outPutPath] = plotFig(varUsed, time, lat_f, startMonth)
    outPutPath = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/proj3_PaperFig/v0.0/');
    auto_mkdir(outPutPath)

    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(1);

    latRange = 90; % Latitude range

    varNames = {['d', '\itLW', '\rm_{up}'], ['d', '\itR_{Heating}'], ['d', '\itR_{Heating}'], 'ts RadEfect', 'q RadEfect', 'alb RadEfect'};
    yLabel = {'W m^{-2}', 'W m^{-2}', 'W m^{-2}', 'W m^{-2}', 'W m^{-2}', 'W m^{-2}'};

    sizeVarUsed = size(varUsed);
    sizeLon = sizeVarUsed(1); sizeLat = sizeVarUsed(2); sizeTime = sizeVarUsed(3); sizeVar = sizeVarUsed(4);
    varUsedYearly = zeros(sizeLon, sizeLat, sizeTime / 12, sizeVar);

    for varNum = 1:sizeVar
        varUsedYearly(:, :, :, varNum) = monthlyToYearly(varUsed(:, :, :, varNum));
    end

    % mask only -90-90 land.
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

    if startMonth ~= 1
        timeYearly = timeYearly(1:end - 1);
    end

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
    timeSer = [1980 1985 1990 1995 2000 2005 2010 2015 2020 2025 2035 2045 2055 2065 2075 2085 2095 2105];

    if startMonth ~= 1
        timeSer = [1980 1985 1990 1995 2000 2002 2004 2006 2008 2010 2012 2014 2018 2025 2035 2045 2055 2065 2075 2085 2095 2105];
    end

    char_timeSer = cellstr(string(timeSer));
    % fig set
    f_row = 1; f_col = 1; % 设置画图的行列
    set(0, 'DefaultFigureVisible', 'on')
    ss = get(0, 'ScreenSize');
    coef_amplify = 1.5;
    % h = figure('Position', [ss(4) / 2 - 100 ss(3) / 35 ss(3) / 5 * f_col * coef_amplify (ss(4) - 80) / 5 * f_row * coef_amplify]);
    % h = figure('Position', [-1051 527 672 439]);
    h = figure('Position', [-1075 435 889 569]);
    % clf reset;
    set(h, 'Color', [1 1 1]);
    % zero line
    zeroLine = zeros(length(timeYearly), 1);

    varUsedTemp1 = squeeze(varUsedYearly_weightMean(:, 3));
    varUsedTemp2 = squeeze(varUsedYearly_weightMean(:, 2));
    varUsedTemp3 = squeeze(varUsedYearly_weightMean(:, 1));

    lineWdth = 2;
    p1=plot(timeYearly, varUsedTemp1', 'color', '#F08212', 'LineWidth', lineWdth);% 橘黄
    hold on
    p2=plot(timeYearly, varUsedTemp2', 'color', '#C2162E', 'LineWidth', lineWdth);% 深红
    hold on
    p3=plot(timeYearly, varUsedTemp3', 'color', '#4450A1', 'LineWidth', lineWdth);% 蓝色
    hold on
    plot(timeYearly, zeroLine', 'k--', 'LineWidth', lineWdth)
    hold on;
    % set x axes
    xticks(timeSer)
    datetick('x', 'yyyy', 'keepticks'); % 用日期的格式显示横坐标
    xlim([timeYearly(1) timeYearly(end)])
    xlabel('time line')
    xticklabels(char_timeSer)

    ylabel(yLabel{2})
    % set y axes
    ymax = 5;
    yticks([-5 -4 -3 -2 -1 0 1 2 3 4 5])
    yticklabels({'-5','-4','-3','-2','-1','0','1','2','3','4','5'})

    ylim([-ymax ymax])

    ax = gca;
    ax.FontName='Microsoft YaHei';% Microsoft YaHei 'Time New Roman'
    ax.XMinorTick = 'off'; ax.YMinorTick = 'off'; % 开启次刻度线
    ax.XAxis.MinorTickValues = (1980:1:2105);
    ax.TickLength = [0.015 0.01]; %刻度线长度      set(gca,'ticklength', [0.02 0.01]);
    ax.XColor = 'k'; ax.YColor = 'k'; % 设置刻度线颜色
    ax.FontSize = 18;
    ax.LineWidth = 1.5;


    % 去掉上边框和右边框的刻度
    box off
    xtick = get(gca, 'XTick');
    ytick = get(gca, 'YTick');
    line([xtick(1), xtick(end)], [ytick(end) ytick(end)], 'Color', 'k', 'LineWidth', 1.5)
    line([timeYearly(end), timeYearly(end)], [ytick(end) ytick(1)], 'Color', 'k', 'LineWidth', 1.5)

    % ylabel({Level{i}, ' Rad Anomaly (Wm^{-2})'}, 'Fontsize', 10)
    lgd = legend([p1 p2 p3],['d', '\itR', '\rm_{Heat}(kernel calc) cc=', num2str(cc_weightMean{3})], ['d', '\itR', '\rm_{Heat} cc=', num2str(cc_weightMean{2})], ['d', '\itLW', '\rm_{up}'], 'Fontsize', 16, 'Location', 'northwest', 'NumColumns', 1); %'FontWeight', 'bold',
    legend('boxoff')%删除图例背景和轮廓
    % lgd_inf = get(lgd);
    % text(lgd_inf.Position(1) + 0.15, lgd_inf.Position(2) + 0.235, ['cc= ', num2str(cc_weightMean{2})], 'Fontsize', 10, 'FontWeight', 'bold', 'Interpreter', 'none', 'Units', 'normalized')
    hold on
end

