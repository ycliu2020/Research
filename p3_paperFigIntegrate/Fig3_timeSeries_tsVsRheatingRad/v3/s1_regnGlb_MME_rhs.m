%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-08-31 17:00:15
% LastEditTime : 2021-06-23 21:53:43
% LastEditors  : Please set LastEditors
% Description  : MME result of 时间序列图
% FilePath     : /code/p3_paperFigIntegrate/Fig3_timeSeries_tsVsRheatingRad/v3/s1_regnGlb_MME_rhs.m
% note : 这里的rhs 实际上是rhs Kern, 这是为了和
%%---------------------------------------------------------
clc; clear; tic;
nowpath = pwd;
% load mask map
run '/home/liuyc/lib/tools/matlab/myTools/autoScript/preLoadVar.m'

timeType = 'monthly'; % monthly or yearly
MMEType='MME1';
MMENum=str2num(MMEType(end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment
for exmNum = 1:2 %1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    % ERA5 number
    exmNum_ERA5 = exmNum + 3;
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % MMEPath
    MMEPath = fullfile(level.path_MME, MMEType, level.time1{exmNum}); %/data1/liuyincheng/CMIP6-process/2000-2014/
    dvarsPath = fullfile(MMEPath, level.process3{2});
    dradEffectPath = fullfile(MMEPath, level.process3{6});
    figDataPath = fullfile(MMEPath, 'FigData/');
    %% load and read
    load([dradEffectPath, 'global_vars.mat']) % 'lon_f', 'lat_f', 'timeEssmble', 'time', 'plev_k', 'readme', 'timeseries', 'MME_Models'
    ntime = length(time.date);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load MME and save one var
    load([figDataPath, 'regionalTsRHeating.mat']) %  'drhsMask', 'drhsKernMask', 'dlwUPKernMask', 'dlwUPRawDataMask', 'dhfssMask', 'dhflsMask' .world .CHNeast .EURwest .USeast
    % % cal the lwUp caled by kernel
    % dlwUPKernMask.world = -dR_tsMask.world;
    % dlwUPKernMask.CHNeast = -dR_tsMask.CHNeast;
    % dlwUPKernMask.USeast = -dR_tsMask.USeast;
    % dlwUPKernMask.EURwest = -dR_tsMask.EURwest;

    regionCode = {'world', 'CHNeast', 'EURwest', 'USeast'}; %'CHNeast', 'EURwest', 'USeast'
    varibalCode = {'drhsMask', 'dlwUPKernMask'};

    varUsed = zeros(nlonf, nlatf, ntime, length(regionCode), length(varibalCode), 2); % dim:[lon, lat, time, region, varibal, dataset]
    sizeVarUsed = size(varUsed);
    sizeLon = sizeVarUsed(1); sizeLat = sizeVarUsed(2); sizeTime = sizeVarUsed(3);
    sizeRegion = sizeVarUsed(4); sizeVar = sizeVarUsed(5); sizeData = sizeVarUsed(6);

    for regionNum = 1:sizeRegion

        for varNum = 1:sizeVar
            varUsed(:, :, :, varNum, regionNum, 2) = eval([varibalCode{varNum}, '.', regionCode{regionNum}]);
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ERA5 Path
    [~, ~, tLin, vars] = obsParameters('ERA5');
    ERA5Path = fullfile('/data2/liuyincheng/Observe-process', tLin.time{exmNum_ERA5}, 'ERA5/');
    ERA5figDataPath = fullfile(ERA5Path, 'FigData/');
    % load ERA5 and add to one var
    load([ERA5figDataPath, 'regionalTsRHeating.mat']) % 'drhsMask', 'drhsKernMask', 'dlwUPKernMask', 'dlwUPRawDataMask', 'dhfssMask', 'dhflsMask'  .world .CHNeast  .EURwest .USeast

    for regionNum = 1:sizeRegion

        for varNum = 1:sizeVar
            varUsed(:, :, :, varNum, regionNum, 1) = eval([varibalCode{varNum}, '.', regionCode{regionNum}]);
        end

    end

    % add Res
    varUsed(:, :, :, 3, :, :) = varUsed(:, :, :, 1, :, :) - varUsed(:, :, :, 2, :, :);
    sizeVar = sizeVar + 1;

    % global annual mean
    varUsedYearly_weightMean = glbAnlMean(varUsed, sizeLon, sizeLat, sizeTime, sizeVar, sizeRegion, sizeData);
    % set x axis
    startMonth = 1;
    time.vec = datevec(time.date);
    timeYearly = unique(time.vec(:, 1));

    if startMonth ~= 1
        timeYearly = timeYearly(1:end - 1);
    end

    % global mean
    varUsed_weightMean = glbMean(varUsed, sizeLon, sizeLat, sizeTime, sizeVar, sizeRegion, sizeData);

    if strcmp(timeType, 'monthly') == 1
        varUsed_plot = varUsed_weightMean;
        time_plot = time.date;
    elseif strcmp(timeType, 'yearly') == 1
        varUsed_plot = varUsedYealy_weightMean;
        time_plot = timeYearly;
    end

    % cal cc of areaMeanLatWeight
    [cc_weightMean, pp_weightMean] = calCC(varUsed_plot, sizeVar, sizeRegion, sizeData);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plot Part
    plotPart(exmNum, timeType, time_plot, varUsed_plot, cc_weightMean, sizeLon, sizeLat, sizeTime, sizeVar, sizeRegion, sizeData);
    title_txt = {['Level: ', num2str(toaSfc{2})], ['Data: ', level.time1{exmNum}(6:end - 11), ' ', level.time1{exmNum}(end - 9:end - 1),', MME Type: ', MMEType]};
    sgtitle(title_txt, 'FontWeight', 'normal') %,'Units', 'normalized', 'position', [0.5, 0.98]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% save Part
    mPath.uniOutput = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/figOfSketch/TimeSeries/RheatingVsLW',timeType); %['dRTs_', lower(mlabels.level)],
    mPath.Output = fullfile(mPath.uniOutput);
    auto_mkdir(mPath.Output)
    figName = [MMEType,'VsERA5_regnGlb_rhs_', tLin.time{exmNum_ERA5}];

    figPath = [mPath.Output, '/', figName, '.eps'];
    export_fig(gcf, figPath, '-r600', '-cmyk')
    % saveas(gcf, figurePath)
    % save_png(figurePath)%high resolution
    % close gcf

end

eval(['cd ', nowpath]);
t = toc; disp(t)

function varUsedYearly_weightMean = glbAnlMean(varUsed, sizeLon, sizeLat, sizeTime, sizeVar, sizeRegion, sizeData)
    run '/home/liuyc/lib/tools/matlab/myTools/autoScript/preLoadVar.m'

    varUsedYearly = zeros(sizeLon, sizeLat, sizeTime / 12, sizeVar, sizeRegion, sizeData);

    for regionNum = 1:sizeRegion

        for varNum = 1:sizeVar

            for dataNum = 1:sizeData
                varUsedYearly(:, :, :, varNum, regionNum, dataNum) = monthlyToYearly(varUsed(:, :, :, varNum, regionNum, dataNum));
            end

        end

    end

    sizeVarUsedYearly = size(varUsedYearly);
    varUsedYearly_weightMean = zeros(sizeVarUsedYearly(3), sizeVar, sizeRegion, sizeData);

    for regionNum = 1:sizeRegion

        for varNum = 1:sizeVar

            for timeNum = 1:sizeVarUsedYearly(3)

                for dataNum = 1:sizeData
                    varUsedYearly_weightMean(timeNum, varNum, regionNum, dataNum) = areaMeanLatWeight(varUsedYearly(:, :, timeNum, varNum, regionNum, dataNum), lat_f);
                end

            end

        end

    end

end

function varUsed_weightMean = glbMean(varUsed, sizeLon, sizeLat, sizeTime, sizeVar, sizeRegion, sizeData)
    run '/home/liuyc/lib/tools/matlab/myTools/autoScript/preLoadVar.m'
    varUsed_weightMean = zeros(sizeTime, sizeVar, sizeRegion, sizeData);

    for regionNum = 1:sizeRegion

        for varNum = 1:sizeVar

            for timeNum = 1:sizeTime

                for dataNum = 1:sizeData
                    varUsed_weightMean(timeNum, varNum, regionNum, dataNum) = areaMeanLatWeight(varUsed(:, :, timeNum, varNum, regionNum, dataNum), lat_f);
                end

            end

        end

    end

end

function [cc_weightMean, pp_weightMean] = calCC(varUsedYearly_weightMean, sizeVar, sizeRegion, sizeData)
    cc_weightMean = cell(sizeVar, sizeRegion, sizeData);
    pp_weightMean = cc_weightMean;

    for dataNum = 1:sizeData

        for regionNum = 1:sizeRegion

            for varNum = 1:sizeVar
                [cc0, pp0] = corrcoef(varUsedYearly_weightMean(:, 2, regionNum, dataNum), varUsedYearly_weightMean(:, varNum, regionNum, dataNum), 'Rows', 'complete');
                cc_weightMean{varNum, regionNum, dataNum} = roundn(cc0(1, 2), -2); %保留俩位小数
                pp_weightMean{varNum, regionNum, dataNum} = pp0(1, 2); % confidence interval
            end

        end

    end

end

function [] = plotPart(exmNum, timeType, time_plot, varUsed_plot, cc_weightMean, sizeLon, sizeLat, sizeTime, sizeVar, sizeRegion, sizeData)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plot Part
    radName = {['\Delta', '\itRH,', '\rm cc='], ['-', '\itK', '\rm_{Ts}', '\cdot\Delta', '\itT', '\rm_s'], ['-\Delta\itTHF','\rm+\itRes']};
    unitName = {'W m^{-2}'};
    dataName = {'ERA5', 'AMIP'};
    regionName = {['Global'], ['Eastern CHN'], ['Western EUR'], ['Eastern US']};
    set(0, 'defaultfigurecolor', 'w'); %设置画布底色为白色
    set(0, 'DefaultFigureVisible', 'on')

    ss = get(0, 'ScreenSize');
    h = figure('Position', [-1071 89 934 918]); %[离左边缘 离下边缘 自身宽 自身高][ss(4) / 2 - 100 ss(3) / 35 ss(3) / 5 * 2.5, (ss(4) - 80) / 5 * 4]

    zeroLine = zeros(length(time_plot), 1);
    subNum = 1;

    for regionNum = 1:4

        for dataNum = 1:2 % ERA5 or CMIP6 Data

            subplot_yc(4, 2, regionNum, dataNum);
            hold on
            lineWdth = 1.5;
            % plot(time_plot, squeeze(varUsed_plot(:, 1, regionNum, dataNum)), 'color', '#F08212', 'LineWidth', lineWdth) % 计算值
            p3 = plot(time_plot, squeeze(varUsed_plot(:, 3, regionNum, dataNum)), 'color', '#A8A8A8', 'LineWidth', 1); % res
            p1 = plot(time_plot, squeeze(varUsed_plot(:, 1, regionNum, dataNum)), 'color', '#BC144A', 'LineWidth', lineWdth); % 真实值值
            p2 = plot(time_plot, squeeze(varUsed_plot(:, 2, regionNum, dataNum)), 'color', '#30388B', 'LineWidth', lineWdth); % 长波向上辐射
            plot(time_plot, zeroLine', 'k--', 'LineWidth', 1.5)

            % x axis
            if exmNum == 1
                timeSer = [2002 2006 2010 2014];
            else
                timeSer = [1982 1990 1998 2006 2014 2020 2025 2035 2045 2055 2065 2075 2085 2095 2105];
            end

            xticks(timeSer)
            char_timeSer = cellstr(string(timeSer));

            if strcmp(timeType, 'monthly')
                timeNum = datenum(timeSer, 6, 15);
                datetick('x', 'yyyy', 'keepticks'); % 用日期的格式显示横坐标
                xticks(timeNum)
            end

            xticklabels([])

            if regionNum == 4
                xticklabels(char_timeSer)
            end

            xlim([time_plot(1) time_plot(end)])

            % y axis
            xlim([time_plot(1) time_plot(end)])
            ymax = 10;
            ylim([-ymax ymax])

            if dataNum == 2 && exmNum == 1 && strcmp(timeType, 'yearly')
                ylim([-5 ymax])
            elseif strcmp(timeType, 'monthly')
                ylim([-20 20])

                if dataNum == 1
                    ylim([-40 40])
                end

            end

            if regionNum == 1
                ylim([-5 5])
                if exmNum == 2
                    ylim([-8 8])
                end

            end

            % y label
            if dataNum == 1
                yLabel = [regionName{regionNum}; unitName];
                ylabel(yLabel, 'FontName', 'Arial', 'Fontsize', 17)
            end

            % 上标题
            if regionNum == 1
                text(0.45, 1.15, dataName{dataNum}, 'FontName', 'Arial', 'Fontsize', 17, 'units', 'normalized');
            end

            % 添加序号
            if regionNum == 4
                text(0.475, -0.2, ['(', char(96 + subNum), ')'], 'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize', 15, 'units', 'normalized');
            else
                text(0.475, -0.1, ['(', char(96 + subNum), ')'], 'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize', 15, 'units', 'normalized');
            end

            subNum = subNum + 1;

            % legend
            gcaGet = get(gca);
            axPosit = gcaGet.Position;
            axPosit = axPosit + [0.0 0.14 -0.04 -0.16]; %(左边缘, 下边缘, 长度, 宽度)
            legend([p2, p1, p3], {radName{2}, [radName{1}, num2str(cc_weightMean{1, regionNum, dataNum})], radName{3}}, 'FontName', 'Arial', 'Fontsize', 9, 'Location', 'none', 'position', axPosit, 'NumColumns', 2); %'Location', 'best''FontWeight', 'bold',
            legend('boxoff') %删除图例背景和轮廓

            % ax set
            ax = gca;
            ax.FontName = 'Arial'; % Arial 'Time New Roman'
            ax.XMinorTick = 'on'; ax.YMinorTick = 'on'; % 开启次刻度线
            % ax.XAxis.MinorTickValues = timeNumFull;
            ax.TickLength = [0.03 0.02]; %刻度线长度      set(gca,'ticklength', [0.02 0.01]);
            ax.XColor = 'k'; ax.YColor = 'k'; % 设置刻度线颜色
            ax.FontSize = 12;
            ax.LineWidth = 1.5;

        end

    end

end

% % detrend
% for dataNum = 1:sizeData

%     for regionNum = 1:sizeRegion

%         for varNum = 1:sizeVar
%             varUsedYearly_weightMean(:, varNum, regionNum, dataNum) = detrend(varUsedYearly_weightMean(:, varNum, regionNum, dataNum));
%         end

%     end

% end
