%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-08-31 17:00:15
% LastEditTime : 2021-04-09 20:24:49
% LastEditors  : Please set LastEditors
% Description  : MME result of 时间序列图
% FilePath     : /code/p3_paperFigIntegrate/Fig3_timeSeries_tsVsRheatingRad/s1_plotRegional_containWorld.m
% note : 统一用startmonth=3 开始计算
%%---------------------------------------------------------
clc; clear; tic;
nowpath = pwd;
% load mask map
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat') % load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat') % load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')

toaSfc = {'toa', 'sfc'};
lon_k = 0:2.5:357.5; nlonk = length(lon_k); % kernel lat lon
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read data
% the number of MME experiment
MMECode = 'MME1';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment
for exmNum_MME = 1:2 %1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    % ERA5 number
    exmNum_ERA5 = exmNum_MME + 3;
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum_MME);
    % MME Path
    MMEPath = ['/data1/liuyincheng/CMIP6-process/z_ensembleMean/MME1/', level.time1{exmNum_MME}]; %/data1/liuyincheng/CMIP6-process/z_ensembleMean/MME1/CMIP/amip_1980-2014
    dvarsPath = fullfile(MMEPath, level.process3{2});
    dradEffectPath = fullfile(MMEPath, level.process3{6});
    figDataPath = fullfile(MMEPath, 'FigData/');

    % ERA5 Path
    [readme, level, tLin, vars] = obsParameters('ERA5');
    ERA5Path = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum_ERA5}, 'ERA5/');
    ERA5figDataPath = fullfile(ERA5Path, 'FigData/');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% load and read
    load([dradEffectPath, 'global_vars.mat']) % 'lon_f', 'lat_f', 'timeEssmble', 'time', 'plev_k', 'readme', 'timeseries', 'MME_Models'
    ntime = length(time.date);

    % load MME and save one var
    load([figDataPath, 'regionalTsRHeating.mat']) % dlwUPMask, drhsKernMask, drhsMask  .world .CHNeast  .EURwest .USeast

    regionCode = {'world', 'CHNeast', 'EURwest', 'USeast'}; %'CHNeast', 'EURwest', 'USeast'
    varibalCode = {'drhsKernMask', 'drhsMask', 'dlwUPMask'};
    varUsed = zeros(nlonf, nlatf, ntime, length(regionCode), length(varibalCode), 2); % dim:[lon, lat, time, region, varibal, dataset]
    sizeVarUsed = size(varUsed);
    sizeLon = sizeVarUsed(1); sizeLat = sizeVarUsed(2); sizeTime = sizeVarUsed(3);
    sizeRegion = sizeVarUsed(4); sizeVar = sizeVarUsed(5); sizeData = sizeVarUsed(6);

    for regionNum = 1:sizeRegion

        for varNum = 1:sizeVar
            varUsed(:, :, :, varNum, regionNum, 2) = eval([varibalCode{varNum}, '.', regionCode{regionNum}]);
        end

    end

    % load ERA5 and add to one var
    load([ERA5figDataPath, 'regionalTsRHeating.mat']) % dlwUPMask, drhsKernMask, drhsMask  .world .CHNeast  .EURwest .USeast

    for regionNum = 1:sizeRegion

        for varNum = 1:sizeVar
            varUsed(:, :, :, varNum, regionNum, 1) = eval([varibalCode{varNum}, '.', regionCode{regionNum}]);
        end

    end

    %% compute areaMeanLatWeight
    % method: 先求年平均,在进行纬向加权平均
    varUsedYearly = zeros(sizeLon, sizeLat, sizeTime / 12, sizeVar, sizeRegion, sizeData);

    for dataNum = 1:sizeData

        for regionNum = 1:sizeRegion

            for varNum = 1:sizeVar
                varUsedYearly(:, :, :, varNum, regionNum, dataNum) = monthlyToYearly(varUsed(:, :, :, varNum, regionNum, dataNum));
            end

        end

    end

    sizeVarUsedYearly = size(varUsedYearly);
    varUsedYearly_weightMean = zeros(sizeVarUsedYearly(3), sizeVar, sizeRegion, sizeData);

    for dataNum = 1:sizeData

        for regionNum = 1:sizeRegion

            for varNum = 1:sizeVar

                for timeNum = 1:sizeVarUsedYearly(3)
                    varUsedYearly_weightMean(timeNum, varNum, regionNum, dataNum) = areaMeanLatWeight(varUsedYearly(:, :, timeNum, varNum, regionNum, dataNum), lat_f);
                end

            end

        end

    end

    % detrend
    for dataNum = 1:sizeData

        for regionNum = 1:sizeRegion

            for varNum = 1:sizeVar
                varUsedYearly_weightMean(:, varNum, regionNum, dataNum) = detrend(varUsedYearly_weightMean(:, varNum, regionNum, dataNum));
            end

        end

    end

    % cal cc of areaMeanLatWeight
    cc_weightMean = cell(sizeVar, sizeRegion, sizeData);
    pp_weightMean = cc_weightMean;

    for dataNum = 1:sizeData

        for regionNum = 1:sizeRegion

            for varNum = 1:sizeVar
                [cc0, pp0] = corrcoef(varUsedYearly_weightMean(:, 3, regionNum, dataNum), varUsedYearly_weightMean(:, varNum, regionNum, dataNum), 'Rows', 'complete');
                cc_weightMean{varNum, regionNum, dataNum} = roundn(cc0(1, 2), -2); %保留俩位小数
                pp_weightMean{varNum, regionNum, dataNum} = pp0(1, 2); % confidence interval
            end

        end

    end

    % x axis: time
    startMonth = 1;
    time.vec = datevec(time.date);
    timeYearly = unique(time.vec(:, 1));

    if startMonth ~= 1
        timeYearly = timeYearly(1:end - 1);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plot Part
    radName = {['d', '\itR', '\rm_{heating,Kern} cc='], ['d', '\itR', '\rm_{heating} cc='], ['d', '\itLW', '\rm_{up} ']};
    unitName = {'W m^{-2}'};
    dataName = {'ERA5', 'AMIP'};
    regionName = {['World'], ['Eastern CHN'], ['Western EUR'], ['Eastern US']};
    set(0, 'defaultfigurecolor', 'w'); %设置画布底色为白色
    set(0, 'DefaultFigureVisible', 'on')

    ss = get(0, 'ScreenSize');
    h = figure('Position', [-1071 89 934 918]); %[离左边缘 离下边缘 自身宽 自身高][ss(4) / 2 - 100 ss(3) / 35 ss(3) / 5 * 2.5, (ss(4) - 80) / 5 * 4]

    zeroLine = zeros(length(timeYearly), 1);
    subNum = 1;

    for regionNum = 1:4

        for dataNum = 1:2 % ERA5 or CMIP6 Data

            subplot_yc(4, 2, regionNum, dataNum);
            hold on
            lineWdth = 2.5;
            plot(timeYearly, squeeze(varUsedYearly_weightMean(:, 1, regionNum, dataNum)), 'color', '#F08212', 'LineWidth', lineWdth) % 计算值
            plot(timeYearly, squeeze(varUsedYearly_weightMean(:, 2, regionNum, dataNum)), 'color', '#BC144A', 'LineWidth', lineWdth) % 真实值值
            plot(timeYearly, squeeze(varUsedYearly_weightMean(:, 3, regionNum, dataNum)), 'color', '#30388B', 'LineWidth', lineWdth) % 长波向上辐射
            plot(timeYearly, zeroLine', 'k--', 'LineWidth', 1.5)

            % x axis
            if exmNum_MME == 1
                timeSer = [2002 2006 2010 2014];
            else
                timeSer = [1980 1985 1990 1995 2000 2005 2010 2015 2020 2025 2035 2045 2055 2065 2075 2085 2095 2105];

            end

            char_timeSer = cellstr(string(timeSer));

            xticks(timeSer)
            xticklabels([])

            if regionNum == 4
                xticklabels(char_timeSer)
            end

            % y axis
            xlim([timeYearly(1) timeYearly(end)])
            ymax = 10;
            ylim([-ymax ymax])

            if dataNum == 2 && exmNum_MME == 1
                ylim([-5 ymax])
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
            legend({[radName{1}, num2str(cc_weightMean{1, regionNum, dataNum})], [radName{2}, num2str(cc_weightMean{2, regionNum, dataNum})], radName{3}}, 'FontName', 'Arial', 'Fontsize', 9, 'Location', 'none', 'position', axPosit, 'NumColumns', 2); %'Location', 'best''FontWeight', 'bold',
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% save Part
    mPath.uniOutput = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/proj3_PaperFig/v0.3/Fig3_regionalLandMeanTimeSeries/'); %['dRTs_', lower(mlabels.level)],
    mPath.Output = fullfile(mPath.uniOutput);
    auto_mkdir(mPath.Output)
    figName = ['Fig3_', 'AMIP-MME1AndERA5_containWorld', '_', tLin.time{exmNum_ERA5}];

    figPath = [mPath.Output, '/', figName, '.eps'];
    export_fig(gcf, figPath, '-r600', '-cmyk')
    % saveas(gcf, figurePath)
    % save_png(figurePath)%high resolution
    % close gcf

end

eval(['cd ', nowpath]);
t = toc; disp(t)
