%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-08-31 17:00:15
% LastEditTime : 2021-05-06 14:33:42
% LastEditors  : Please set LastEditors
% Description  : 同时画时间序列和相关性分布图
% FilePath     : /code/p3_paperFigIntegrate/Fig6_7_timeSeriANS_tsVsVarsRad/s1_MME_tsVsVarsRad_plotRegional.m
%
%%---------------------------------------------------------
clc; clear; tic;
nowpath = pwd;
run '/home/liuyc/lib/tools/matlab/myTools/autoScript/preLoadVar.m'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read data
% the number of MME experiment
MMECode = 'MME1';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment
detrendCode = 1;

for exmNum_MME = 4:4 %1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    %CAMS-CSM1-0 didn't have sfc clear sky radiation, delete it
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum_MME);
    % MMEPath
    MMEPath = ['/data1/liuyincheng/CMIP6-process/z_ensembleMean/MME1/', level.time1{exmNum_MME}]; %/data1/liuyincheng/CMIP6-process/2000-2014/
    dvarsPath = fullfile(MMEPath, level.process3{2});
    dradEffectPath = fullfile(MMEPath, level.process3{6});
    figDataPath = fullfile(MMEPath, 'FigData/');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% load and read
    load([dradEffectPath, 'global_vars.mat']) % 'lon_f', 'lat_f', 'timeEssmble', 'time', 'plev_k', 'readme', 'timeseries', 'MME_Models'
    ntime = length(time.date);

    % load MME and save one var
    load([figDataPath, 'regionalVarsRad_sfc_cld.mat']) % 'dR_cloudMask', 'dR_residualMask', 'dR_albMask', 'dR_husMask', 'dR_taMask', 'dR_tsMask'
    load([figDataPath, 'regionalTsRHeating.mat']) % 'drhsMask', 'drhsKernMask', 'dlwUPMask', 'dlwUPRawDataMask', 'dhfssMask', 'dhflsMask' .world .CHNeast .EURwest .USeast
    % dlwUPMask.world = -dR_tsMask.world;
    % dlwUPMask.CHNeast = -dR_tsMask.CHNeast;
    % dlwUPMask.USeast = -dR_tsMask.USeast;
    % dlwUPMask.EURwest = -dR_tsMask.EURwest;
    regionCode = {'world', 'CHNeast', 'EURwest', 'USeast'}; %'CHNeast', 'EURwest', 'USeast'
    varibalCode = {'dlwUPMask', 'dR_husMask', 'dR_taMask', 'dR_albMask', 'dR_residualMask', 'dR_cloudMask', 'dhfssMask', 'dhflsMask'};

    varUsed = zeros(nlonf, nlatf, ntime, length(varibalCode), length(regionCode));
    sizeVarUsed = size(varUsed);
    sizeLon = sizeVarUsed(1); sizeLat = sizeVarUsed(2); sizeTime = sizeVarUsed(3);
    sizeVar = sizeVarUsed(4); sizeRegion = sizeVarUsed(5);

    for regionNum = 1:sizeRegion

        for varNum = 1:sizeVar
            varUsed(:, :, :, varNum, regionNum) = eval([varibalCode{varNum}, '.', regionCode{regionNum}]);
        end

    end

    % varUsed(:,:,:,5,:)=varUsed(:,:,:,5,:)-varUsed(:,:,:,7,:)-varUsed(:,:,:,8,:);% Rnet-Rvars-LH-SH
    varUsed(:, :, :, 5, :) = varUsed(:, :, :, 1, :) - varUsed(:, :, :, 2, :) - varUsed(:, :, :, 3, :) - varUsed(:, :, :, 4, :) - varUsed(:, :, :, 6, :); % Rts-Rvars
    % varUsed(:,:,:,5,:)=varUsed(:,:,:,1,:)-varUsed(:,:,:,2,:)-varUsed(:,:,:,3,:)-varUsed(:,:,:,4,:)-varUsed(:,:,:,6,:)-varUsed(:,:,:,7,:)-varUsed(:,:,:,8,:);% Rts-Rvars-LH-SH
    sizeVar = 6;
    %% compute areaMeanLatWeight
    % method: 先求年平均,在进行纬向加权平均
    varUsedYearly = zeros(sizeLon, sizeLat, sizeTime / 12, sizeVar, sizeRegion);

    for regionNum = 1:sizeRegion

        for varNum = 1:sizeVar
            varUsedYearly(:, :, :, varNum, regionNum) = monthlyToYearly(varUsed(:, :, :, varNum, regionNum));
        end

    end

    sizeVarUsedYearly = size(varUsedYearly);
    varUsedYearly_weightMean = zeros(sizeVarUsedYearly(3), sizeVar, sizeRegion);

    for regionNum = 1:sizeRegion

        for varNum = 1:sizeVar

            for timeNum = 1:sizeVarUsedYearly(3)
                varUsedYearly_weightMean(timeNum, varNum, regionNum) = areaMeanLatWeight(varUsedYearly(:, :, timeNum, varNum, regionNum), lat_f);
            end

        end

    end

    % detrend
    varUsedYearly_weightMean2 = varUsedYearly_weightMean;

    for regionNum = 1:sizeRegion

        for varNum = 1:sizeVar
            varUsedYearly_weightMean2(:, varNum, regionNum, 2) = detrend(varUsedYearly_weightMean(:, varNum, regionNum));
        end

    end

    % cal cc of areaMeanLatWeight (有必要再算)
    sizeWeightMean = size(varUsedYearly_weightMean);
    cc_weightMean = cell(1, sizeWeightMean(2));
    pp_weightMean = cell(1, sizeWeightMean(2));

    for regionNum = 1:sizeRegion

        [cc0, pp0] = corrcoef(varUsedYearly_weightMean2(:, 5, regionNum, 2), varUsedYearly_weightMean2(:, 6, regionNum, 2), 'Rows', 'complete');
        cc_cld(regionNum) = roundn(cc0(1, 2), -2); %保留俩位小数
    end

    % for varNum = 1:sizeWeightMean(2)
    %     [cc0, pp0] = corrcoef(varUsedYearly_weightMean(:, 1), varUsedYearly_weightMean(:, varNum), 'Rows', 'complete');
    %     cc_weightMean{varNum} = roundn(cc0(1, 2), -2); %保留俩位小数
    %     pp_weightMean{varNum} = pp0(1, 2); % confidence interval
    % end

    % x axis: time
    startMonth = 1;
    time.vec = datevec(time.date);
    timeYearly = unique(time.vec(:, 1));

    if startMonth ~= 1
        timeYearly = timeYearly(1:end - 1);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plot Part
    varNames = {['d', '\itLW', '\rm_{up}'], ['d', '\itR', '\rm_{WV}'], ['d', '\itR', '\rm_{Ta}'], ['d', '\itR', '\rm_{Alb}'], ['\rm{Residual}'], ['d', '\itR', '\rm_{cld}']};
    unitName = {'W m^{-2}'};
    regionName = {['Global land'], ['Eastern CHN'], ['Western EUR'], ['Eastern US']};

    detrendSwich = {'no detrend', 'Detrend'};

    set(0, 'defaultfigurecolor', 'w'); %设置画布底色为白色
    set(0, 'DefaultFigureVisible', 'on')
    ss = get(0, 'ScreenSize');
    h = figure('Position', [10 18 1072 719]); %[离左边缘 离下边缘 自身宽 自身高][ss(4) / 2 - 100 ss(3) / 35 ss(3) / 5 * 2.5, (ss(4) - 80) / 5 * 4]

    varColor = {'#F08212', '#4450A1', '#67C1EE', '#90C64D', '#000000', '#C2162E'}; %橘色 深蓝 浅蓝 浅绿 黑色 红色

    % time series
    timeSer = [1985 1990 1995 2000 2005 2010 2015 2025 2035 2045 2055 2065 2075 2085 2095 2105];

    if exmNum_MME == 1
        timeSer = [2000 2002 2006 2010 2014];
    elseif exmNum_MME == 2
        timeSer = [1982 1990 1998 2006 2014 2020 2025 2035 2045 2055 2065 2075 2085 2095 2105];
    end

    char_timeSer = cellstr(string(timeSer));

    zeroLine = zeros(length(timeYearly), 1);
    locDemo = randn(2);

    for regionNum = 1:4 % 'world', 'CHNeast', 'EURwest', 'USeast'
        varUsedYearly_weightMean = varUsedYearly_weightMean2(:, :, :, detrendCode);
        [figCol, figRow] = ind2sub(size(locDemo), regionNum);
        subplot_yc(2, 2, figRow, figCol);
        hold on
        lineWdth = 2.5;

        for varNum = 1:sizeVar
            varUsedTemp = squeeze(varUsedYearly_weightMean(:, varNum, regionNum));
            plot(timeYearly, varUsedTemp', 'color', varColor{varNum}, 'LineWidth', lineWdth)
            hold on
        end

        % zero line
        zeroLine = zeros(length(timeYearly), 1);
        plot(timeYearly, zeroLine', 'k--')
        hold on

        % set x axes
        xticks(timeSer)

        if exmNum_MME == 1
            char_timeSer{1} = [];
        end

        xticklabels([])

        if regionNum == 3 || regionNum == 4
            xticklabels(char_timeSer)
        end

        xlim([timeYearly(1) timeYearly(end)])

        % y axis
        ymax = 10;

        if exmNum_MME == 1
            ymax = 4;

            if regionNum == 1
                ymax = 1.5;
            end

        elseif exmNum_MME == 2
            ymax = 5;
        elseif exmNum_MME == 4 && detrendCode == 2
            ymax = 4;

            if regionNum == 1
                ymax = 2;
            end

        end

        ylim([-ymax ymax])

        % y label
        if regionNum == 1 || regionNum == 3
            yLabel = unitName;
            ylabel(yLabel, 'FontName', 'Arial', 'Fontsize', 17)
        end

        % 上标题
        text(0.34, 1.08, regionName{regionNum}, 'FontName', 'Arial', 'Fontsize', 17, 'units', 'normalized');

        % 添加序号

        if regionNum == 3 || regionNum == 4
            text(0.475, -0.17, ['(', char(96 + regionNum), ')'], 'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize', 15, 'units', 'normalized');
        else
            text(0.475, -0.1, ['(', char(96 + regionNum), ')'], 'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize', 15, 'units', 'normalized');
        end

        % legend (单独放在外面)

        % 去掉上边框和右边框的刻度
        box off
        xtick = get(gca, 'XTick');
        ytick = get(gca, 'YTick');
        line([xtick(1), xtick(end)], [ytick(end) ytick(end)], 'Color', 'k', 'LineWidth', 1.5)
        line([timeYearly(end), timeYearly(end)], [ytick(end) ytick(1)], 'Color', 'k', 'LineWidth', 1.5)

        ax = gca;
        ax.FontName = 'Arial'; % Microsoft YaHei 'Time New Roman'

        ax.XMinorTick = 'off'; ax.YMinorTick = 'off'; % 开启次刻度线
        ax.XAxis.MinorTickValues = (1980:1:2105);
        ax.TickLength = [0.015 0.01]; %刻度线长度      set(gca,'ticklength', [0.02 0.01]);
        ax.XColor = 'k'; ax.YColor = 'k'; % 设置刻度线颜色
        ax.FontSize = 12;
        ax.LineWidth = 1.5;

        % legend
        lgd = legend(varNames, 'Fontsize', 8, 'Location', 'north', 'NumColumns', 3);
        legend('boxoff') %删除图例背景和轮廓
        lgd_inf = get(lgd);

    end

    title_txt = {['Level:', num2str(toaSfc{2})], ['Data: ', level.time1{exmNum_MME}(6:end - 11), ' ', level.time1{exmNum_MME}(end - 9:end - 1)], ['Method: ', detrendSwich{detrendCode}]};
    sgtitle(title_txt, 'FontWeight', 'normal') %,'Units', 'normalized', 'position', [0.5, 0.98]

end

mPath.uniOutput = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/proj3_PaperFig/v0.3/Fig6_7TimeSeries_analysis_landMean_RadEffectAll'); %['dRTs_', lower(mlabels.level)],
mPath.Output1 = fullfile(mPath.uniOutput);
auto_mkdir(mPath.Output1)

% save figures
% figName = [level.time1{exmNum_MME}(6:end - 1), '_', level.model2{mdlNum}, '_', esmName{esmNum}];
% figurePath = [mPath.Output1, '/', figName, '.eps'];
% export_fig(gcf,figurePath,'-r600','-cmyk')

% figName = ['MME1_', level.time1{exmNum_MME}(6:end - 1), '_', detrendSwich{detrendCode}];
% figurePath = [mPath.Output1, '/', figName];
% save_png(figurePath)%high resolution
% close gcf

eval(['cd ', nowpath]);
t = toc; disp(t)
