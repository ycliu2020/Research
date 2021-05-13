%%---------------------------------------------------------
% Author       : LYC
% Date         : 2021-05-02 16:59:49
% LastEditors  : Please set LastEditors
% Description  :
% FilePath     : /code/p1_processObserveData/ERA/plot/radTimeSeries/s3_radTS_ERA_rhsTs.m
%
%%---------------------------------------------------------
run '/home/liuyc/lib/tools/matlab/myTools/autoScript/preLoadVar.m'
chooseData='ERA5';
[readme, level, tLin, vars] = obsParameters(chooseData);

detrendCode = 1;
detrendSwich = {'no detrend', 'Detrend'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% different time series, 1mean 2000-03 to 2018-02(18*12). 2 mean 200207-201706(15*12)
for exmNum = 5:5 % '200003-201802', '200207-201706', '200003-201402', '200001-201412', '198001-201412', '200003-202011'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data path
    exmPath = fullfile(level.obsPath, tLin.time{exmNum}, level.dataName);
    varsPath = fullfile(exmPath, level.standVarPath{1}); %rawdata
    dvarsPath = fullfile(exmPath, level.standVarPath{2}); %anomaly
    dvarsTrendPath = fullfile(exmPath, level.standVarPath{3}); %anomaly_trend
    kernelCalPath = fullfile(exmPath, level.standVarPath{4}); % kernelCal
    dradEffectPath = fullfile(exmPath, level.standVarPath{5}); %radEffect
    figDataPath = fullfile(exmPath, level.standVarPath{8});

    % outPutPath = fullfile(level.figOfSketch, level.dataName, 'radTimeSeries'); % '/home/liuyc/Research/P02.Ts_change_research/figure/figOfSketch/';
    % auto_mkdir(outPutPath)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% load and read
    load([dradEffectPath, 'global_vars.mat']) % 'lon_f', 'lat_f', 'timeEssmble', 'time', 'plev_k', 'readme', 'timeseries', 'MME_Models', 'dhfssMask', 'dhflsMask'
    ntime = length(time);
    % load save one var
    load([figDataPath, 'regionalVarsRad_sfc_cld.mat']) % 'dR_cloudMask', 'dR_residualMask', 'dR_albMask', 'dR_husMask', 'dR_taMask', 'dR_tsMask''dhfssMask', 'dhflsMask', 'drhsMask', 'drhsKernMask','dlwUPRawDataMask'
    regionCode = {'world', 'CHNeast', 'EURwest', 'USeast'}; %'CHNeast', 'EURwest', 'USeast'
    varibalCode = {'dR_tsMask', 'drhsMask', 'dR_taMask',...
    'dR_albMask', 'dR_residualMask', 'dR_cloudMask',...
    'dhfssMask', 'dhflsMask'};

    dR_tsMask.world = -dR_tsMask.world;
    dR_tsMask.CHNeast = -dR_tsMask.CHNeast;
    dR_tsMask.USeast = -dR_tsMask.USeast;
    dR_tsMask.EURwest = -dR_tsMask.EURwest;

    varUsed = zeros(nlonf, nlatf, ntime, length(varibalCode), length(regionCode));
    sizeVarUsed = size(varUsed);
    sizeLon = sizeVarUsed(1); sizeLat = sizeVarUsed(2); sizeTime = sizeVarUsed(3);
    sizeVar = sizeVarUsed(4); sizeRegion = sizeVarUsed(5);

    for regionNum = 1:sizeRegion

        for varNum = 1:sizeVar
            varUsed(:, :, :, varNum, regionNum) = eval([varibalCode{varNum}, '.', regionCode{regionNum}]);
        end

    end
    varUsed(:, :, :, 3, :) =  varUsed(:, :, :, 2, :) + varUsed(:, :, :, 7, :) + varUsed(:, :, :, 8, :);% drhs+LH+SH
    disp('plz check if res minus correct value or not!')
    varUsed(:, :, :, 5, :) = varUsed(:, :, :, 1, :) - varUsed(:, :, :, 2, :);% dR_ts-drhs
    sizeVar = 6;
    % global annual mean
    varUsedYearly_weightMean = glbAnlMean(varUsed, sizeLon, sizeLat, sizeTime, sizeVar, sizeRegion);

    % detrend
    varUsedYearly_weightMean2 = varUsedYearly_weightMean;

    for regionNum = 1:sizeRegion

        for varNum = 1:sizeVar
            varUsedYearly_weightMean2(:, varNum, regionNum, 2) = detrend(varUsedYearly_weightMean(:, varNum, regionNum));
        end

    end

    % x axis: time
    startMonth = 1;
    time = datevec(time);
    timeYearly = unique(time(:, 1));

    if startMonth ~= 1
        timeYearly = timeYearly(1:end - 1);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plot Part
    plotPart(exmNum, detrendCode, timeYearly,  varUsedYearly_weightMean2, sizeLon, sizeLat, sizeTime, sizeVar, sizeRegion);
    title_txt = {['Level:', num2str(toaSfc{2})], ['Data: ', level.dataName, ' ', tLin.time{exmNum}], ['Method: ', detrendSwich{detrendCode}]};
    sgtitle(title_txt, 'FontWeight', 'normal') %,'Units', 'normalized', 'position', [0.5, 0.98]

end

function varUsedYearly_weightMean = glbAnlMean(varUsed, sizeLon, sizeLat, sizeTime, sizeVar, sizeRegion)
    run '/home/liuyc/lib/tools/matlab/myTools/autoScript/preLoadVar.m'

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

end

function [] = plotPart(exmNum, detrendCode, timeYearly, varUsedYearly_weightMean2, sizeLon, sizeLat, sizeTime, sizeVar, sizeRegion)
    run '/home/liuyc/lib/tools/matlab/myTools/autoScript/preLoadVar.m'

    varNames = {['-d', '\itR', '\rm_{Ts}'], ['dRHeating'], ['dRHeating+LH+SH'], ['d', '\itR', '\rm_{Alb}'], ['\rm{Residual}'], ['d', '\itR', '\rm_{cld}']};
    unitName = {'W m^{-2}'};
    regionName = {['Global land'], ['Eastern CHN'], ['Western EUR'], ['Eastern US']};


    set(0, 'defaultfigurecolor', 'w'); %设置画布底色为白色
    set(0, 'DefaultFigureVisible', 'on')
    ss = get(0, 'ScreenSize');
    h = figure('Position', [10 18 1072 719]); %[离左边缘 离下边缘 自身宽 自身高][ss(4) / 2 - 100 ss(3) / 35 ss(3) / 5 * 2.5, (ss(4) - 80) / 5 * 4]
    % varColor = {'#F08212',  'none', '#67C1EE',  'none', '#000000', 'none'}; %橘色 深蓝 浅蓝 浅绿 黑色 红色'#67C1EE', #4450A1
    varColor = {'#F08212', '#4450A1', 'none', 'none', '#000000', 'none'}; %橘色 深蓝 浅蓝 浅绿 黑色 红色'#67C1EE', #4450A1

    % time series
    timeSer = [1985 1990 1995 2000 2005 2010 2015 2025 2035 2045 2055 2065 2075 2085 2095 2105];

    if exmNum == 1
        timeSer = [2000 2002 2006 2010 2014];
    elseif exmNum == 5
        timeSer = [1982 1990 1998 2006 2014 2020 2025 2035 2045 2055 2065 2075 2085 2095 2105];
    end

    char_timeSer = cellstr(string(timeSer));

    zeroLine = zeros(length(timeYearly), 1);
    locDemo = randn(2);
    % detrendCode=1;
    varUsedYearly_weightMean = varUsedYearly_weightMean2(:, :, :, detrendCode);
    for regionNum = 1:4 % 'world', 'CHNeast', 'EURwest', 'USeast'
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

        if exmNum == 1
            char_timeSer{1} = [];
        end

        xticklabels([])

        if regionNum == 3 || regionNum == 5
            xticklabels(char_timeSer)
        end

        xlim([timeYearly(1) timeYearly(end)])

        % y axis
        ymax = 10;


        if exmNum == 5
            if regionNum==1
                ymax = 4;

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

        if regionNum == 3 || regionNum == 5
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


end
