%%---------------------------------------------------------
% Author       : LYC
% Date         : 2021-04-22 09:41:39
% LastEditors  : Please set LastEditors
% Description  :
% FilePath     : /code/p3_paperFigIntegrate/Fig3_timeSeries_tsVsRheatingRad/Fig_test/timeSeries_LHSH_chinaEast_annual.m
%
%%---------------------------------------------------------
clc; clear; tic;
nowpath = pwd;
dbstop if error
% load global vars
run '/home/liuyc/lib/tools/matlab/myTools/autoScript/preLoadVar.m'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment
for exmNum = 2:2
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    %1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    %CAMS-CSM1-0 didn't have sfc clear sky radiation, delete it
    % exmPath
    exmPath = ['/data1/liuyincheng/CMIP6-process/', level.time1{exmNum}]; %/data1/liuyincheng/cmip6-process/2000-2014/
    mPath.uniOutput = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/proj3_PaperFig/v0.3/Fig0_CMIP6_200003-201402_world_testRes/'); %['dRTs_', lower(mlabels.level)],
    mPath.Output = fullfile(mPath.uniOutput);
    auto_mkdir(mPath.Output)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model
    tempCount = 1;

    for mdlNum = 1:length(level.model2) %1:length(level.model2)
        temp = MME_Models.name;

        if sum(strcmp(level.model2{mdlNum}, temp)) == 0
            disp(['jump model']);
            continue
        end

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
        for esmNum = specificNum:specificNum %1:1 %length(esmName) % note that r1i1p1 sometime not the first folder
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% load and read
            esmPath = fullfile(mdlPath, esmName{esmNum, 1});

            dvarsPath = fullfile(esmPath, level.process3{2});
            dradEffectPath = fullfile(esmPath, level.process3{6});
            figDataPath = fullfile(esmPath, 'FigData/');

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% load and read
            load([dradEffectPath, 'global_vars.mat']) % 'lon_f', 'lat_f', 'timeEssmble', 'time', 'plev_k', 'readme', 'timeseries', 'MME_Models'

            % load MME and save one var
            load([figDataPath, 'regionalTsRHeating.mat']) % dlwUPMask, drhsKernMask, drhsMask  .world .CHNeast  .EURwest .USeast
            ntime = length(drhsMask.world);

            % index chose
            regionCode = {'CHN', 'EURwest', 'USeast'}; %'CHNeast', 'EURwest', 'USeast'
            varibalCode = {'dhfssMask', 'dhflsMask', 'drhsMask'};
            varUsed = zeros(nlonf, nlatf, ntime, length(varibalCode), length(regionCode)); % dim:[lon, lat, time, region, varibal, dataset]
            sizeVarUsed = size(varUsed);
            sizeLon = sizeVarUsed(1); sizeLat = sizeVarUsed(2); sizeTime = sizeVarUsed(3);
            sizeVar = sizeVarUsed(4); sizeRegion = sizeVarUsed(5);

            for regionNum = 1:sizeRegion

                for varNum = 1:sizeVar
                    varUsed(:, :, :, varNum, regionNum) = eval([varibalCode{varNum}, '.', regionCode{regionNum}]);
                end

            end

            temp = varUsed(:, :, :, 3, :);
            varUsed(:, :, :, 3, :) = varUsed(:, :, :, 1, :) + varUsed(:, :, :, 2, :); % LH+SH
            varUsed = -varUsed;
            varUsed(:, :, :, 4, :) = temp;
            sizeVar = 4;
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

            % x axis: time
            startMonth = 1;
            time.vec = datevec(time.date);
            timeYearly = unique(time.vec(:, 1));

            if startMonth ~= 1
                timeYearly = timeYearly(1:end - 1);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% plot Part
            % radName = {['d', '\itR', '\rm_{heating,Kern} cc='], ['d', '\itR', '\rm_{heating} cc='], ['d', '\itLW', '\rm_{up} ']};
            unitName = {'W m^{-2}'};
            regionName = {['CHN'], ['Western CHN'], ['Eastern US']};
            varibalName = {'SH', 'LH', 'SH+LH', 'Rheating'};

            set(0, 'defaultfigurecolor', 'w'); %设置画布底色为白色
            set(0, 'DefaultFigureVisible', 'on')

            ss = get(0, 'ScreenSize');
            h = figure('Position', [-1075 188 934 764]); %[离左边缘 离下边缘 自身宽 自身高][ss(4) / 2 - 100 ss(3) / 35 ss(3) / 5 * 2.5, (ss(4) - 80) / 5 * 4]

            zeroLine = zeros(length(timeYearly), 1);
            subNum = 1;

            for varNum = 1:4 % LH SH LH+SH

                for regionNum = 1:3 % 'CHNeast', 'world'

                    subplot_yc(4, 3, varNum, regionNum);
                    hold on
                    lineWdth = 1.5;
                    plot(timeYearly, squeeze(varUsedYearly_weightMean(:, varNum, regionNum)), 'color', '#BC144A', 'LineWidth', lineWdth) % 真实值值
                    plot(timeYearly, zeroLine', 'k--', 'LineWidth', 1.5)

                    % x axis
                    if exmNum == 1
                        timeSer = [1985 1990 1995 2002 2006 2010 2014 2025 2035 2045 2055 2065 2075 2085 2095 2105];
                    elseif exmNum == 2
                        timeSer = [1980 1990 2000 2010 2015 2020 2025 2035 2045 2055 2065 2075 2085 2095 2105];
                    end

                    char_timeSer = cellstr(string(timeSer));

                    xlim([timeYearly(1) timeYearly(end)]);

                    xticks(timeSer)
                    xticklabels([])

                    if varNum == length(varibalName)
                        xticklabels(char_timeSer)
                    end

                    % y axis
                    ymax = 20;
                    ylim([-ymax ymax])

                    % if regionNum == length(regionName)
                    %     ylim([-5 5])
                    % end

                    % y label
                    if regionNum == 1
                        yLabel = [varibalName{varNum}; unitName];
                        ylabel(yLabel, 'FontName', 'Arial', 'Fontsize', 17)
                    end

                    % 上标题
                    if varNum == 1
                        text(0.45, 1.15, regionName{regionNum}, 'FontName', 'Arial', 'Fontsize', 17, 'units', 'normalized');
                    end

                    % 添加序号
                    if varNum == length(varibalName)
                        text(0.475, -0.2, ['(', char(96 + subNum), ')'], 'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize', 15, 'units', 'normalized');
                    else
                        text(0.475, -0.1, ['(', char(96 + subNum), ')'], 'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize', 15, 'units', 'normalized');
                    end

                    subNum = subNum + 1;

                    % % legend
                    % gcaGet = get(gca);
                    % axPosit = gcaGet.Position;
                    % axPosit = axPosit + [0.0 0.17 -0.04 -0.18]; %(左边缘, 下边缘, 长度, 宽度)
                    % legend({[radName{1}, num2str(cc_weightMean{1, regionNum, dataNum})], [radName{2}, num2str(cc_weightMean{2, regionNum, dataNum})], radName{3}}, 'FontName', 'Arial', 'Fontsize', 9, 'Location', 'none', 'position', axPosit, 'NumColumns', 2); %'Location', 'best''FontWeight', 'bold',
                    % legend('boxoff') %删除图例背景和轮廓

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

            sgtitle({['Model:', level.model2{mdlNum} ', Ensemble: ', esmName{esmNum}, ', Era: ', level.time1{exmNum}(1:end - 11)]}, 'FontName', 'Arial', 'Fontsize', 18)

        end

    end

end
