%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-08-31 17:00:15
% LastEditTime : 2021-05-04 21:48:46
% LastEditors  : Please set LastEditors
% Description  : 同时画时间序列和相关性分布图
% FilePath     : /code/p3_paperFigIntegrate/Fig6_7_timeSeriANS_tsVsVarsRad/s1_singleModel_timeSeriANS_tsVsVarsRad.m
%
%%---------------------------------------------------------
clc; clear; tic;
nowpath = pwd;
% load mask map
run '/home/liuyc/lib/tools/matlab/myTools/autoScript/preLoadVar.m'

MME1Path = '/data1/liuyincheng/CMIP6-process/z_ensembleMean/MME1/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment
for exmNum = [1] %1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    %CAMS-CSM1-0 didn't have sfc clear sky radiation, delete it
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % exmPath
    exmPath = ['/data1/liuyincheng/CMIP6-process/', level.time1{exmNum}]; %/data1/liuyincheng/CMIP6-process/2000-2014/
    MMEPath = [MME1Path, level.time1{exmNum}];
    load([MMEPath, '/radEffect/global_vars.mat']) % lat_f, lon_f, MME_Models, plev_k, readme, time, timeEssmble, timeseries

    mPath.uniOutput = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/proj3_PaperFig/v0.3/Fig6_7TimeSeries_analysis_landMean_RadEffectAll', level.time1{exmNum}); %['dRTs_', lower(mlabels.level)],
    mPath.Output1 = fullfile(mPath.uniOutput);
    auto_mkdir(mPath.Output1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model
    for mdlNum = 3:3%length(level.model2)

        % model path
        mdlName = level.model2{mdlNum};

        if ismember(mdlName, MME_Models.name) ~= 1
            continue
        end

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
            % data path
            varsPath = fullfile(esmPath, level.process3{1}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/rawdata_regrid
            dvarsPath = fullfile(esmPath, level.process3{2}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/anomaly
            dvarsTrendPath = fullfile(esmPath, level.process3{3}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/anomaly_trend
            kernelPath = fullfile(esmPath, level.process3{5}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/kernelsCal
            dradEffectPath = fullfile(esmPath, level.process3{6}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/radEffect/
            dnonLocalCldPath = fullfile(esmPath, level.process3{8}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/non_localCld/
            dvsTsEffectPath = fullfile(esmPath, level.process3{9}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/vsTsEffect/
            load([dradEffectPath, 'global_vars.mat']) % lat_f lon_f time plev_k readme
            ntime = length(time.date);

            % cal surface temperature
            load([varsPath, 'ts.mat']) % ts
            load([dvarsPath, 'dts.mat']) %dts
            dts = autoRegrid3(lon_k, lat_k, time.date, dts, lon_f, lat_f, time.date);
            ts = autoRegrid3(lon_k, lat_k, time.date, ts, lon_f, lat_f, time.date);

            % cal rhs(RHeating)
            load([dvarsPath, 'drlds.mat']) % surface_downwelling_longwave_flux_in_air
            load([dvarsPath, 'drsds.mat']) % surface_downwelling_shortwave_flux_in_air
            load([dvarsPath, 'drsus.mat']) % surface_upwelling_shortwave_flux_in_air
            dR_swnet = drsds - drsus; % sfc net shortwave flux
            drhs = drlds + dR_swnet; % equilibrium equation's RHS, nearly equal to sfc upward rad
            drhs = autoRegrid3(lon_k, lat_k, time.date, drhs, lon_f, lat_f, time.date);

            % cloud radRffect
            load([dradEffectPath, 'dR_cloud.mat'])

            % TotalEffect
            load([dradEffectPath, 'dradEffect_sfc_cld.mat']) % 'totalEffect', 'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', 'taEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect', 'mainEffect'
            load([dradEffectPath, 'model_dradEffect.mat']) % 'dR_allsky', 'l_rad', 's_rad', 'dR_clr', 'readme_realradEfect'
            load([dradEffectPath, 'dR_residual_cld_sfc.mat']) % dR_residual_cld_sfc
            load([dradEffectPath, 'dR_residual_cld_toa.mat']) % dR_residual_cld_toa
            dR_netSfc = squeeze(dR_net_cld(:, :, :, 1));
            dR_netTOA = squeeze(dR_net_cld(:, :, :, 2));

            % vars tsEffect
            load([dvsTsEffectPath, 'dTs_x_sfc.mat']) % dTs_alb, dTs_cloud, dTs_hus, dTs_ta, dTs_residual

            %
            load([dvarsPath, 'dhfls.mat']) % Surface Upward Latent Heat Flux
            load([dvarsPath, 'dhfss.mat']) % Surface Upward Sensible Heat Flux
            % unite define a vector component which is positive when directed downward
            dhfls = autoRegrid3(lon_k, lat_k, time.date, dhfls, lon_f, lat_f, time.date);
            dhfss = autoRegrid3(lon_k, lat_k, time.date, dhfss, lon_f, lat_f, time.date);

            dhFlux = -dhfls - dhfss; % LH+SH

            varUsed = zeros(nlonf, nlatf, ntime, 5);

            % fig1
            dRheating = dR_cloud_sfc + husEffect + taEffect + albEffect + dR_residual_cld_sfc;

            regionCode = {'world', 'CHNeast', 'EURwest', 'USeast'}; %'CHNeast', 'EURwest', 'USeast'
            regionOrder = {'world', 'china east', 'EUR west', 'USA east'};
            varibalCode = {'dlwUPMask', 'dR_husMask', 'dR_taMask', 'dR_albMask', 'dR_residualMask', 'dR_cloudMask'};

            varUsed(:, :, :, 1) = -tsEffect;
            varUsed(:, :, :, 2) = husEffect;
            varUsed(:, :, :, 3) = taEffect;
            varUsed(:, :, :, 4) = albEffect;
            varUsed(:, :, :, 5) = dR_residual_cld_sfc;
            varUsed(:, :, :, 6) = dR_cloud_sfc;
            % varUsed(:, :, :, 7) = dhfls.mat+dhfss;

            % varUsed(:, :, :, 7) = dR_netSfc; %sum(varUsed(:,:,:,1:6),4);
            % varUsed(:, :, :, 8) = tsEffect+dRheating;
            % varUsed(:, :, :, 9) = dRheating;

            % transfor to yearly Data
            sizeVarUsed = size(varUsed);
            sizeLon = sizeVarUsed(1); sizeLat = sizeVarUsed(2); sizeTime = sizeVarUsed(3); sizeVar = sizeVarUsed(4); sizeRegion = length(regionCode);

            varUsedYearly = zeros(sizeLon, sizeLat, sizeTime / 12, sizeVar);

            for varNum = 1:sizeVar
                varUsedYearly(:, :, :, varNum) = monthlyToYearly(squeeze(varUsed(:, :, :, varNum)));
            end

            % compute each region
            varUsedYearly_temp = varUsedYearly;
            sizeVarUsedYearly = size(varUsedYearly);
            varUsedYearly_weightMean = zeros(sizeVarUsedYearly(3), sizeVar, sizeRegion);

            for regionNum = 1:sizeRegion
                % mask part.
                areaStr = regionOrder{regionNum};

                for varNum = 1:sizeVar
                    [varUsedYearly_temp(:, :, :, varNum), ~, ~] = maskArea(squeeze(varUsedYearly(:, :, :, varNum)), lat_f, latRange, -latRange, areaStr);
                end

                % cal areaMeanLatWeight

                for varNum = 1:sizeVar

                    for timeNum = 1:sizeVarUsedYearly(3)
                        varUsedYearly_weightMean(timeNum, varNum, regionNum) = areaMeanLatWeight(squeeze(squeeze(varUsedYearly_temp(:, :, timeNum, varNum))), lat_f);
                    end

                end

            end
            % varUsedYearly_weightMean(:, 5, :)= varUsedYearly_weightMean(:, 1, :)- (sum(varUsedYearly_weightMean,2)-varUsedYearly_weightMean(:, 5, :)-varUsedYearly_weightMean(:, 1, :));
            % detrend
            for regionNum = 1:sizeRegion

                for varNum = 1:sizeVar
                    varUsedYearly_weightMean(:,varNum, regionNum) = detrend(varUsedYearly_weightMean(:, varNum, regionNum));

                end
            end
            time.vec = datevec(time.date);
            timeYearly = unique(time.vec(:, 1));

            % % cal cc of areaMeanLatWeight
            % sizeWeightMean = size(varUsedYearly_weightMean);
            % cc_weightMean = cell(1, sizeWeightMean(2));
            % pp_weightMean = cell(1, sizeWeightMean(2));

            % for varNum = 1:sizeWeightMean(2)
            %     [cc0, pp0] = corrcoef(varUsedYearly_weightMean(:, 1), varUsedYearly_weightMean(:, varNum), 'Rows', 'complete');
            %     cc_weightMean{varNum} = roundn(cc0(1, 2), -2); %保留俩位小数
            %     pp_weightMean{varNum} = pp0(1, 2); % confidence interval
            % end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% plot Part
            varNames = {['d', '\itLW', '\rm_{up}'], ['d', '\itR', '\rm_{WV}'], ['d', '\itR', '\rm_{Ta}'], ['d', '\itR', '\rm_{Alb}'], ['d', '\itR', '\rm_{Residual}'], ['d', '\itR', '\rm_{cld}']};
            unitName = {'W m^{-2}'};
            regionName = {['Global land'], ['Western EUR']; ['Eastern CHN'], ['Eastern US']};
            set(0, 'defaultfigurecolor', 'w'); %设置画布底色为白色
            set(0, 'DefaultFigureVisible', 'on')
            ss = get(0, 'ScreenSize');
            h = figure('Position', [10 18 1072 719]); %[离左边缘 离下边缘 自身宽 自身高][ss(4) / 2 - 100 ss(3) / 35 ss(3) / 5 * 2.5, (ss(4) - 80) / 5 * 4]

            varColor = {'#F08212', '#4450A1', '#67C1EE', '#90C64D', '#000000', '#C2162E'}; %橘色 深蓝 浅蓝 浅绿 黑色 红色

            % time series
            timeSer = [1985 1990 1995 2000 2005 2010 2015 2025 2035 2045 2055 2065 2075 2085 2095 2105];

            if exmNum == 1
                timeSer = [2000 2002 2006 2010 2014];
            elseif exmNum == 2
                timeSer = [1982 1990 1998 2006 2014 2020 2025 2035 2045 2055 2065 2075 2085 2095 2105];
            end

            char_timeSer = cellstr(string(timeSer));

            zeroLine = zeros(length(timeYearly), 1);
            locDemo = randn(2);

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

                if regionNum == 3 || regionNum == 4
                    xticklabels(char_timeSer)
                end

                xlim([timeYearly(1) timeYearly(end)])

                % y axis
                ymax = 10;

                if exmNum == 1
                    ymax = 8;
                elseif exmNum == 2
                    ymax = 5;
                end

                if regionNum == 1 && exmNum == 1
                    ymax = 2;
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

            title_txt = {['Level:', num2str(toaSfc{2}), ', Era: ', level.time1{exmNum}(1:end - 10), ' ', level.time1{exmNum}(end - 9:end - 1)], ['Model:', level.model2{mdlNum} ', Ensemble: ', esmName{esmNum}], [num2str(latRange), 'N-', num2str(latRange), 'S land, global mean '], ''};
            sgtitle(title_txt, 'FontWeight', 'normal') %,'Units', 'normalized', 'position', [0.5, 0.98]


            % text(lgd_inf.Position(1) - 0.12, lgd_inf.Position(2) + 0.085, ['cc= ', num2str(cc_weightMean{2})], 'FontWeight', 'normal', 'Interpreter', 'none', 'Units', 'normalized')
            % hold on

            % save figures
            % figName = [level.time1{exmNum}(6:end - 1), '_', level.model2{mdlNum}, '_', esmName{esmNum}];
            % figurePath = [mPath.Output1, '/', figName, '.eps'];
            % export_fig(gcf,figurePath,'-r600','-cmyk')

            % figName = [level.time1{exmNum}(6:end - 1), '_', level.model2{mdlNum}, '_', esmName{esmNum}];
            % figurePath = [mPath.Output1, '/', figName];
            % save_png(figurePath)%high resolution
            % close gcf

        end

    end

end

eval(['cd ', nowpath]);
t = toc; disp(t)
