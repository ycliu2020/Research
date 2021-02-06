%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-08-31 17:00:15
% LastEditTime : 2021-01-14 16:37:25
% LastEditors  : Please set LastEditors
% Description  : 同时画时间序列和相关性分布图
% FilePath     : /code/p2_processCMIP6Data/s2.radEffTrend/timeSeriAnalysis/timeSeriANS_radClose.m
%
%%---------------------------------------------------------
clc; clear; tic;
nowpath = pwd;
dbstop if error

% load mask map
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')

latRange = 90; % Latitude range
lon1 = [2.5 357.5]; lat1 = [-latRange + 1 latRange - 1]; % world area
toaSfc = {'toa', 'sfc'};
lon_k = 0:2.5:357.5; nlonk = length(lon_k); % kernel lat lon
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment
for exmNum = 2:2%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    %CAMS-CSM1-0 didn't have sfc clear sky radiation, delete it
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % exmPath
    exmPath = ['/data1/liuyincheng/cmip6-process/', level.time1{exmNum}]; %/data1/liuyincheng/cmip6-process/2000-2014/
    mPath.uniOutput = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/proj2_cmip6Result/TimeSeries_analysis/radClosure', level.time1{exmNum}); %['dRTs_', lower(mlabels.level)],
    mPath.Output = fullfile(mPath.uniOutput);
    auto_mkdir(mPath.Output)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model
    for mdlNum = 3:3%length(level.model2)
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
        for esmNum = specificNum:specificNum%1:1%length(esmName)% note that r1i1p1 sometime not the first folder
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
            dvsTsEffectPath = fullfile(esmPath, level.process3{9}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/vsTsEffect/
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
            load([dradEffectPath, 'dradEfect_sfc_clr.mat'])% 'totalEffect', 'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', 'taEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect', 'mainEffect'
            load([dradEffectPath, 'real_dradEfect.mat'])% 'dR_allsky', 'l_rad', 's_rad', 'dR_clr', 'readme_realradEfect'
            load([dradEffectPath, 'dR_residual_cld_sfc.mat'])% dR_residual_cld_sfc
            load([dradEffectPath, 'dR_residual_clr_sfc.mat'])% dR_residual_cld_sfc
            load([dradEffectPath, 'dR_residual_cld_toa.mat'])% dR_residual_cld_toa
            dR_net_cld_sfc = squeeze(dR_allsky(:, :, :, 1));
            dR_net_cld_toa = squeeze(dR_allsky(:, :, :, 2));

            dR_net_clr_sfc = squeeze(dR_clr(:, :, :, 1));
            dR_net_clr_toa = squeeze(dR_clr(:, :, :, 2));
            % vars tsEffect
            % load([dvsTsEffectPath, 'dTs_x_sfc.mat'])% dTs_alb, dTs_cloud, dTs_hus, dTs_ta, dTs_residual

            %
            load([dvarsPath, 'dhfls.mat'])% Surface Upward Latent Heat Flux
            load([dvarsPath, 'dhfss.mat'])% Surface Upward Sensible Heat Flux
            % unite define a vector component which is positive when directed downward
            dhfls = autoRegrid3(lon_k, lat_k, time.date, dhfls, lon_f, lat_f, time.date);
            dhfss = autoRegrid3(lon_k, lat_k, time.date, dhfss, lon_f, lat_f, time.date);

            dhFlux = -dhfls - dhfss; % LH+SH

            varUsed = zeros(nlonf, nlatf, ntime, 2);

            % fig1
            dRheating=dR_cloud_sfc+husEffect+taEffect+albEffect+dR_residual_cld_sfc;
            
            varUsed(:, :, :, 2) = dR_net_clr_sfc;
            varUsed(:, :, :, 1) = tsEffect+husEffect+taEffect+albEffect;
            varUsed(:, :, :, 3) = varUsed(:, :, :, 2)-varUsed(:, :, :, 1);
    
            sizeVarUsed = size(varUsed);
            sizeLon = sizeVarUsed(1); sizeLat = sizeVarUsed(2); sizeTime = sizeVarUsed(3); sizeVar = sizeVarUsed(4);

            varNames = { ['d','\itR','\rm_{net}(kernel calc)'], ['d','\itR','\rm_{net}'], 'Residual'};%, 'dR_{albedo}', 'dR_{residual}', 'dR_{cloud}'
            yLabel = {'W m^{-2}', 'W m^{-2}', 'W m^{-2}', 'W m^{-2}', 'W m^{-2}', 'W m^{-2}', 'W m^{-2}', 'W m^{-2}', 'W m^{-2}'};
            varColor = { '#F08212', '#C2162E', '#000000', '#4450A1',  '#90C64D', '#000000', '#C2162E'}; %橘色 深红 黑色 深蓝 浅蓝 浅绿 红色

 
            varUsedYearly=varUsed;
            areaStr = 'world';
            for varNum = 1:sizeVar
                [varUsedYearly(:,:,:,varNum), ~, ~] = maskArea(squeeze(varUsedYearly(:,:,:,varNum)), lat_f, latRange, -latRange, areaStr);
            end

            % % method 1 cal areaMeanLatWeight
            sizeVarUsedYearly = size(varUsedYearly);
            varUsedYearly_weightMean = zeros(sizeVarUsedYearly(3), sizeVar);

            for varNum = 1:sizeVar

                for timeNum = 1:sizeVarUsedYearly(3)
                    varUsedYearly_weightMean(timeNum, varNum) = areaMeanLatWeight(squeeze(squeeze(varUsedYearly(:, :, timeNum, varNum))), lat_f);
                end

            end

            time.vec = datevec(time.date);
            timeYearly = unique(time.vec(:, 1));
            timeYearly=time.date;
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
            timeSer = [1980 1985 1990 1995 2000 2005 2010 2015 2025 2035 2045 2055 2065 2075 2085 2095 2105];
            char_timeSer = cellstr(string(timeSer));
            timeNum = zeros(length(timeSer),1);
            for timeCount = 1:length(timeSer)
                timeNum(timeCount,1) =datenum(timeSer(timeCount),01,01);
            end
            timeNum(1)=timeYearly(1);
            timeNum(8)=timeYearly(end);
            % fig set
            set(0, 'DefaultFigureVisible', 'on')
            ss = get(0, 'ScreenSize');
            coef_amplify = 1.5;
            % h = figure('Position', [ss(4) / 2 - 100 ss(3) / 35 ss(3) / 5 * f_col * coef_amplify (ss(4) - 80) / 5 * f_row * coef_amplify]);
            h = figure('Position', [-1075 435 889 569]);
            % clf reset;
            set(h, 'Color', [1 1 1]);

            pPlot=cell(1,3);
            for varNum = 1:sizeWeightMean(2)
                varUsedTemp = squeeze(varUsedYearly_weightMean(:, varNum));
                pPlot{varNum}=plot(timeYearly, varUsedTemp', 'color', varColor{varNum}, 'LineWidth', 2);
                hold on
            end

            % zero line
            zeroLine = zeros(length(timeNum), 1);
            plot(timeNum, zeroLine', '--k','LineWidth', 2)
            hold on;

            % set x axes
            datetick('x', 'yyyy', 'keepticks'); % 用日期的格式显示横坐标
            xlim([timeYearly(1) timeYearly(end)])
            xlabel('time line')
            xticks(timeNum)
            xticklabels(char_timeSer)

            ylabel(yLabel{2})
            % set y axes
            ymax = 15;

            if exmNum == 1 
                ymax = 6;
            elseif  exmNum == 2
                ymax = 5;
                yticks([-5 -4 -3 -2 -1 0 1 2 3 4 5])
                yticklabels({'-5','-4','-3','-2','-1','0','1','2','3','4','5'})
            end

            ylim([-ymax ymax])

            ax = gca;
            ax.FontName='Microsoft YaHei';% Microsoft YaHei 'Time New Roman'
            ax.XMinorTick = 'off'; ax.YMinorTick = 'off'; % 开启次刻度线
            ax.TickLength = [0.015 0.01]; %刻度线长度      set(gca,'ticklength', [0.02 0.01]);
            ax.XColor = 'k'; ax.YColor = 'k'; % 设置刻度线颜色
            ax.FontSize = 18;
            ax.LineWidth = 1.5;
            title_txt = {['Level:', num2str(toaSfc{2}), ', Era: ', level.time1{exmNum}(1:end - 10),' ',level.time1{exmNum}(end - 9:end-1)], ['Model:', level.model2{mdlNum} ', Ensemble: ', esmName{esmNum}], [num2str(latRange),'N-',num2str(latRange),'S land, global mean '], ''};
            title(title_txt, 'FontWeight', 'normal')
            
            % 去掉上边框和右边框的刻度
            box off
            xtick = get(gca,'XTick');
            ytick = get(gca,'YTick');
            line([xtick(1),xtick(end)],[ytick(end) ytick(end)],'Color','k','LineWidth', 1.5)
            line([timeYearly(end),timeYearly(end)],[ytick(end) ytick(1)],'Color','k','LineWidth', 1.5)

            
            % lgdTxtTmp(1:6)={',cc='};
            % cc_weightMeanTxt=cellfun(@num2str,cc_weightMean,'UniformOutput',false);
            % lgdTxt=cellfun(@strcat,varNames,lgdTxtTmp,cc_weightMeanTxt,'UniformOutput',false);
            % lgd = legend(lgdTxt, 'Fontsize', 8, 'Location', 'northwest', 'NumColumns', 3);
            lgd = legend([pPlot{1} pPlot{2} pPlot{3}], varNames, 'Fontsize', 16, 'Location', 'northwest', 'NumColumns', 1);
            legend('boxoff')%删除图例背景和轮廓
            lgd_inf = get(lgd);
            % text(lgd_inf.Position(1) - 0.12, lgd_inf.Position(2) + 0.085, ['cc= ', num2str(cc_weightMean{2})], 'FontWeight', 'normal', 'Interpreter', 'none', 'Units', 'normalized')
            % hold on
            

            % save figures
            figName = [level.time1{exmNum}(1:end - 1), '_', level.model2{mdlNum}, '_', esmName{esmNum}];
            figurePath = [mPath.Output, '/', figName,'.eps'];
 
            export_fig(gcf,figurePath,'-r600','-cmyk')

            % save_png(figurePath)%high resolution
            % print(gcf,figurePath,'-depsc2','-r600')
            % close gcf

        end

    end

end

eval(['cd ', nowpath]);
t = toc; disp(t)