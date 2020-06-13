%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-09 15:52:00
% LastEditTime : 2020-06-12 17:05:17
% LastEditors  : LYC
% Description  :
% FilePath     : /Research/p2_processCMIP6Data/s4.radFeedbakIntensity/s1_calRadFeedIntens.m
%
%%---------------------------------------------------------
% 计算全球各种量的反馈强度
% 云, 温度 水汽 反照率
% 注意此脚本保存应该都要加上mask类型的后缀('nomask'or'maskLand')
clear; clc; tic;
% pause(9999)
% model settings
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')% ????????????????word_mapx(:),word_mapy(:)
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')

p_3 = 88.75; % Latitude range
lon1 = [2.5 357.5]; lat1 = [-p_3 + 1 p_3 - 1]; % world area
toaSfc = {'toa', 'sfc'};
maskLandSW = 'nomask'; %{'nomask', 'maskLand'};
areaNum = 1; % world land

for p_1 = 1:4%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    nowpath = pwd;
    %CAMS-CSM1-0 didn't have sfc clear sky radiation, delete it
    [readme, Experiment, level, tLin, mPlev, vars] = modelParameters(p_1);

    % inputPath
    inputPath = ['/data1/liuyincheng/cmip6-process/', level.time1{p_1}]; %/data1/liuyincheng/cmip6-process/amip_2000-2014/

    % model loop
    dvarsPath = [inputPath, level.model2{1}, '/', level.process3{2}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
    load([dvarsPath, 'global_vars.mat'])% latf lonf time plevf readme
    lat = 88.75:-2.5:-88.75; nlat = length(lat);
    lon = lonf; nlon = length(lon);
    nlonf = length(lonf); nlatf = length(latf);
    periodDtsAssemble = zeros(length(level.model2), 2);
    lamda_globalAssemble = zeros(5, length(level.model2), 2); %(vars, model, toaSfc)
    lamda_gridAssemble = zeros(nlon, nlat, 5, length(level.model2), 2); %(lon,lat,vars,model,toaSfc)
    dradFeedbackAssemblePath = [inputPath, 'z_assembleData/', level.process3{7}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
    figTestPath = [inputPath, 'z_assembleData/figTest/']; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
    rmdir([inputPath, 'z_assembleData'],'s')
    auto_mkdir(dradFeedbackAssemblePath)
    auto_mkdir(figTestPath)

    for level1 = 1:length(level.model2)
        % load dvars
        dvarsPath = [inputPath, level.model2{level1}, '/', level.process3{2}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
        varsPath = [inputPath, level.model2{level1}, '/', level.process3{1}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/rawdata
        dradEffectPath = [inputPath, level.model2{level1}, '/', level.process3{6}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
        dradFeedbackPath = [inputPath, level.model2{level1}, '/', level.process3{8}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
        auto_mkdir(dradFeedbackPath)
        % vars.lamda={'cloud','wv','ta','alb'};
        %% Part1: load data
        load([dvarsPath, 'global_vars.mat'])% latf lonf time plevf readme
        load([dvarsPath, 'dts.mat'])%
        % regrid dts to 144x72(unite grids)
        [Xlon, Ylat, Ttime] = meshgrid(lat, lon, time);
        [Xlonf, Ylatf, Ttimef] = meshgrid(latf, lonf, time);
        dts = interp3(Xlonf, Ylatf, Ttimef, dts, Xlon, Ylat, Ttime);

        % sky level
        for skyLevel = 1:2% 1 mean toa, 2 mean sfc
            load([dradEffectPath, ['dradEfect_', toaSfc{skyLevel}, '_cld.mat']])%albEffect, husEffect, taEffect, mainEffect, totalEffect, tsEffect, wvlwEffect, wvswEffect
            load([dradEffectPath, ['dR_cloud_', toaSfc{skyLevel}, '.mat']])%dR_cloud_toa

            % add to one var
            tEffect = taEffect + tsEffect;
            varUsed(:, :, :, 1) = dts; %144x72x1800x2

            if skyLevel == 1
                varUsed(:, :, :, 2) = dR_cloud_toa;
            elseif skyLevel == 2
                varUsed(:, :, :, 2) = dR_cloud_sfc;
            end

            varUsed(:, :, :, 3) = husEffect;
            varUsed(:, :, :, 4) = taEffect;
            varUsed(:, :, :, 5) = albEffect;
            varUsedSize = size(varUsed);
            varUsedSize_time = varUsedSize(3);
            varUsedSize_vars = varUsedSize(4);
            % land mask (optional)
            if strcmp(maskLandSW, 'maskLand') == 1

                for var_id = 1:varUsedSize_vars
                    tempVar = squeeze(varUsed(:, :, :, var_id));
                    [tempVar] = maskLand(tempVar, lat, p_3, -p_3, areaNum);
                    varUsed(:, :, :, var_id) = tempVar;
                end

            end

            % global Zonal weighted average (time, vars)
            jiaquan = cosd(lat);
            wei = ones(144, 72); %格点纬度加权

            for latiNum = 1:72
                wei(:, latiNum) = wei(:, latiNum) * jiaquan(latiNum); %格点相对大小
            end

            varUsed = varUsed .* wei;
            varWorldMean = zeros(varUsedSize_time, varUsedSize_vars);

            for var_id = 1:varUsedSize_vars

                for timeNum = 1:varUsedSize_time
                    varWorldMean(timeNum, var_id) = nansum(nansum(varUsed(:, :, timeNum, var_id))) / nansum(nansum(wei));
                end

            end

            dts_global = varWorldMean(:, 1);
            save([dradFeedbackPath, ['dts_groble_', maskLandSW, '.mat']], 'dts_global')

            %% Part2: cal period deltaTs (whole time line)

            X = [ones(size(time)) time]; %dts
            [b, bint, r, rint, stats] = regress(dts_global, X); % 1 mean ts
            k1_ts = b(2); % slop of time series K/时间的单位
            periodDts = k1_ts * (time(end) - time(1));
            readme.maskStatement = 'nomask mean land+sea; mask land mean no sea';
            save([dradFeedbackPath, 'global_vars.mat'], 'lonf', 'latf', 'time', 'plevf', 'readme', 'timeseries', 'modelname')
            save([dradFeedbackPath, ['periodDts_', maskLandSW, '.mat']], 'k1_ts', 'periodDts')
            periodDtsAssemble(level1, 1) = periodDts;
            periodDtsAssemble(level1, 2) = k1_ts;

            %% Part3: cal globalmean lamda_x({'ts','cloud','wv','ta','alb'})
            % regression to cal global lamda(lamda_global save inside for; lamda_globalAssemble save outside for)
            lamda_global = zeros(varUsedSize_vars, 1);
            X = [ones(size(varWorldMean(:, 1))) varWorldMean(:, 1)]; %dts

            for var_id = 1:varUsedSize_vars
                [b, bint, r, rint, stats] = regress(varWorldMean(:, var_id), X);
                lamda_global(var_id) = b(2);
            end

            lamda_globalAssemble(:, level1, skyLevel) = lamda_global; % save every model result

            %% Part4: cal every grid lamda_x({'ts','cloud','wv','ta','alb'})
            % grid regresssion
            lamda_grid = zeros(nlon, nlat, varUsedSize_vars);

            for var_id = 1:varUsedSize_vars

                for latNum = 1:nlat

                    for lonNum = 1:nlon
                        varUsed_temp1 = squeeze(squeeze(varUsed(lonNum, latNum, :, 1)));
                        varUsed_temp2 = squeeze(squeeze(varUsed(lonNum, latNum, :, var_id)));
                        X = [ones(size(varUsed_temp1)) varUsed_temp1];
                        [b, bint, r, rint, stats] = regress(varUsed_temp2, X);
                        lamda_grid(lonNum, latNum, var_id) = b(2);
                    end

                end

            end

            lamda_gridAssemble(:, :, :, level1, skyLevel) = lamda_grid;
            readme.meanOfRow = 'row1-5 ts,cloud,wv,ta,alb';
            save([dradFeedbackPath, 'global_vars.mat'], 'lonf', 'latf', 'time', 'plevf', 'readme', 'timeseries', 'modelname')
            save([dradFeedbackPath, ['lamda_global_', maskLandSW, '_', toaSfc{skyLevel}, '.mat']], 'lamda_global')
            save([dradFeedbackPath, ['lamda_grid_', maskLandSW, '_', toaSfc{skyLevel}, '.mat']], 'lamda_grid')

        end

        disp([level.model2{level1}, ' model is done!'])
        disp(' ')
    end

    % save assemble data
    lamda_globalAssembleMean = squeeze(mean(lamda_globalAssemble, 2));
    lamda_gridAssembleMean = squeeze(mean(lamda_gridAssemble, 4));
    readme.meanOfdim = '(lon,lat,vars,modelname,skylevel)';
    readme.meanOfPeriodDtsAssemble = '(periodDts,k1_ts)';
    save([dradFeedbackAssemblePath, 'global_vars.mat'], 'lonf', 'latf', 'time', 'plevf', 'readme', 'timeseries', 'modelname')
    save([dradFeedbackAssemblePath, ['lamda_globalAssemble_', maskLandSW, '.mat']], 'lamda_globalAssemble', 'lamda_globalAssembleMean')
    save([dradFeedbackAssemblePath, ['lamda_gridAssemble_', maskLandSW, '.mat']], 'lamda_gridAssemble', 'lamda_gridAssembleMean')
    save([dradFeedbackAssemblePath, ['periodDtsAssemble_', maskLandSW, '.mat']], 'periodDtsAssemble')

    %% Part5: regress and plot lamda_cloud vs period_dts
    lamda_cld_temp = squeeze(lamda_globalAssemble(2, :, :));
    periodDtsVs = periodDtsAssemble(:, 1);

    for skyLevel = 1:2

        lamda_cld = squeeze(lamda_cld_temp(:, skyLevel));
        X = [ones(size(lamda_cld)) lamda_cld];
        [b, bint, r, rint, stats] = regress(periodDtsVs, X);
        k1_cld = b(2);
        k2_cld = b(1);
        save([dradFeedbackAssemblePath, ['k_cld_', maskLandSW, '_', toaSfc{skyLevel}, '.mat']], 'k1_cld', 'k2_cld')
        %plot
        k = polyfit(lamda_cld, periodDtsVs, 1); % 一元一次拟合
        yfit = polyval(k, lamda_cld);
        yresid = periodDtsVs - yfit; %将残差值计算为有符号数的向量：
        SSresid = sum(yresid.^2); %计算残差的平方并相加，以获得残差平方和：
        SStotal = (length(periodDtsVs) - 1) * var(periodDtsVs); %通过将观测次数减 1(自由度) 再乘以 y 的方差，计算 y 的总平方和：
        rsq = 1 - SSresid / SStotal;%计算r2
        
        plot(lamda_cld, periodDtsVs, 'o', lamda_cld, yfit, '-')
        xlabel('\lambda_{cloud}');
        ylabel('\DeltaTsg');    
        title(['Experient:', level.time1{p_1}(1:end-1), ' level:', toaSfc{skyLevel}],'Interpreter', 'none');
        text(0.65, 0.95,['r2= ',num2str(rsq)],'units','normalized');    
        text(0.65, 0.9,['y=',num2str(k1_cld),'x+',num2str(k2_cld)],'units','normalized');    
        figurename = [figTestPath, level.time1{p_1}(1:end-1), '_', toaSfc{skyLevel}, '.png'];
        saveas(gcf, figurename)
        close gcf

    end

    disp([level.time1{p_1}, ' era is done!'])
    disp(' ')
    clear varUsed
end

t = toc; disp(t)
% plot grid figure
% figure
% pcolor(lon, lat, lamda_gridAssemble1');
% shading flat
% colorbar
% caxis([-2 2])
% hold on
% colormap(mycolor(100))
% plot(world_mapx, world_mapy, 'k')
% caxis([-3 3])
