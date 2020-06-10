%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-09 15:52:00
% LastEditTime : 2020-06-10 20:44:32
% LastEditors  : LYC
% Description  :
% FilePath     : /Research/p2_processCMIP6Data/s4.radFeedbakIntensity/s1_calRadFeedIntens.m
%
%%---------------------------------------------------------
% 计算全球各种量的反馈强度
% 云, 温度 水汽 反照率
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
    inputPath = ['/data1/liuyincheng/cmip6-process/', level.time1{p_1}]; %/data1/liuyincheng/cmip6-process/2000-2014/
    % model loop
    k_grid_assmble = zeros(144, 72, length(level.model2));
    lamda_grobleAssmble = zeros(1, length(level.model2));

    for level1 = 1:length(level.model2)
        % load dvars
        dvarsPath = [inputPath, level.model2{level1}, '/', level.process3{2}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
        varsPath = [inputPath, level.model2{level1}, '/', level.process3{1}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/rawdata
        dradEffectPath = [inputPath, level.model2{level1}, '/', level.process3{6}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/

        for skyLevel = 1:2% 1 mean toa, 2 mean sfc
            % vars.lamda={'cloud','wv','ta','alb'};
            load([dvarsPath, 'global_vars.mat'])% latf lonf time plevf readme
            load([dvarsPath, 'dts.mat'])%
            load([dradEffectPath, ['dradEfect_', toaSfc{skyLevel}, '_cld.mat']])%albEffect, husEffect, taEffect, mainEffect, totalEffect, tsEffect, wvlwEffect, wvswEffect
            load([dradEffectPath, ['dR_cloud_', toaSfc{skyLevel}, '.mat']])%dR_cloud_toa

            % regrid dts to 144x72(unite grids)
            lat = 88.75:-2.5:-88.75; nlat = length(lat);
            lon = lonf; nlon = length(lon);
            nlonf = length(lonf); nlatf = length(latf);
            [Xlon, Ylat, Ttime] = meshgrid(lat, lon, time);
            [Xlonf, Ylatf, Ttimef] = meshgrid(latf, lonf, time);
            dts = interp3(Xlonf, Ylatf, Ttimef, dts, Xlon, Ylat, Ttime);
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
            varUsedSize_time=varUsedSize(3);
            varUsedSize_vars=varUsedSize(4);
            % mask(world land)
            if strcmp(maskLandSW, 'maskLand') == 1
                for var_id = 1:varUsedSize_vars
                    tempVar = squeeze(varUsed(:, :, :, var_id));
                    [tempVar] = maskLand(tempVar, lat, p_3, -p_3, areaNum);
                    varUsed(:, :, :, var_id) = tempVar;
                end
            end

            % grobal Zonal weighted average (time, vars)
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

            % regression to cal grobal lamda{'cloud','wv','ta','alb'}
            lamda_grobal=zeros(varUsedSize_vars,1);
            lamda_grobleAssmble =zeros(varUsedSize_vars,length(level.model2));
            X = [ones(size(varWorldMean(:, 1))) varWorldMean(:, 1)];%dts
            for var_id = 1:varUsedSize_vars
                [b, bint, r, rint, stats] = regress(varWorldMean(:, var_id), X);
                lamda_grobal(var_id) = b(2);

            end
            lamda_grobleAssmble(level1) = lamda_grobal;


            % grid regress
            k_grid = zeros(nlon, nlat);
            for latNum = 1:nlat

                for lonNum = 1:nlon
                    varUsed_temp1 = squeeze(squeeze(varUsed(lonNum, latNum, :, 1)));
                    varUsed_temp2 = squeeze(squeeze(varUsed(lonNum, latNum, :, 2)));
                    X = [ones(size(varUsed_temp1)) varUsed_temp1];
                    [b, bint, r, rint, stats] = regress(varUsed_temp2, X);
                    k_grid(lonNum, latNum) = b(2);
                end

            end

            k_grid_assmble(:, :, level1) = k_grid;
        end

        k_grid_assmble1 = mean(k_grid_assmble, 3);
        lamda_grobleAssmble1 = mean(lamda_grobleAssmble);
    end

end

% plot grid figure
% figure
% pcolor(lon, lat, k_grid_assmble1');
% shading flat
% colorbar
% caxis([-2 2])
% hold on
% colormap(mycolor(100))
% plot(world_mapx, world_mapy, 'k')
% caxis([-3 3])
