% 计算全球各种量的反馈值
% 温度 水汽 反照率, 云
clear; clc; tic;
% pause(9999)
% model settings
load('/home/lyc/lib/tools/matlab/map/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/lyc/lib/tools/matlab/map/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/lyc/lib/tools/matlab/map/02.world_map/mat_file/correct_worldmap.mat')% ????????????????word_mapx(:),word_mapy(:)
load('/home/lyc/lib/tools/matlab/map/01.china_map/mat_file/mask14472.mat')
areaNum = 1; % world land
p_3 = 88.75; % Latitude range
lon1 = [2.5 357.5]; lat1 = [-p_3 + 1 p_3 - 1]; % world area

for p_1 = 5:5%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    nowpath = pwd;
    %CAMS-CSM1-0 didn't have sfc clear sky radiation, delete it
    [readme, Experiment, level, tLin, mPlev, vars] = modelParameters(p_1);

    % inputPath
    inputPath = ['/data1/liuyincheng/cmip6-process/', level.time1{p_1}]; %/data1/liuyincheng/cmip6-process/2000-2014/
    % model loop
    k_grid_assmble = zeros(144, 72, length(level.model2));
    k_groble_assmble = zeros(1, length(level.model2));

    for level1 = 1:length(level.model2)
        % load dvars
        dvarsPath = [inputPath, level.model2{level1}, '/', level.process3{2}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
        varsPath = [inputPath, level.model2{level1}, '/', level.process3{1}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/rawdata
        dradEffectPath = [inputPath, level.model2{level1}, '/', level.process3{6}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/

        load([dvarsPath, 'global_vars.mat'])% latf lonf time plevf readme
        load([dvarsPath, 'dts.mat'])%
        load([dradEffectPath, 'dradEfect_toa_cld.mat'])%albEffect, husEffect, mainEffect, taEffect, totalEffect, tsEffect, wvlwEffect, wvswEffect
        load([dradEffectPath, 'dR_cloud_toa.mat'])%dR_cloud_toa

        %regrid 144x72(unite grids)
        lat = 88.75:-2.5:-88.75; nlat = length(lat);
        lon = lonf; nlon = length(lon);
        nlonf = length(lonf); nlatf = length(latf);
        [Xlon, Ylat, Ttime] = meshgrid(lat, lon, time);
        [Xlonf, Ylatf, Ttimef] = meshgrid(latf, lonf, time);
        dts = interp3(Xlonf, Ylatf, Ttimef, dts, Xlon, Ylat, Ttime);
        % 1 mean
        tEffect = taEffect + tsEffect;
        varUsed(:, :, :, 1) = dts; %144x72x1800x2
        varUsed(:, :, :, 2) = dR_cloud_toa;
        % mask(world land)
        % for var_id = 1:2
        %     tempVar=squeeze(varUsed(:,:,:,var_id));
        %     [tempVar] = maskLand(tempVar, lat, p_3, -p_3, areaNum);
        %     varUsed(:,:,:,var_id)=tempVar;
        % end
        % Zonal weighted average(vars, time)
        jiaquan = cosd(lat);
        wei = ones(144, 72); %格点纬度加权

        for latiNum = 1:72
            wei(:, latiNum) = wei(:, latiNum) * jiaquan(latiNum); %格点相对大小
        end

        varUsed = varUsed .* wei;
        varUsedSize = size(varUsed);
        varWorldMean = zeros(varUsedSize(3), varUsedSize(4));

        for var_id = 1:varUsedSize(4)

            for timeNum = 1:length(varUsed)
                varWorldMean(timeNum, var_id) = nansum(nansum(varUsed(:, :, timeNum, var_id))) / nansum(nansum(wei));
            end

        end

        % regression
        % grobal mean
        X = [ones(size(varWorldMean(:, 1))) varWorldMean(:, 1)];
        [b, bint, r, rint, stats] = regress(varWorldMean(:, 2), X);
        k_grobal = b(2);
        k_groble_assmble(level1) = k_grobal;
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
    k_groble_assmble1 = mean(k_groble_assmble);
    % plot grid figure
    figure
    pcolor(lon, lat, k_grid_assmble1');
    shading flat
    colorbar
    caxis([-2 2])
    hold on
    colormap(mycolor(100))
    plot(world_mapx, world_mapy, 'k')
    caxis([-3 3])

end
