%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-15 10:44:34
% LastEditTime : 2020-06-15 16:09:07
% LastEditors  : LYC
% Description  : cal dRvars effect on Ts : according Ts=dRx/Kts
% FilePath     : /Research/p2_processCMIP6Data/s3.radEffTrend/s4p1_calTsEffect.m
% Note         : only cal cloud, hus, ta, alb on sfc/toa, if need more vars, add later.
%%---------------------------------------------------------
clear; clc; tic;
p_3 = 88.75; % Latitude range
lon1 = [2.5 357.5]; lat1 = [-p_3 + 1 p_3 - 1]; % world area
toaSfc = {'toa', 'sfc'};

for p_1 = [2 4]%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    nowpath = pwd;
    [readme, Experiment, level, tLin, mPlev, vars] = modelParameters(p_1);
    % inputPath
    inputPath = ['/data1/liuyincheng/cmip6-process/', level.time1{p_1}]; %/data1/liuyincheng/cmip6-process/amip_2000-2014/
    dvarsPath = [inputPath, level.model2{1}, '/', level.process3{6}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
    load([dvarsPath, 'global_vars.mat'])% lat lon time plevf readme
    nlat = length(lat); nlon = length(lon); ntime = length(time);
    startmonth = 1;

    for level1 = 1:length(level.model2)
        % data path
        dvarsPath = [inputPath, level.model2{level1}, '/', level.process3{2}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
        dvarsTrendPath = [inputPath, level.model2{level1}, '/', level.process3{3}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly_trend
        kernelPath = [inputPath, level.model2{level1}, '/', level.process3{5}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/kernelsCal
        varsPath = [inputPath, level.model2{level1}, '/', level.process3{1}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/rawdata
        dradEffectPath = [inputPath, level.model2{level1}, '/', level.process3{6}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
        dnonLocalCldPath = [inputPath, level.model2{level1}, '/', level.process3{8}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/non_localCld/
        vsTsEffectPath = [inputPath, level.model2{level1}, '/', level.process3{9}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/vsTsEffect/
        vsTsEffectTrendPath = [inputPath, level.model2{level1}, '/', level.process3{10}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/vsTsEffectTrend/
        auto_mkdir(vsTsEffectPath)
        auto_mkdir(vsTsEffectTrendPath)

        varUsed = zeros(nlon, nlat, ntime, 4); % dR
        varKerUsed = varUsed; % dR/kernels

        for skyLevel = 1:2% 1 mean toa, 2 mean sfc\
            % load dRx and load to one var
            load([dradEffectPath, ['dradEfect_', toaSfc{skyLevel}, '_cld.mat']])%albEffect, husEffect, taEffect, mainEffect, totalEffect, tsEffect, talwEffect, taswEffect
            load([dradEffectPath, ['dR_cloud_', toaSfc{skyLevel}, '.mat']])%dR_cloud_toa

            if strcmp(toaSfc{skyLevel}, 'toa') == 1
                varUsed(:, :, :, 1) = dR_cloud_toa;
            elseif strcmp(toaSfc{skyLevel}, 'sfc') == 1
                varUsed(:, :, :, 1) = dR_cloud_sfc;
            end

            varUsed(:, :, :, 2) = husEffect;
            varUsed(:, :, :, 3) = taEffect;
            varUsed(:, :, :, 4) = albEffect;
            varUsedSize = size(varUsed);
            varUsedNum = varUsedSize(4);
            %
            % load ts_kernel
            load([kernelPath, '/kernels_', toaSfc{skyLevel}, '_cld'], 'ts_lwkernel');
            load([kernelPath, '/global_vars.mat']); % latf,lonf,time
            % regrid ts_lwkernel to 144x72(unite grids)
            kernelTime = 1:12;
            [Xlon, Ylat, Ttime] = meshgrid(lat, lon, kernelTime);
            [Xlonf, Ylatf, Ttimef] = meshgrid(latf, lonf, kernelTime);
            ts_lwkernel = interp3(Xlonf, Ylatf, Ttimef, ts_lwkernel, Xlon, Ylat, Ttime);
            % extend to the whole time series
            startmonth = 1;

            if startmonth ~= 1
                tempK = zeros(nlon, nlat, 12);
                tempK(1:12 - startmonth + 1) = ts_lwkernel(:, :, startmonth:12);
                tempK(12 - startmonth + 2:12) = ts_lwkernel(:, :, 1:startmonth - 1);
                ts_lwkernel = tempK;
            end

            nyear = ntime / 12;
            ts_lwkernelSeries = repmat(ts_lwkernel, [1 1 nyear]);
            % cal dRs/kernel_ts
            for varNum = 1:varUsedNum
                varKerUsed(:, :, :, varNum) = squeeze(varUsed(:, :, :, varNum)) ./ ts_lwkernelSeries;
            end

            % save...
            dTs_cloud = squeeze(varKerUsed(:, :, :, 1));
            dTs_hus = squeeze(varKerUsed(:, :, :, 2));
            dTs_ta = squeeze(varKerUsed(:, :, :, 3));
            dTs_alb = squeeze(varKerUsed(:, :, :, 4));
            save([vsTsEffectPath, 'dTs_x_', toaSfc{skyLevel}, '.mat'], 'dTs_cloud', 'dTs_hus', 'dTs_ta', 'dTs_alb')

            % cal trend
            [trendm_dTs_cld, trends_dTs_cld, trendyr_dTs_cld, ~, ~] = autoCalTrend(dTs_cloud, nlon, nlat, time, startmonth);
            [trendm_dTs_hus, trends_dTs_hus, trendyr_dTs_hus, ~, ~] = autoCalTrend(dTs_hus, nlon, nlat, time, startmonth);
            [trendm_dTs_ta, trends_dTs_ta, trendyr_dTs_ta, ~, ~] = autoCalTrend(dTs_ta, nlon, nlat, time, startmonth);
            [trendm_dTs_alb, trends_dTs_alb, trendyr_dTs_alb, ~, ~] = autoCalTrend(dTs_alb, nlon, nlat, time, startmonth);
            % save
            save([vsTsEffectTrendPath, 'trend_dTs_x_', toaSfc{skyLevel}, '.mat'], ...
                'trendm_dTs_cld', 'trends_dTs_cld', 'trendyr_dTs_cld', ...
                'trendm_dTs_hus', 'trends_dTs_hus', 'trendyr_dTs_hus', ...
                'trendm_dTs_ta', 'trends_dTs_ta', 'trendyr_dTs_ta', ...
                'trendm_dTs_alb', 'trends_dTs_alb', 'trendyr_dTs_alb')

        end

        % load dtsg

        disp([level.model2{level1}, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{p_1}, ' era is done!'])
    disp(' ')

end

t = toc; disp(t)
