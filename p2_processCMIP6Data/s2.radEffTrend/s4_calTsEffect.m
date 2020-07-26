%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-15 10:44:34
% LastEditTime : 2020-07-23 15:34:11
% LastEditors  : LYC
% Description  : cal dRvars effect on Ts : according Ts=dRx/Kts
% FilePath     : /code/p2_processCMIP6Data/s2.radEffTrend/s4_calTsEffect.m
% Note         : only cal cloud, hus, ta, alb on sfc/toa, if need more vars, add later.
%%---------------------------------------------------------
clear; clc; tic;
latRange = 88.75; % Latitude range
lon1 = [2.5 357.5]; lat1 = [-latRange + 1 latRange - 1]; % world area
toaSfc = {'toa', 'sfc'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment
for exmNum = 1:3%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    nowpath = pwd;
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % exmPath
    exmPath = ['/data1/liuyincheng/cmip6-process/', level.time1{exmNum}]; %/data1/liuyincheng/cmip6-process/amip_2000-2014/
    disp([level.time1{exmNum}, ' era start!'])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model
    for mdlNum = 1:length(level.model2)
        % model path
        mdlPath = fullfile(exmPath, level.model2{mdlNum});
        eval(['cd ', mdlPath]);
        disp(' ')
        disp([level.model2{mdlNum}, ' model start!'])
        % ensemble member path
        esmName = getPath_fileName(mdlPath, '.');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensemble member
        for esmNum = 1:length(esmName)
            esmPath = fullfile(mdlPath, esmName{esmNum, 1});
            eval(['cd ', esmPath]);
            % input and outputpath
            dvarsPath = fullfile(esmPath, level.process3{2}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
            dvarsTrendPath = fullfile(esmPath, level.process3{3}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly_trend
            kernelPath = fullfile(esmPath, level.process3{5}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/kernelsCal
            varsPath = fullfile(esmPath, level.process3{1}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/rawdata
            dradEffectPath = fullfile(esmPath, level.process3{6}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
            dnonLocalCldPath = fullfile(esmPath, level.process3{8}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/non_localCld/
            % outPath
            vsTsEffectPath = fullfile(esmPath, level.process3{9}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/vsTsEffect/
            vsTsEffectTrendPath = fullfile(esmPath, level.process3{10}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/vsTsEffectTrend/
            auto_mkdir(vsTsEffectPath)
            auto_mkdir(vsTsEffectTrendPath)

            load([dradEffectPath, 'global_vars.mat'])% lat_f lon_f time plev_k readme
            nlatf = length(lat_f); nlonf = length(lon_f); ntime = length(time.date);
            startmonth = 1;
            varUsed = zeros(nlonf, nlatf, ntime, 4); % dR
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
                load([kernelPath, '/global_vars.mat']); % lat_k,lon_k,time
                % regrid ts_lwkernel to 144x72(unite grids)
                kernelTime = 1:12;
                ts_lwkernel = autoRegrid3(lon_k, lat_k, kernelTime, ts_lwkernel, lon_f, lat_f, kernelTime);
                % extend to the whole time series
                startmonth = 1;
                if startmonth ~= 1
                    tempK = zeros(nlonf, nlatf, 12);
                    tempK(1:12 - startmonth + 1) = ts_lwkernel(:, :, startmonth:12);
                    tempK(12 - startmonth + 2:12) = ts_lwkernel(:, :, 1:startmonth - 1);
                    ts_lwkernel = tempK;
                end

                nyear = ntime / 12;
                ts_lwkernelSeries = repmat(ts_lwkernel, [1 1 nyear]);
                % cal dRs/kernel_ts
                for varNum = 1:varUsedNum
                    varKerUsed(:, :, :, varNum) = squeeze(varUsed(:, :, :, varNum)) ./ -ts_lwkernelSeries;
                end

                % save...
                dTs_cloud = squeeze(varKerUsed(:, :, :, 1));
                dTs_hus = squeeze(varKerUsed(:, :, :, 2));
                dTs_ta = squeeze(varKerUsed(:, :, :, 3));
                dTs_alb = squeeze(varKerUsed(:, :, :, 4));
                save([vsTsEffectPath, 'dTs_x_', toaSfc{skyLevel}, '.mat'], 'dTs_cloud', 'dTs_hus', 'dTs_ta', 'dTs_alb')

                % cal trend
                [trendm_dTs_cld, trends_dTs_cld, trendyr_dTs_cld, ~, ~] = autoCalTrend(dTs_cloud, nlonf, nlatf, time.date, startmonth);
                [trendm_dTs_hus, trends_dTs_hus, trendyr_dTs_hus, ~, ~] = autoCalTrend(dTs_hus, nlonf, nlatf, time.date, startmonth);
                [trendm_dTs_ta, trends_dTs_ta, trendyr_dTs_ta, ~, ~] = autoCalTrend(dTs_ta, nlonf, nlatf, time.date, startmonth);
                [trendm_dTs_alb, trends_dTs_alb, trendyr_dTs_alb, ~, ~] = autoCalTrend(dTs_alb, nlonf, nlatf, time.date, startmonth);
                % save
                save([vsTsEffectTrendPath, 'trend_dTs_x_', toaSfc{skyLevel}, '.mat'], ...
                    'trendm_dTs_cld', 'trends_dTs_cld', 'trendyr_dTs_cld', ...
                    'trendm_dTs_hus', 'trends_dTs_hus', 'trendyr_dTs_hus', ...
                    'trendm_dTs_ta', 'trends_dTs_ta', 'trendyr_dTs_ta', ...
                    'trendm_dTs_alb', 'trends_dTs_alb', 'trendyr_dTs_alb')
            end
            disp([esmName{esmNum,1}, ' ensemble is done!'])
        end
        disp([level.model2{mdlNum}, ' model is done!'])
        disp(' ')
    end
    disp([level.time1{exmNum}, ' era is done!'])
    disp(' ')
end

t = toc; disp(t)
