%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-07-29 09:58:38
% LastEditTime : 2021-04-02 10:46:11
% LastEditors  : Please set LastEditors
% Description  : 相比旧方法, 这里计算地表平均温度时不直接用dts来计算, 而是从ts出发, 先算出总的纬向加权, 然后在进行求anomaly等一系列操作
% FilePath     : /code/p2_processCMIP6Data/s4.nonLocalCld/s1_cal_dtsg_newMehod.m
%
%%---------------------------------------------------------

% 云, 温度 水汽 反照率
% 注意此脚本保存应该都要加上mask类型的后缀('nomask'or'maskLand')
clear; clc; tic;
% pause(9999)
% model settings
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')

latRange = 88.75; % Latitude range
lon1 = [2.5 357.5]; lat1 = [-latRange + 1 latRange - 1]; % world area
toaSfc = {'toa', 'sfc'};
maskLandSW = 'nomask'; %{'nomask', 'maskLand'};
areaNum = 1; % world land
figTestPath = '/data1/liuyincheng/CMIP6-process/z_assembleData/figTest/';
auto_mkdir(figTestPath)
% load lamda_cloud
load('/data1/liuyincheng/CMIP6-process/z_globalVar/lamda_cloud.mat')
lon_k = 0:2.5:357.5; nlonk = length(lon_k);
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);
startMonth = 1;

for exmNum = 1:4%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    nowpath = pwd;
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);

    % exmPath
    exmPath = fullfile('/data1/liuyincheng/CMIP6-process', level.time1{exmNum}); %/data1/liuyincheng/CMIP6-process/amip_2000-2014/
    dvarsPath = fullfile(exmPath, level.model2{1}, 'r1i1p1f1', level.process3{2}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/anomaly
    load([dvarsPath, 'global_vars.mat'])% lat_k lon_k time plev_k readme
    ntime = length(time.date);
    % model loop
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
            % path
            dvarsPath = fullfile(esmPath, level.process3{2}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/anomaly
            dvarsTrendPath = fullfile(esmPath, level.process3{3}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/anomaly_trend
            varsPath = fullfile(esmPath, level.process3{1}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/rawdata
            dradEffectPath = fullfile(esmPath, level.process3{6}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/radEffect/
            dnonCloudPath = fullfile(esmPath, level.process3{8}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/non_localCld/
            auto_mkdir(dnonCloudPath)
            load([varsPath, 'global_vars.mat'])% lat_k lon_k time plev_k readme
            load([varsPath, 'ts.mat'])%
            % regrid ts to 144x72(unite grids)
            ts = autoRegrid3(lon_k, lat_k, time.date, ts, lon_f, lat_f, time.date);
            %% Part1: cal tsg series and DeltaTsg
            % 1.1 cal tsg
            % land mask (optional)
            if strcmp(maskLandSW, 'maskLand') == 1

                for var_id = 1:varUsedSize_vars
                    [ts] = maskLand(ts, lat_f, latRange, -latRange, areaNum);
                end

            end

            % global Zonal weighted average (time, vars)
            jiaquan = cosd(lat_f);
            wei = ones(144, 72); %格点纬度加权

            for latiNum = 1:72
                wei(:, latiNum) = wei(:, latiNum) * jiaquan(latiNum); %格点相对大小
            end

            tsUsed = ts .* wei;
            tsg = zeros(ntime, 1); % global mean temperture anomaly time series

            for timeNum = 1:ntime
                tsg(timeNum) = nansum(nansum(tsUsed(:, :, timeNum))) / nansum(nansum(wei));
            end

            readmeDtsg = 'global Zonal weighted average surface temperture';
            save([varsPath, ['tsg_', maskLandSW, '.mat']], 'tsg', 'readmeDtsg')
            % cal dtsg
            [dtsg, tsg_clim] = monthlyAnomaly(144, 73, 73, time.date, tsg, startMonth); %[nlongitude,nlatitude,time,var,startMonth]
            save([dvarsPath, ['dtsg_', maskLandSW, '.mat']], 'dtsg', 'tsg_clim')

            % 1.2 cal period DeltaTsg (whole time line)
            % method_ori
            % X = [ones(size(time.date)) time.date]; %ts
            % [b, bint, r, rint, stats] = regress(tsg, X);
            % trend_dtsg = b(2); % slop of time series K/时间的单位
            % trendUnit = 'K/day';
            % DeltaTsg = trend_dtsg * (time.date(end) - time.date(1));

            % method_new
            trends_dtsg = zeros(4, 1);
            [~, trendm_dtsg, ~, ~] = detrend_yc(dtsg, time.date, startMonth);
            % seasons mean trend (unite:per day)
            for seaNum = 1:4
                if seaNum ~= 4
                    trends_dtsg(seaNum) = squeeze(mean(trendm_dtsg(seaNum * 3:seaNum * 3 + 2,1)));
                elseif seaNum == 4
                    seaNum1 = [1, 2, 12];
                    trends_dtsg(seaNum) = squeeze(mean(trendm_dtsg(seaNum1,1)));
                end
            end
            trendyr_dtsg = squeeze(mean(trendm_dtsg(:,1)));
            trendUnit = 'K/day';
            DeltaTsg = trendyr_dtsg * (time.date(end) - time.date(1));

            DeltaTsgMeaning = 'whole period ts change';
            readme.maskStatement = 'nomask mean land+sea; mask land mean no sea';
            save([dnonCloudPath, 'global_vars.mat'], 'lon_f', 'lat_f', 'time', 'plev_k', 'readme', 'timeseries', 'modelname')
            save([dvarsTrendPath, 'trend_dtsg.mat'], 'trendyr_dtsg','trends_dtsg','trendm_dtsg', 'trendUnit')
            save([dvarsTrendPath, ['DeltaTsg_', maskLandSW, '.mat']], 'DeltaTsg', 'DeltaTsgMeaning')
            disp([esmName{esmNum, 1}, ' ensemble is done!'])
        end

        disp([level.model2{mdlNum}, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{exmNum}, ' era is done!'])
    disp(' ')

end

t = toc; disp(t)
