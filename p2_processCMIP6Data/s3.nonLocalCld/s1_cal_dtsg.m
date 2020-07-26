%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-07-22 20:44:45
% LastEditTime : 2020-07-22 20:51:56
% LastEditors  : LYC
% Description  :
% FilePath     : /code/p2_processCMIP6Data/s3.nonLocalCld/s1_cal_dtsg.m
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
figTestPath = '/data1/liuyincheng/cmip6-process/z_assembleData/figTest/';
auto_mkdir(figTestPath)
% load lamda_cloud
load('/data1/liuyincheng/cmip6-process/z_globalVar/lamda_cloud.mat')
lon_k = 0:2.5:357.5; nlonk = length(lon_k);
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);

for exmNum = 1:4%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    nowpath = pwd;
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);

    % exmPath
    exmPath = ['/data1/liuyincheng/cmip6-process/', level.time1{exmNum}]; %/data1/liuyincheng/cmip6-process/amip_2000-2014/
    dvarsPath = [exmPath, level.model2{1}, '/r1i1p1f1/', level.process3{2}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
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
            dvarsPath = fullfile(esmPath, level.process3{2}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
            dvarsTrendPath = fullfile(esmPath, level.process3{3}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly_trend
            varsPath = fullfile(esmPath, level.process3{1}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/rawdata
            dradEffectPath = fullfile(esmPath, level.process3{6}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
            dnonCloudPath = fullfile(esmPath, level.process3{8}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/non_localCld/
            auto_mkdir(dnonCloudPath)
            load([dvarsPath, 'global_vars.mat'])% lat_k lon_k time plev_k readme
            load([dvarsPath, 'dts.mat'])%
            % regrid dts to 144x72(unite grids)
            dts = autoRegrid3(lon_k, lat_k, time.date, dts, lon_f, lat_f, time.date);
            %% Part1: cal dtsg series and DeltaTsg
            % 1.1 cal dtsg
            % land mask (optional)
            if strcmp(maskLandSW, 'maskLand') == 1

                for var_id = 1:varUsedSize_vars
                    [dts] = maskLand(dts, lat_f, latRange, -latRange, areaNum);
                end

            end

            % global Zonal weighted average (time, vars)
            jiaquan = cosd(lat_f);
            wei = ones(144, 72); %格点纬度加权

            for latiNum = 1:72
                wei(:, latiNum) = wei(:, latiNum) * jiaquan(latiNum); %格点相对大小
            end

            tsUsed = dts .* wei;
            dtsg = zeros(ntime, 1); % global mean temperture anomaly time series

            for timeNum = 1:ntime
                dtsg(timeNum) = nansum(nansum(tsUsed(:, :, timeNum))) / nansum(nansum(wei));
            end

            readmeDtsg = 'global Zonal weighted average surface temperture';
            save([dvarsPath, ['dtsg_', maskLandSW, '.mat']], 'dtsg', 'readmeDtsg')

            % 1.2 cal period DeltaTsg (whole time line)
            X = [ones(size(time.date)) time.date]; %dts
            [b, bint, r, rint, stats] = regress(dtsg, X);
            trend_dtsg = b(2); % slop of time series K/时间的单位
            trendUnit = 'K/day';
            DeltaTsg = trend_dtsg * (time.date(end) - time.date(1));
            DeltaTsgMeaning = 'whole period ts change';
            readme.maskStatement = 'nomask mean land+sea; mask land mean no sea';
            save([dnonCloudPath, 'global_vars.mat'], 'lon_k', 'lat_k', 'time', 'plev_k', 'readme', 'timeseries', 'modelname')
            save([dvarsTrendPath, 'trend_dtsg.mat'], 'trend_dtsg', 'trendUnit')
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
