%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-09 15:52:00
% LastEditTime : 2020-07-22 20:40:32
% LastEditors  : LYC
% Description  : 计算lamda_cloud 和 deltaTsg的线性关系(k1和k2)
% FilePath     : /code/p2_processCMIP6Data/s3.nonLocalCld/s1p1_cal1_nonLocalCld_Tsgk1k2.m
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
    dvarsPath = [exmPath, level.model2{1}, '/', level.process3{2}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
    load([dvarsPath, 'global_vars.mat'])% lat_k lon_k time plev_k readme
    ntime = length(time.date);
    DeltaTsgAssemble = zeros(length(level.model2), 2); % global mean temperture change during whole period
    lamdaCloudAssemble = zeros(length(level.model2), 1);
    dnonCloudAssemblePath = [exmPath, 'z_assembleData/', level.process3{8}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/non_localCld/
    auto_mkdir(dnonCloudAssemblePath)
    countm=1;
    % model loop
    for mdlNum = 1:length(level.model2)
        % model path
        mdlPath = fullfile(exmPath, level.model2{mdlNum});
        eval(['cd ', mdlPath]);
        disp(' ')
        disp([level.model2{mdlNum}, ' model start!'])

        % path
        dvarsPath = [exmPath, level.model2{mdlNum}, '/', level.process3{2}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
        dvarsTrendPath = [exmPath, level.model2{mdlNum}, '/', level.process3{3}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly_trend
        varsPath = [exmPath, level.model2{mdlNum}, '/', level.process3{1}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/rawdata
        dradEffectPath = [exmPath, level.model2{mdlNum}, '/', level.process3{6}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
        dnonCloudPath = [exmPath, level.model2{mdlNum}, '/', level.process3{8}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/non_localCld/
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
        DeltaTsgAssemble(mdlNum, 1) = DeltaTsg;
        DeltaTsgAssemble(mdlNum, 2) = trend_dtsg;
        if isempty(lamda_cloud.modelVar(strcmp(lamda_cloud.modelName, level.model2{mdlNum}) == 1))==1
            lamdaCloudAssemble(mdlNum) = nan;
        else
            lamdaCloudAssemble(mdlNum) = lamda_cloud.modelVar(strcmp(lamda_cloud.modelName, level.model2{mdlNum}) == 1);
            existModelName{countm} = level.model2{mdlNum}; % record used model name 
            countm=countm+1;
        end

        disp([level.model2{mdlNum}, ' model is done!'])
        disp(' ')
    end
    % save assemble data
    readme.meanOfDeltaTsgAssemble = '(DeltaTsg,trend_dtsg)';
    save([dnonCloudAssemblePath, 'global_vars.mat'], 'lon_k', 'lat_k', 'time', 'plev_k', 'readme', 'timeseries', 'modelname')
    save([dnonCloudAssemblePath, ['DeltaTsgAssemble_', maskLandSW, '.mat']], 'DeltaTsgAssemble')
    save([dnonCloudAssemblePath, 'lamdaCloudAssemble.mat'], 'lamdaCloudAssemble','existModelName')
    %% Part2: cal k1, k2 about lamda_cloud and DeltaTsAssemble(DeltaTsAssemble=k1*lamda_cloud+k2) and plot
    DeltaTsgVs = DeltaTsgAssemble(:, 1);

    DeltaTsgVs=DeltaTsgVs(~isnan(lamdaCloudAssemble));
    lamdaCloudAssemble=lamdaCloudAssemble(~isnan(lamdaCloudAssemble));
    X = [ones(size(lamdaCloudAssemble)) lamdaCloudAssemble];
    [b, bint, r, rint, stats] = regress(DeltaTsgVs, X);
    k1_cld = b(2);
    k2_cld = b(1);
    save([dnonCloudAssemblePath, 'k_cld.mat'], 'k1_cld', 'k2_cld','existModelName')
    clear existModelName
    %plot
    k = polyfit(lamdaCloudAssemble, DeltaTsgVs, 1); % 一元一次拟合
    yfit = polyval(k, lamdaCloudAssemble);
    yresid = DeltaTsgVs - yfit; %将残差值计算为有符号数的向量：
    SSresid = sum(yresid.^2); %计算残差的平方并相加，以获得残差平方和：
    SStotal = (length(DeltaTsgVs) - 1) * var(DeltaTsgVs); %通过将观测次数减 1(自由度) 再乘以 y 的方差，计算 y 的总平方和：
    rsq = 1 - SSresid / SStotal; %计算r2

    plot(lamdaCloudAssemble, DeltaTsgVs, 'o', lamdaCloudAssemble, yfit, '-')
    xlabel('\lambda_{cloud}');
    ylabel('\DeltaTsg');
    title(['Experient:', level.time1{exmNum}(1:end - 1)], 'Interpreter', 'none');
    text(0.65, 0.95, ['r2= ', num2str(rsq)], 'units', 'normalized');
    text(0.65, 0.9, ['y=', num2str(k1_cld), 'x+', num2str(k2_cld)], 'units', 'normalized');
    figurename = [figTestPath, level.time1{exmNum}(1:end - 1), '.png'];
    saveas(gcf, figurename)
    close gcf
    
    disp([level.time1{exmNum}, ' era is done!'])
    disp(' ')

end

t = toc; disp(t)
