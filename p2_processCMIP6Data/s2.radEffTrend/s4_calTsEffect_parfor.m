%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-15 10:44:34
% LastEditTime : 2020-11-08 15:53:43
% LastEditors  : LYC
% Description  : cal dRvars effect on Ts : according Ts=dRx/Kts
% FilePath     : /code/p2_processCMIP6Data/s2.radEffTrend/s4_calTsEffect_parfor.m
% Note         : only cal cloud, hus, ta, alb on sfc/toa, if need more vars, add later.
%%---------------------------------------------------------
clear; clc; tic;
nowpath = pwd;

% 开启并行环境
poolobj = gcp('nocreate'); % If no pool,  create new one.

if isempty(poolobj)
    MyPar = parpool(10);
else
    MyPar = gcp('nocreate');
    disp('Already initialized'); %说明并行环境已经启动。
end

%% 预选读取所有的路径
exm1 = 1; exm2 = 4;
mdl1 = 1; mdl2 = 'end';
esm1 = 1; esm2 = 'end';
esmPath_assmble = cell(2, 1);
mdlPath_assmble = esmPath_assmble;
exmNum_assmble = zeros(2, 1);
mdlNum_assmble = exmNum_assmble;
esmNum_assmble = exmNum_assmble;
esmCount = 0;

for exmNum = exm1:exm2%1 mean amip 2000; 2 mean amip 1980; 3 means ssp245, 4 means ssp370, 6 abrupt-4xCO2_150years
    % model parameters
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % experiment path (tLin:1740)
    inputPath = '/data1/liuyincheng/CMIP6-process/';
    exmPath_all = cell(1, length(Experiment));

    for varNum = 1:length(Experiment)
        exmPath_all{varNum} = fullfile(inputPath, level.time1{exmNum});
    end

    exmPath = exmPath_all{exmNum};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model
    if strcmp(mdl2, 'end')
        mdlR = length(level.model2);
    else
        mdlR=mdl2;
    end

    % mdlPath_assmble=cell((mdl2-mdl1+1)*(exm2-exm1+1),1);
    for mdlNum = mdl1:mdlR% model number
        % model path
        mdlPath = fullfile(exmPath, level.model2{mdlNum});
        % ensemble member path
        esmName = getPath_fileName(mdlPath, '.');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensemble member
        if strcmp(esm2, 'end')
            esmR = length(esmName);
        end

        for esmNum = esm1:esmR
            esmCount = esmCount + 1;
            esmPath = fullfile(mdlPath, esmName{esmNum, 1});
            esmPath_assmble{esmCount, 1} = esmPath;
            mdlPath_assmble{esmCount, 1} = mdlPath;
            exmNum_assmble(esmCount, 1) = exmNum;
            mdlNum_assmble(esmCount, 1) = mdlNum;
            esmNum_assmble(esmCount, 1) = esmNum;
        end

    end

end

parfor esmNum = 1:length(esmPath_assmble)
    esmFun(mdlPath_assmble{esmNum}, esmPath_assmble{esmNum}, exmNum_assmble(esmNum), mdlNum_assmble(esmNum), esmNum_assmble(esmNum));
end

delete(MyPar)

eval(['cd ', nowpath]);

t = toc; tickTok(t)

function [] = esmFun(mdlPath, esmPath, exmNum, mdlNum, esmNum)

    latRange = 88.75; % Latitude range
    lon1 = [2.5 357.5]; lat1 = [-latRange + 1 latRange - 1]; % world area
    toaSfc = {'toa', 'sfc'};
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);

    % ensemble member path
    esmName = getPath_fileName(mdlPath, '.');
    disp([level.model2{mdlNum}, ', ', level.time1{exmNum}, ', ', esmName{esmNum, 1}, ' ensemble start!'])

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
    varUsed = zeros(nlonf, nlatf, ntime, 5); % dR
    varKerUsed = varUsed; % dR/kernels

    for skyLevel = 1:2% 1 mean toa, 2 mean sfc\
        % load dRx and load to one var
        load([dradEffectPath, ['dradEffect_', toaSfc{skyLevel}, '_cld.mat']])%albEffect, husEffect, taEffect, mainEffect, totalEffect, tsEffect, talwEffect, taswEffect
        load([dradEffectPath, 'dR_cloud.mat'])%dR_cloud_toa
        load([dradEffectPath, ['dR_residual_cld_', toaSfc{skyLevel}, '.mat']])%dR_residual_cld_sfc/dR_residual_cld_toa

        if strcmp(toaSfc{skyLevel}, 'toa') == 1
            varUsed(:, :, :, 1) = dR_cloud_toa;
            varUsed(:, :, :, 5) = dR_residual_cld_toa;
        elseif strcmp(toaSfc{skyLevel}, 'sfc') == 1
            varUsed(:, :, :, 1) = dR_cloud_sfc;
            varUsed(:, :, :, 5) = dR_residual_cld_sfc;
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
        dTs_residual = squeeze(varKerUsed(:, :, :, 5));
        save([vsTsEffectPath, 'dTs_x_', toaSfc{skyLevel}, '.mat'], 'dTs_cloud', 'dTs_hus', 'dTs_ta', 'dTs_alb', 'dTs_residual')

        % cal trend
        [trendm_dTs_cld, trends_dTs_cld, trendyr_dTs_cld, ~, ~] = autoCalTrend(dTs_cloud, nlonf, nlatf, time.date, startmonth);
        [trendm_dTs_hus, trends_dTs_hus, trendyr_dTs_hus, ~, ~] = autoCalTrend(dTs_hus, nlonf, nlatf, time.date, startmonth);
        [trendm_dTs_ta, trends_dTs_ta, trendyr_dTs_ta, ~, ~] = autoCalTrend(dTs_ta, nlonf, nlatf, time.date, startmonth);
        [trendm_dTs_alb, trends_dTs_alb, trendyr_dTs_alb, ~, ~] = autoCalTrend(dTs_alb, nlonf, nlatf, time.date, startmonth);
        [trendm_dTs_residual, trends_dTs_residual, trendyr_dTs_residual, ~, ~] = autoCalTrend(dTs_residual, nlonf, nlatf, time.date, startmonth);
        % save
        save([vsTsEffectTrendPath, 'trend_dTs_x_', toaSfc{skyLevel}, '.mat'], ...
            'trendm_dTs_cld', 'trends_dTs_cld', 'trendyr_dTs_cld', ...
            'trendm_dTs_hus', 'trends_dTs_hus', 'trendyr_dTs_hus', ...
            'trendm_dTs_ta', 'trends_dTs_ta', 'trendyr_dTs_ta', ...
            'trendm_dTs_alb', 'trends_dTs_alb', 'trendyr_dTs_alb', ...
            'trendm_dTs_residual', 'trends_dTs_residual', 'trendyr_dTs_residual')
        clear -regexp ^trendm ^trends ^trendyr
    end
    clear varKerUsed varUsed 

    disp([level.model2{mdlNum}, ', ', level.time1{exmNum}, ', ', esmName{esmNum, 1}, ' ensemble is done!'])
    disp(' ')

end
