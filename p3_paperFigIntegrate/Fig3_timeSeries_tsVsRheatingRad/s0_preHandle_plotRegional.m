%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-08-31 17:00:15
% LastEditTime : 2021-04-08 20:51:14
% LastEditors  : Please set LastEditors
% Description  : MME result of 时间序列图
% FilePath     : /code/p3_paperFigIntegrate/Fig3_timeSeries_tsVsRheatingRad/s0_preHandle_plotRegional_fig3.m
% note : 统一用startmonth=3 开始计算
%%---------------------------------------------------------
clc; clear; tic;
nowpath = pwd;
% load mask map
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat') % load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat') % load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')

toaSfc = {'toa', 'sfc'};
lon_k = 0:2.5:357.5; nlonk = length(lon_k); % kernel lat lon
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CMIP6 part
% the number of MME experiment
MMENum = 'MME1';
% experiment
for exmNum = 1:2 %1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    %CAMS-CSM1-0 didn't have sfc clear sky radiation, delete it
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % MMEPath
    MMEPath = ['/data1/liuyincheng/CMIP6-process/z_ensembleMean/MME1/', level.time1{exmNum}]; %/data1/liuyincheng/CMIP6-process/z_ensembleMean/MME1/CMIP/amip_1980-2014
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% load and read
    % data path
    dvarsPath = fullfile(MMEPath, level.process3{2}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/anomaly
    dradEffectPath = fullfile(MMEPath, level.process3{6}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/radEffect/
    load([dradEffectPath, 'global_vars.mat']) % 'lon_f', 'lat_f', 'timeEssmble', 'time', 'plev_k', 'readme', 'timeseries', 'MME_Models'
    ntime = length(time.date);

    % cal rhs(RHeating)
    load([dvarsPath, 'drlds.mat']) % surface_downwelling_longwave_flux_in_air
    load([dvarsPath, 'drsds.mat']) % surface_downwelling_shortwave_flux_in_air
    load([dvarsPath, 'drsus.mat']) % surface_upwelling_shortwave_flux_in_air
    dR_swnet = drsds - drsus; % sfc net shortwave flux
    drhs = drlds + dR_swnet; % equilibrium equation's RHS, nearly equal to sfc upward rad
    drhs = autoRegrid3(lon_k, lat_k, time.date, drhs, lon_f, lat_f, time.date);

    % cloud radRffect
    load([dradEffectPath, 'dR_cloud.mat'])

    % radEffect of models
    load([dradEffectPath, 'dradEffect_sfc_cld.mat']) % 'totalEffect', 'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', 'taEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect', 'mainEffect'
    drhsKern = albEffect + husEffect + taEffect + dR_cloud_sfc;

    time = time.date;
    ntime = length(time);
    startMonth = 1;

    % mask CHN US EUR
    areaStr = {'world', 'china east', 'USA east', 'EUR west'};
    latRange = 90;
    % world
    [drhsMask.world, ~, ~] = maskArea(drhs, lat_f, latRange, -latRange, 'world');
    [drhsKernMask.world, ~, ~] = maskArea(drhsKern, lat_f, latRange, -latRange, 'world');
    [dlwUPMask.world, ~, ~] = maskArea(-tsEffect, lat_f, latRange, -latRange, 'world');
    % china east
    [drhsMask.CHNeast, ~, ~] = maskArea(drhs, lat_f, latRange, -latRange, 'china east');
    [drhsKernMask.CHNeast, ~, ~] = maskArea(drhsKern, lat_f, latRange, -latRange, 'china east');
    [dlwUPMask.CHNeast, ~, ~] = maskArea(-tsEffect, lat_f, latRange, -latRange, 'china east');
    % US east
    [drhsMask.USeast, ~, ~] = maskArea(drhs, lat_f, latRange, -latRange, 'USA east');
    [drhsKernMask.USeast, ~, ~] = maskArea(drhsKern, lat_f, latRange, -latRange, 'USA east');
    [dlwUPMask.USeast, ~, ~] = maskArea(-tsEffect, lat_f, latRange, -latRange, 'USA east');
    % EUR west
    [drhsMask.EURwest, ~, ~] = maskArea(drhs, lat_f, latRange, -latRange, 'EUR west');
    [drhsKernMask.EURwest, ~, ~] = maskArea(drhsKern, lat_f, latRange, -latRange, 'EUR west');
    [dlwUPMask.EURwest, ~, ~] = maskArea(-tsEffect, lat_f, latRange, -latRange, 'EUR west');

    % save value
    outPutPath = fullfile(MMEPath, 'FigData/');
    auto_mkdir(outPutPath)
    save([outPutPath, 'regionalTsRHeating.mat'], 'drhsMask', 'drhsKernMask', 'dlwUPMask')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ERA5 part
[readme, level, tLin, vars] = obsParameters('ERA5');

for exmNum = 4:5 % 200001-201412 198001-201412
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %read data
    varsPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{1}); %rawdata
    dvarsPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{2}); %anomaly
    dvarsTrendPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{3}); %anomaly_trend
    kernelCalPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{4}); % kernelCal
    radEffectPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{5}); %radEffect
    dradTrendPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{6}); %/data1/liuyincheng/cmip6-proces/aimp_2000-2014/MRI-ESM2-0/ensemble/radEffect_trend/

    load([dvarsTrendPath, 'global_vars.mat']) % % 'lon_f', 'lat_f', 'lon_k', 'lat_k', 'plev_k', 'time'
    load([radEffectPath, 'dradEffect_sfc_cld.mat']) % load albEffect, husEffect, mainEffect, taEffect, taOnlyEffect, taOnlyEffect2, tasEffect, tasEffect2, totalEffect, tsEffect, wvlwEffect, wvswEffect
    load([radEffectPath, 'dR_cloud.mat'])
    load([dvarsPath, 'drhs.mat']) % load drhs
    nlonf = length(lon_f);
    nlatf = length(lat_f);
    drhsKern = albEffect + husEffect + taEffect + dR_cloud_sfc;

    % regrid
    drhs = autoRegrid3(lon_k, lat_k, time, drhs, lon_f, lat_f, time);

    % mask CHN US EUR
    areaStr = {'world', 'china east', 'USA east', 'EUR west'};
    latRange = 90;
    % world
    [drhsMask.world, ~, ~] = maskArea(drhs, lat_f, latRange, -latRange, 'world');
    [drhsKernMask.world, ~, ~] = maskArea(drhsKern, lat_f, latRange, -latRange, 'world');
    [dlwUPMask.world, ~, ~] = maskArea(-tsEffect, lat_f, latRange, -latRange, 'world');
    % china east
    [drhsMask.CHNeast, ~, ~] = maskArea(drhs, lat_f, latRange, -latRange, 'china east');
    [drhsKernMask.CHNeast, ~, ~] = maskArea(drhsKern, lat_f, latRange, -latRange, 'china east');
    [dlwUPMask.CHNeast, ~, ~] = maskArea(-tsEffect, lat_f, latRange, -latRange, 'china east');
    % US east
    [drhsMask.USeast, ~, ~] = maskArea(drhs, lat_f, latRange, -latRange, 'USA east');
    [drhsKernMask.USeast, ~, ~] = maskArea(drhsKern, lat_f, latRange, -latRange, 'USA east');
    [dlwUPMask.USeast, ~, ~] = maskArea(-tsEffect, lat_f, latRange, -latRange, 'USA east');
    % EUR west
    [drhsMask.EURwest, ~, ~] = maskArea(drhs, lat_f, latRange, -latRange, 'EUR west');
    [drhsKernMask.EURwest, ~, ~] = maskArea(drhsKern, lat_f, latRange, -latRange, 'EUR west');
    [dlwUPMask.EURwest, ~, ~] = maskArea(-tsEffect, lat_f, latRange, -latRange, 'EUR west');

    % save value
    outPutPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5','FigData/');
    auto_mkdir(outPutPath)
    save([outPutPath, 'regionalTsRHeating.mat'], 'drhsMask', 'drhsKernMask', 'dlwUPMask')
end

eval(['cd ', nowpath]);
t = toc; disp(t)
