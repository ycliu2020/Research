%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-08-31 17:00:15
% LastEditTime : 2021-05-10 14:45:37
% LastEditors  : Please set LastEditors
% Description  : MME result of 时间序列图
% FilePath     : /code/p3_paperFigIntegrate/Fig3_timeSeries_tsVsRheatingRad/s0_preHandle_plotRegional_MMEERA.m
% note :
%%---------------------------------------------------------
clc; clear; tic;
nowpath = pwd;
% load mask map
run '/home/liuyc/lib/tools/matlab/myTools/autoScript/preLoadVar.m'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CMIP6 part
% experiment
% for exmNum = 4:4 %1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
%     %CAMS-CSM1-0 didn't have sfc clear sky radiation, delete it
%     [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
%     % MMEPath
%     MMEPath = ['/data1/liuyincheng/CMIP6-process/z_ensembleMean/MME1/', level.time1{exmNum}]; %/data1/liuyincheng/CMIP6-process/z_ensembleMean/MME1/CMIP/amip_1980-2014
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% load and read
%     % data path
%     dvarsPath = fullfile(MMEPath, level.process3{2}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/anomaly
%     dradEffectPath = fullfile(MMEPath, level.process3{6}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/radEffect/
%     load([dradEffectPath, 'global_vars.mat']) % 'lon_f', 'lat_f', 'timeEssmble', 'time', 'plev_k', 'readme', 'timeseries', 'MME_Models'
%     ntime = length(time.date);

%     % cal rhs(RHeating)
%     load([dvarsPath, 'drlds.mat']) % surface_downwelling_longwave_flux_in_air
%     load([dvarsPath, 'drsds.mat']) % surface_downwelling_shortwave_flux_in_air
%     load([dvarsPath, 'drsus.mat']) % surface_upwelling_shortwave_flux_in_air
%     load([dvarsPath, 'drlus.mat']) % surface_upwelling_longwave_flux_in_air
%     %LH and SH
%     load([dvarsPath, 'dhfss.mat']) % Surface Upward Sensible Heat Flux
%     load([dvarsPath, 'dhfls.mat']) % Surface Upward Latent Heat Flux

%     dR_swnet = drsds - drsus; % sfc net shortwave flux
%     drhs = drlds + dR_swnet; % equilibrium equation's RHS, nearly equal to sfc upward rad
%     drhs = autoRegrid3(lon_k, lat_k, time.date, drhs, lon_f, lat_f, time.date);
%     dlwUPRawData = autoRegrid3(lon_k, lat_k, time.date, drlus, lon_f, lat_f, time.date);
%     dhfss = autoRegrid3(lon_k, lat_k, time.date, dhfss, lon_f, lat_f, time.date);
%     dhfls = autoRegrid3(lon_k, lat_k, time.date, dhfls, lon_f, lat_f, time.date);

%     % cloud radRffect
%     load([dradEffectPath, 'dR_cloud.mat'])

%     % radEffect of models
%     load([dradEffectPath, 'dradEffect_sfc_cld.mat']) % 'totalEffect', 'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', 'taEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect', 'mainEffect'
%     drhsKern = albEffect + husEffect + taEffect + dR_cloud_sfc;

%     time = time.date;
%     ntime = length(time);
%     startMonth = 1;

%     % mask CHN US EUR
%     areaStr = {'world', 'china', 'china east', 'USA east', 'EUR west'};
%     latRange = 90;
%     % world
%     [drhsMask.world, ~, ~] = maskArea(drhs, lat_f, latRange, -latRange, 'world');
%     [drhsKernMask.world, ~, ~] = maskArea(drhsKern, lat_f, latRange, -latRange, 'world');
%     [dlwUPKernMask.world, ~, ~] = maskArea(-tsEffect, lat_f, latRange, -latRange, 'world');
%     [dlwUPRawDataMask.world, ~, ~] = maskArea(dlwUPRawData, lat_f, latRange, -latRange, 'world');
%     [dhfssMask.world, ~, ~] = maskArea(dhfss, lat_f, latRange, -latRange, 'world');
%     [dhflsMask.world, ~, ~] = maskArea(dhfls, lat_f, latRange, -latRange, 'world');

%     % china
%     [drhsMask.CHN, ~, ~] = maskArea(drhs, lat_f, latRange, -latRange, 'china');
%     [drhsKernMask.CHN, ~, ~] = maskArea(drhsKern, lat_f, latRange, -latRange, 'china');
%     [dlwUPKernMask.CHN, ~, ~] = maskArea(-tsEffect, lat_f, latRange, -latRange, 'china');
%     [dlwUPRawDataMask.CHN, ~, ~] = maskArea(dlwUPRawData, lat_f, latRange, -latRange, 'china');
%     [dhfssMask.CHN, ~, ~] = maskArea(dhfss, lat_f, latRange, -latRange, 'china');
%     [dhflsMask.CHN, ~, ~] = maskArea(dhfls, lat_f, latRange, -latRange, 'china');

%     % china east
%     [drhsMask.CHNeast, ~, ~] = maskArea(drhs, lat_f, latRange, -latRange, 'china east');
%     [drhsKernMask.CHNeast, ~, ~] = maskArea(drhsKern, lat_f, latRange, -latRange, 'china east');
%     [dlwUPKernMask.CHNeast, ~, ~] = maskArea(-tsEffect, lat_f, latRange, -latRange, 'china east');
%     [dlwUPRawDataMask.CHNeast, ~, ~] = maskArea(dlwUPRawData, lat_f, latRange, -latRange, 'china east');
%     [dhfssMask.CHNeast, ~, ~] = maskArea(dhfss, lat_f, latRange, -latRange, 'china east');
%     [dhflsMask.CHNeast, ~, ~] = maskArea(dhfls, lat_f, latRange, -latRange, 'china east');

%     % US east
%     [drhsMask.USeast, ~, ~] = maskArea(drhs, lat_f, latRange, -latRange, 'USA east');
%     [drhsKernMask.USeast, ~, ~] = maskArea(drhsKern, lat_f, latRange, -latRange, 'USA east');
%     [dlwUPKernMask.USeast, ~, ~] = maskArea(-tsEffect, lat_f, latRange, -latRange, 'USA east');
%     [dlwUPRawDataMask.USeast, ~, ~] = maskArea(dlwUPRawData, lat_f, latRange, -latRange, 'USA east');
%     [dhfssMask.USeast, ~, ~] = maskArea(dhfss, lat_f, latRange, -latRange, 'USA east');
%     [dhflsMask.USeast, ~, ~] = maskArea(dhfls, lat_f, latRange, -latRange, 'USA east');

%     % EUR west
%     [drhsMask.EURwest, ~, ~] = maskArea(drhs, lat_f, latRange, -latRange, 'EUR west');
%     [drhsKernMask.EURwest, ~, ~] = maskArea(drhsKern, lat_f, latRange, -latRange, 'EUR west');
%     [dlwUPKernMask.EURwest, ~, ~] = maskArea(-tsEffect, lat_f, latRange, -latRange, 'EUR west');
%     [dlwUPRawDataMask.EURwest, ~, ~] = maskArea(dlwUPRawData, lat_f, latRange, -latRange, 'EUR west');
%     [dhfssMask.EURwest, ~, ~] = maskArea(dhfss, lat_f, latRange, -latRange, 'EUR west');
%     [dhflsMask.EURwest, ~, ~] = maskArea(dhfls, lat_f, latRange, -latRange, 'EUR west');

%     % save value
%     outPutPath = fullfile(MMEPath, 'FigData/');
%     auto_mkdir(outPutPath)
%     save([outPutPath, 'regionalTsRHeating.mat'], 'drhsMask', 'drhsKernMask', 'dlwUPKernMask', 'dlwUPRawDataMask', 'dhfssMask', 'dhflsMask')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ERA5 part
dataName = 'ERA5';
[readme, level, tLin, vars] = obsParameters(dataName);
obsPath = '/data2/liuyincheng/Observe-process/';
for exmNum = 4:5 % 200001-201412 198001-201412
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %read data
    varsPath = fullfile(obsPath, tLin.time{exmNum}, dataName, level.standVarPath{1}); %rawdata
    dvarsPath = fullfile(obsPath, tLin.time{exmNum}, dataName, level.standVarPath{2}); %anomaly
    dvarsTrendPath = fullfile(obsPath, tLin.time{exmNum}, dataName, level.standVarPath{3}); %anomaly_trend
    kernelCalPath = fullfile(obsPath, tLin.time{exmNum}, dataName, level.standVarPath{4}); % kernelCal
    radEffectPath = fullfile(obsPath, tLin.time{exmNum}, dataName, level.standVarPath{5}); %radEffect
    dradTrendPath = fullfile(obsPath, tLin.time{exmNum}, dataName, level.standVarPath{6}); %/data1/liuyincheng/cmip6-proces/aimp_2000-2014/MRI-ESM2-0/ensemble/radEffect_trend/

    load([dvarsTrendPath, 'global_vars.mat']) % % 'lon_f', 'lat_f', 'lon_k', 'lat_k', 'plev_k', 'time'
    load([radEffectPath, 'dradEffect_sfc_cld.mat']) % load albEffect, husEffect, mainEffect, taEffect, taOnlyEffect, taOnlyEffect2, tasEffect, tasEffect2, totalEffect, tsEffect, wvlwEffect, wvswEffect
    load([radEffectPath, 'dR_cloud.mat'])
    load([dvarsPath, 'drhs.mat']) % load drhs
    load([dvarsPath, 'dstrd.mat']) % load Surface thermal radiation downwards
    load([dvarsPath, 'dstr.mat']) % load surface_net_upward_longwave_flux
    load([dvarsPath, 'dslhf.mat']) % load surface_upward_latent_heat_flux
    load([dvarsPath, 'dsshf.mat']) % load surface_upward_sensible_heat_flux
    nlonf = length(lon_f);
    nlatf = length(lat_f);
    drhsKern = albEffect + husEffect + taEffect + dR_cloud_sfc;
    dlwUPRawData = dstrd - dstr;
    
    % regrid
    drhs = autoRegrid3(lon_k, lat_k, time, drhs, lon_f, lat_f, time);
    dlwUPRawData = autoRegrid3(lon_k, lat_k, time, dlwUPRawData, lon_f, lat_f, time);
    dhfss = autoRegrid3(lon_k, lat_k, time, dsshf, lon_f, lat_f, time);
    dhfls = autoRegrid3(lon_k, lat_k, time, dslhf, lon_f, lat_f, time);

    % mask CHN US EUR
    areaStr = {'world', 'china', 'china east', 'USA east', 'EUR west'};
    latRange = 90;
    % world
    [drhsMask.world, ~, ~] = maskArea(drhs, lat_f, latRange, -latRange, 'world');
    [drhsKernMask.world, ~, ~] = maskArea(drhsKern, lat_f, latRange, -latRange, 'world');
    [dlwUPKernMask.world, ~, ~] = maskArea(-tsEffect, lat_f, latRange, -latRange, 'world');
    [dlwUPRawDataMask.world, ~, ~] = maskArea(dlwUPRawData, lat_f, latRange, -latRange, 'world');
    [dhfssMask.world, ~, ~] = maskArea(dhfss, lat_f, latRange, -latRange, 'world');
    [dhflsMask.world, ~, ~] = maskArea(dhfls, lat_f, latRange, -latRange, 'world');

    % china
    [drhsMask.CHN, ~, ~] = maskArea(drhs, lat_f, latRange, -latRange, 'china');
    [drhsKernMask.CHN, ~, ~] = maskArea(drhsKern, lat_f, latRange, -latRange, 'china');
    [dlwUPKernMask.CHN, ~, ~] = maskArea(-tsEffect, lat_f, latRange, -latRange, 'china');
    [dlwUPRawDataMask.CHN, ~, ~] = maskArea(dlwUPRawData, lat_f, latRange, -latRange, 'china');
    [dhfssMask.CHN, ~, ~] = maskArea(dhfss, lat_f, latRange, -latRange, 'china');
    [dhflsMask.CHN, ~, ~] = maskArea(dhfls, lat_f, latRange, -latRange, 'china');

    % china east
    [drhsMask.CHNeast, ~, ~] = maskArea(drhs, lat_f, latRange, -latRange, 'china east');
    [drhsKernMask.CHNeast, ~, ~] = maskArea(drhsKern, lat_f, latRange, -latRange, 'china east');
    [dlwUPKernMask.CHNeast, ~, ~] = maskArea(-tsEffect, lat_f, latRange, -latRange, 'china east');
    [dlwUPRawDataMask.CHNeast, ~, ~] = maskArea(dlwUPRawData, lat_f, latRange, -latRange, 'china east');
    [dhfssMask.CHNeast, ~, ~] = maskArea(dhfss, lat_f, latRange, -latRange, 'china east');
    [dhflsMask.CHNeast, ~, ~] = maskArea(dhfls, lat_f, latRange, -latRange, 'china east');

    % US east
    [drhsMask.USeast, ~, ~] = maskArea(drhs, lat_f, latRange, -latRange, 'USA east');
    [drhsKernMask.USeast, ~, ~] = maskArea(drhsKern, lat_f, latRange, -latRange, 'USA east');
    [dlwUPKernMask.USeast, ~, ~] = maskArea(-tsEffect, lat_f, latRange, -latRange, 'USA east');
    [dlwUPRawDataMask.USeast, ~, ~] = maskArea(dlwUPRawData, lat_f, latRange, -latRange, 'USA east');
    [dhfssMask.USeast, ~, ~] = maskArea(dhfss, lat_f, latRange, -latRange, 'USA east');
    [dhflsMask.USeast, ~, ~] = maskArea(dhfls, lat_f, latRange, -latRange, 'USA east');

    % EUR west
    [drhsMask.EURwest, ~, ~] = maskArea(drhs, lat_f, latRange, -latRange, 'EUR west');
    [drhsKernMask.EURwest, ~, ~] = maskArea(drhsKern, lat_f, latRange, -latRange, 'EUR west');
    [dlwUPKernMask.EURwest, ~, ~] = maskArea(-tsEffect, lat_f, latRange, -latRange, 'EUR west');
    [dlwUPRawDataMask.EURwest, ~, ~] = maskArea(dlwUPRawData, lat_f, latRange, -latRange, 'EUR west');
    [dhfssMask.EURwest, ~, ~] = maskArea(dhfss, lat_f, latRange, -latRange, 'EUR west');
    [dhflsMask.EURwest, ~, ~] = maskArea(dhfls, lat_f, latRange, -latRange, 'EUR west');

    % save value
    outPutPath = fullfile(obsPath, tLin.time{exmNum}, dataName, 'FigData/');
    auto_mkdir(outPutPath)
    save([outPutPath, 'regionalTsRHeating.mat'], 'drhsMask', 'drhsKernMask', 'dlwUPKernMask', 'dlwUPRawDataMask', 'dhfssMask', 'dhflsMask')
end

eval(['cd ', nowpath]);
t = toc; disp(t)
