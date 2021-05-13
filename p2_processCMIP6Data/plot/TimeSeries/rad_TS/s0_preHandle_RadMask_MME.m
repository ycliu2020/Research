%%---------------------------------------------------------
% Author       : LYC
% Date         : 2021-04-12 21:00:50
% LastEditors  : Please set LastEditors
% Description  :
% FilePath     : /code/p2_processCMIP6Data/plot/TimeSeries/rad_TS/s0_preHandle_RadMask_MME.m
%
%%---------------------------------------------------------
clc; clear; tic;
nowpath = pwd;
run '/home/liuyc/lib/tools/matlab/myTools/autoScript/preLoadVar.m'

MMEType='MME3';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CMIP6 part
% experiment
for exmNum = [1 2 4] %1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    %CAMS-CSM1-0 didn't have sfc clear sky radiation, delete it
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % MMEPath
    MMEPath= fullfile(level.path_MME, MMEType, level.time1{exmNum});%/data1/liuyincheng/CMIP6-process/2000-2014/
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% load and read
    % data path
    dvarsPath = fullfile(MMEPath, level.process3{2}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/anomaly
    dradEffectPath = fullfile(MMEPath, level.process3{6}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/radEffect/
    load([dradEffectPath, 'global_vars.mat']) % 'lon_f', 'lat_f', 'timeEssmble', 'time', 'plev_k', 'readme', 'timeseries', 'MME_Models'
    ntime = length(time.date);

    % rad Effect
    load([dradEffectPath, 'dR_cloud.mat']) % dR_cloud_sfc, dR_cloud_toa
    load([dradEffectPath, 'dR_residual_cld_sfc.mat']) % dR_residual_cld_sfc, lw_residual_cld_sfc, sw_residual_cld_sfc
    load([dradEffectPath, 'dradEffect_sfc_cld.mat']) % albEffect, husEffect, taEffect, tsEffect
    % time = time.date;
    ntime = length(time.date);
    startMonth = 1;

    % cal rhs(RHeating)
    load([dvarsPath, 'drlds.mat']) % surface_downwelling_longwave_flux_in_air
    load([dvarsPath, 'drsds.mat']) % surface_downwelling_shortwave_flux_in_air
    load([dvarsPath, 'drsus.mat']) % surface_upwelling_shortwave_flux_in_air
    load([dvarsPath, 'drlus.mat']) % surface_upwelling_longwave_flux_in_air
    
    % LH and SH
    load([dvarsPath, 'dhfss.mat']) % Surface Upward Sensible Heat Flux
    load([dvarsPath, 'dhfls.mat']) % Surface Upward Latent Heat Flux

    dR_swnet = drsds - drsus; % sfc net shortwave flux
    drhs = drlds + dR_swnet; % equilibrium equation's RHS, nearly equal to sfc upward rad
    drhs = autoRegrid3(lon_k, lat_k, time.date, drhs, lon_f, lat_f, time.date);
    drhsKern = albEffect + husEffect + taEffect + dR_cloud_sfc;

    dlwUPRawData = autoRegrid3(lon_k, lat_k, time.date, drlus, lon_f, lat_f, time.date);
    dlwUPKern = -tsEffect;

    dhfss = autoRegrid3(lon_k, lat_k, time.date, dhfss, lon_f, lat_f, time.date);
    dhfls = autoRegrid3(lon_k, lat_k, time.date, dhfls, lon_f, lat_f, time.date);

    run regnMask.m

    % save value
    outPutPath = fullfile(MMEPath, 'FigData/');
    auto_mkdir(outPutPath)
    save([outPutPath, 'regionalVarsRad_sfc_cld.mat'], 'dR_cloudMask','dR_residualMask','dR_albMask','dR_husMask','dR_taMask','dR_tsMask','dhfssMask','dhflsMask',...
    'drhsMask', 'drhsKernMask', 'dlwUPRawDataMask','dlwUPKernMask')
end

eval(['cd ', nowpath]);
t = toc; disp(t)
