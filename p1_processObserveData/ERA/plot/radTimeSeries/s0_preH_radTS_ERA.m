%%---------------------------------------------------------
% Author       : LYC
% Date         : 2021-05-04 17:27:46
% LastEditors  : Please set LastEditors
% Description  : 
% FilePath     : /code/p1_processObserveData/ERA/plot/radTimeSeries/s0_preH_radTS_ERA.m
%  
%%---------------------------------------------------------
run '/home/liuyc/lib/tools/matlab/myTools/autoScript/preLoadVar.m'

chooseData='ERA5';
[readme, level, tLin, vars] = obsParameters(chooseData);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% different time series, 1mean 2000-03 to 2018-02(18*12). 2 mean 200207-201706(15*12)
for exmNum = 1:5 % '200003-201802', '200207-201706', '200003-201402','200001-201412','198001-201412','200003-202011'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data path
    exmPath = fullfile(level.obsPath, tLin.time{exmNum}, level.dataName);
    varsPath = fullfile(exmPath, level.standVarPath{1}); %rawdata
    dvarsPath = fullfile(exmPath, level.standVarPath{2}); %anomaly
    dradEffectPath = fullfile(exmPath, level.standVarPath{5}); %radEffect
    
    % rad Effect
    load([dradEffectPath, 'global_vars.mat']) % dR_cloud_sfc, dR_cloud_toa
    load([dradEffectPath, 'dR_cloud.mat']) % dR_cloud_sfc, dR_cloud_toa
    load([dradEffectPath, 'dR_residual_cld_sfc.mat']) % dR_residual_cld_sfc, lw_residual_cld_sfc, sw_residual_cld_sfc
    load([dradEffectPath, 'dradEffect_sfc_cld.mat']) % albEffect, husEffect, taEffect, tsEffect
    ntime = length(time);

    % LH and SH
    load([dvarsPath, 'dsshf.mat']) % Surface Upward Sensible Heat Flux
    load([dvarsPath, 'dslhf.mat']) % Surface Upward Latent Heat Flux

    % model output radiation
    load([dvarsPath, 'dssr.mat']) % surface_net_downward_shortwave_flux
    load([dvarsPath, 'dstrd.mat']) % Surface thermal radiation downwards
    load([dvarsPath, 'dstr.mat']) % surface_net_upward_longwave_flux

    dR_swnet = dssr;
    drhs = dstrd + dR_swnet; % surface radiative heating 
    drhs = autoRegrid3(lon_k, lat_k, time, drhs, lon_f, lat_f, time);
    drhsKern = albEffect + husEffect + taEffect + dR_cloud_sfc;

    drlus = -(dstr-dstrd); % surface_upwelling_longwave_flux_in_air 符号以向上为正
    dlwUPRawData = autoRegrid3(lon_k, lat_k, time, drlus, lon_f, lat_f, time);
    dlwUPKern = -tsEffect;

    dhfss = autoRegrid3(lon_k, lat_k, time, dsshf, lon_f, lat_f, time);
    dhfls = autoRegrid3(lon_k, lat_k, time, dslhf, lon_f, lat_f, time);

    run regnMask.m

    % save value
    outPutPath = fullfile(exmPath, 'FigData/');
    auto_mkdir(outPutPath)
    save([outPutPath, 'regionalVarsRad_sfc_cld.mat'], 'dR_cloudMask','dR_residualMask','dR_albMask','dR_husMask','dR_taMask','dR_tsMask','dhfssMask','dhflsMask',...
    'drhsMask', 'drhsKernMask', 'dlwUPRawDataMask','dlwUPKernMask')
end