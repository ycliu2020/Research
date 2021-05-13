%%---------------------------------------------------------
% Author       : LYC
% Date         : 2021-04-12 21:00:50
% LastEditors  : Please set LastEditors
% Description  :
% FilePath     : /code/p3_paperFigIntegrate/Fig6_7_timeSeriANS_tsVsVarsRad/s0_preHandle_TsVsVarsRad.m
%
%%---------------------------------------------------------
clc; clear; tic;
nowpath = pwd;
run '/home/liuyc/lib/tools/matlab/myTools/autoScript/preLoadVar.m'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CMIP6 part
% the number of MME experiment
MMENum = 'MME1';
% experiment
for exmNum = 4:4 %1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
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
    load([dradEffectPath, 'dR_cloud.mat']) % dR_cloud_sfc, dR_cloud_toa
    load([dradEffectPath, 'dR_residual_cld_sfc.mat']) % dR_residual_cld_sfc, lw_residual_cld_sfc, sw_residual_cld_sfc
    load([dradEffectPath, 'dradEffect_sfc_cld.mat']) % albEffect, husEffect, taEffect, tsEffect

    time = time.date;
    ntime = length(time);
    startMonth = 1;

    % mask CHN US EUR
    areaStr = {'world', 'china east', 'USA east', 'EUR west'};
    latRange = 90;
    % world
    [dR_cloudMask.world, ~, ~] = maskArea(dR_cloud_sfc, lat_f, latRange, -latRange, 'world');
    [dR_residualMask.world, ~, ~] = maskArea(dR_residual_cld_sfc, lat_f, latRange, -latRange, 'world');
    [dR_albMask.world, ~, ~] = maskArea(albEffect, lat_f, latRange, -latRange, 'world');
    [dR_husMask.world, ~, ~] = maskArea(husEffect, lat_f, latRange, -latRange, 'world');
    [dR_taMask.world, ~, ~] = maskArea(taEffect, lat_f, latRange, -latRange, 'world');
    [dR_tsMask.world, ~, ~] = maskArea(tsEffect, lat_f, latRange, -latRange, 'world');

    % china east
    [dR_cloudMask.CHNeast, ~, ~] = maskArea(dR_cloud_sfc, lat_f, latRange, -latRange, 'china east');
    [dR_residualMask.CHNeast, ~, ~] = maskArea(dR_residual_cld_sfc, lat_f, latRange, -latRange, 'china east');
    [dR_albMask.CHNeast, ~, ~] = maskArea(albEffect, lat_f, latRange, -latRange, 'china east');
    [dR_husMask.CHNeast, ~, ~] = maskArea(husEffect, lat_f, latRange, -latRange, 'china east');
    [dR_taMask.CHNeast, ~, ~] = maskArea(taEffect, lat_f, latRange, -latRange, 'china east');
    [dR_tsMask.CHNeast, ~, ~] = maskArea(tsEffect, lat_f, latRange, -latRange, 'china east');

    % US east
    [dR_cloudMask.USeast, ~, ~] = maskArea(dR_cloud_sfc, lat_f, latRange, -latRange, 'USA east');
    [dR_residualMask.USeast, ~, ~] = maskArea(dR_residual_cld_sfc, lat_f, latRange, -latRange, 'USA east');
    [dR_albMask.USeast, ~, ~] = maskArea(albEffect, lat_f, latRange, -latRange, 'USA east');
    [dR_husMask.USeast, ~, ~] = maskArea(husEffect, lat_f, latRange, -latRange, 'USA east');
    [dR_taMask.USeast, ~, ~] = maskArea(taEffect, lat_f, latRange, -latRange, 'USA east');
    [dR_tsMask.USeast, ~, ~] = maskArea(tsEffect, lat_f, latRange, -latRange, 'USA east');

    % EUR west
    [dR_cloudMask.EURwest, ~, ~] = maskArea(dR_cloud_sfc, lat_f, latRange, -latRange, 'EUR west');
    [dR_residualMask.EURwest, ~, ~] = maskArea(dR_residual_cld_sfc, lat_f, latRange, -latRange, 'EUR west');
    [dR_albMask.EURwest, ~, ~] = maskArea(albEffect, lat_f, latRange, -latRange, 'EUR west');
    [dR_husMask.EURwest, ~, ~] = maskArea(husEffect, lat_f, latRange, -latRange, 'EUR west');
    [dR_taMask.EURwest, ~, ~] = maskArea(taEffect, lat_f, latRange, -latRange, 'EUR west');
    [dR_tsMask.EURwest, ~, ~] = maskArea(tsEffect, lat_f, latRange, -latRange, 'EUR west');

    % save value
    outPutPath = fullfile(MMEPath, 'FigData/');
    auto_mkdir(outPutPath)
    save([outPutPath, 'regionalVarsRad_sfc_cld.mat'], 'dR_cloudMask','dR_residualMask','dR_albMask','dR_husMask','dR_taMask','dR_tsMask')
end

eval(['cd ', nowpath]);
t = toc; disp(t)