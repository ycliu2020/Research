%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-09 15:52:00
% LastEditTime : 2020-11-08 15:58:48
% LastEditors  : LYC
% Description  : cal mainly include 1.regrid vars, 2.vars anomly
%                CMIP6 mothly data
%                time:2000.01-2014.12(interval:15*12);1980.01-2014.12(interval:35*12); 2015.01-2099.12(interval:85*12)
%                initial time in amip(432 total): 253 of 432(2000.03);13 of 432(1980.01);
%                initial time in futrue(1032 total): 1 of 1032(2015.01);
%                initial time in amip-hist(1740 total): 1,561 of 1740(2000.03);1,321 of 1740(1980.01);
%                PS: m mean month, s mean season, yr mean year
% FilePath     : /code/p2_processCMIP6Data/s2.radEffTrend/s3_trendRadEff_parfor.m
%
%%---------------------------------------------------------
clear; clc; tic;
nowpath = pwd;

% 开启并行环境
poolobj = gcp('nocreate'); % If no pool,  create new one.

if isempty(poolobj)
    MyPar = parpool(16);
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

    %% parameters and path
    lon_k = 0:2.5:357.5; nlonk = length(lon_k);
    lat_k = 90:-2.5:-90; nlatk = length(lat_k);
    lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
    lon_f = lon_k; nlonf = length(lon_f);
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    startmonth = 1;
    % ensemble member path
    esmName = getPath_fileName(mdlPath, '.');
    disp([level.model2{mdlNum}, ', ', level.time1{exmNum}, ', ', esmName{esmNum, 1}, ' ensemble start!'])

    eval(['cd ', esmPath]);
    % input and outputpath
    radEffectPath = fullfile(esmPath, level.process3{6}); %/data1/liuyincheng/cmip6-process/amip_1980-2014/MRI-ESM2-0/ensemble/radEffect/
    trend_radEffectPath = fullfile(esmPath, level.process3{7}); %/data1/liuyincheng/cmip6-proces/aimp_2000-2014/MRI-ESM2-0/ensemble/radEffect_trend/
    auto_mkdir(trend_radEffectPath)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % cloud SFC rad trend
    %%%%Part1: cal radEffect
    % load radEffect
    load([radEffectPath, 'dCRF.mat'], 'dCRF')% dCRF
    load([radEffectPath, 'global_vars.mat'])% lat_f lon_f time plev_k readme
    load([radEffectPath, 'dradEffect_sfc_cld.mat'])%'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', 'taEffect', 'tasEffect', taOnlyEffect', tasEffect2', 'taOnlyEffect2', 'totalEffect'
    load([radEffectPath, 'dR_tsAtom_cld.mat'])% dR_tsAtom_cld(ts effect on atoms)
    load([radEffectPath, 'dR_nonTs_sfc.mat'])%dR_nonTs_sfc (ta+alb+q+clod effect)
    load([radEffectPath, 'dR_residual_cld_sfc.mat'])%dR_residual_cld_sfc(co2 effect)
    load([radEffectPath, 'dR_cloud.mat'])%dR_residual_cld_sfc(co2 effect)
    clear -regexp ^lw ^sw
    dCRF_sfc = squeeze(dCRF(:, :, :, 1));
    nlonf = length(lon_f); nlatf = length(lat_f);
    timeDate = time.date;
    % cal the trend(10 vars)
    [trendm_dRsfc_CRF, trends_dRsfc_CRF, trendyr_dRsfc_CRF, ~, ~] = autoCalTrend(dCRF_sfc, nlonf, nlatf, timeDate, startmonth);
    [trendm_dRsfc_ta, trends_dRsfc_ta, trendyr_dRsfc_ta, ~, ~] = autoCalTrend(taEffect, nlonf, nlatf, timeDate, startmonth);
    [trendm_dRsfc_taOnly, trends_dRsfc_taOnly, trendyr_dRsfc_taOnly, ~, ~] = autoCalTrend(taOnlyEffect, nlonf, nlatf, timeDate, startmonth);
    [trendm_dRsfc_taOnly2, trends_dRsfc_taOnly2, trendyr_dRsfc_taOnly2, ~, ~] = autoCalTrend(taOnlyEffect2, nlonf, nlatf, timeDate, startmonth);
    [trendm_dRsfc_tas, trends_dRsfc_tas, trendyr_dRsfc_tas, ~, ~] = autoCalTrend(tasEffect, nlonf, nlatf, timeDate, startmonth);
    [trendm_dRsfc_tas2, trends_dRsfc_tas2, trendyr_dRsfc_tas2, ~, ~] = autoCalTrend(tasEffect2, nlonf, nlatf, timeDate, startmonth);
    [trendm_dRatm_tsAtom, trends_dRatm_tsAtom, trendyr_dRatm_tsAtom, ~, ~] = autoCalTrend(dR_tsAtom_cld, nlonf, nlatf, timeDate, startmonth);
    [trendm_dRsfc_nonTs, trends_dRsfc_nonTs, trendyr_dRsfc_nonTs, ~, ~] = autoCalTrend(dR_nonTs_sfc, nlonf, nlatf, timeDate, startmonth);
    [trendm_dRsfc_residual, trends_dRsfc_residual, trendyr_dRsfc_residual, ~, ~] = autoCalTrend(dR_residual_cld_sfc, nlonf, nlatf, timeDate, startmonth);
    [trendm_dRsfc_cloud, trends_dRsfc_cloud, trendyr_dRsfc_cloud, ~, ~] = autoCalTrend(dR_cloud_sfc, nlonf, nlatf, timeDate, startmonth);
    [trendm_dRsfc_q, trends_dRsfc_q, trendyr_dRsfc_q, ~, ~] = autoCalTrend(husEffect, nlonf, nlatf, timeDate, startmonth);
    [trendm_dRsfc_alb, trends_dRsfc_alb, trendyr_dRsfc_alb, ~, ~] = autoCalTrend(albEffect, nlonf, nlatf, timeDate, startmonth);
    [trendm_dRsfc_ts, trends_dRsfc_ts, trendyr_dRsfc_ts, ~, ~] = autoCalTrend(tsEffect, nlonf, nlatf, timeDate, startmonth);
    clear dCRF_sfc 
    clear -regexp ^alb ^ts ^ta ^hus
    readme.aboutTas2 = 'tas2 mean consider 2 levels near sfc';
    save([trend_radEffectPath, 'global_vars.mat'], 'lon_f', 'lat_f', 'time', 'plev_k', 'readme', 'timeseries', 'modelname')
    save([trend_radEffectPath, 'trend_dradEffect_sfc_cld.mat'], ...
        'trendm_dRsfc_CRF', 'trends_dRsfc_CRF', 'trendyr_dRsfc_CRF', ...
        'trendm_dRsfc_ta', 'trends_dRsfc_ta', 'trendyr_dRsfc_ta', ...
        'trendm_dRsfc_taOnly', 'trends_dRsfc_taOnly', 'trendyr_dRsfc_taOnly', ...
        'trendm_dRsfc_taOnly2', 'trends_dRsfc_taOnly2', 'trendyr_dRsfc_taOnly2', ...
        'trendm_dRsfc_tas', 'trends_dRsfc_tas', 'trendyr_dRsfc_tas', ...
        'trendm_dRsfc_tas2', 'trends_dRsfc_tas2', 'trendyr_dRsfc_tas2', ...
        'trendm_dRatm_tsAtom', 'trends_dRatm_tsAtom', 'trendyr_dRatm_tsAtom', ...
        'trendm_dRsfc_nonTs', 'trends_dRsfc_nonTs', 'trendyr_dRsfc_nonTs', ...
        'trendm_dRsfc_residual', 'trends_dRsfc_residual', 'trendyr_dRsfc_residual', ...
        'trendm_dRsfc_cloud', 'trends_dRsfc_cloud', 'trendyr_dRsfc_cloud', ...
        'trendm_dRsfc_q', 'trends_dRsfc_q', 'trendyr_dRsfc_q', ...
        'trendm_dRsfc_alb', 'trends_dRsfc_alb', 'trendyr_dRsfc_alb', ...
        'trendm_dRsfc_ts', 'trends_dRsfc_ts', 'trendyr_dRsfc_ts')
    clear -regexp ^trendm ^trends ^trendyr
    %%%%Part2: cal the model output net radEffect
    load([radEffectPath, 'global_vars.mat'])% lat_f lon_f time plev_k readme
    load([radEffectPath, 'model_dradEffect.mat'])%'lw_net', 'sw_net', 'dR_net_cld', 'dR_net_clr', 'readme_realradEfect'
    modelRadEffect = dR_net_cld(:, :, :, 1);
    nlonf = length(lon_f); nlatf = length(lat_f);
    % cal the trend(10 vars)
    [trendm_dRsfc_modelRad, trends_dRsfc_modelRad, trendyr_dRsfc_modelRad, ~, ~] = autoCalTrend(modelRadEffect, nlonf, nlatf, timeDate, startmonth);
    save([trend_radEffectPath, 'trend_dmodelRadEfect_sfc_cld.mat'], ...
        'trendm_dRsfc_modelRad', 'trends_dRsfc_modelRad', 'trendyr_dRsfc_modelRad')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % cloud TOA rad trend
    %%%%Part1: cal radEffect
    load([radEffectPath, 'global_vars.mat'])% lat_f lon_f time plev_k readme
    load([radEffectPath, 'dCRF.mat'],'dCRF')% dCRF
    load([radEffectPath, 'dradEffect_toa_cld.mat'])%'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', 'taEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect'
    load([radEffectPath, 'dR_nonTs_toa.mat'])%dR_nonTs_toa (ta+alb+q+clod effect)
    load([radEffectPath, 'dR_residual_cld_toa.mat'])%dR_residual_cld_toa(co2 effect)
    load([radEffectPath, 'dR_cloud.mat'])%dR_residual_cld_toa(co2 effect)
    clear -regexp ^lw ^sw
    dCRF_toa = squeeze(dCRF(:, :, :, 2));
    nlonf = length(lon_f); nlatf = length(lat_f);
    timeDate = time.date;
    % cal the trend(10 vars)
    [trendm_dRtoa_CRF, trends_dRtoa_CRF, trendyr_dRtoa_CRF, ~, ~] = autoCalTrend(dCRF_toa, nlonf, nlatf, timeDate, startmonth);
    [trendm_dRtoa_ta, trends_dRtoa_ta, trendyr_dRtoa_ta, ~, ~] = autoCalTrend(taEffect, nlonf, nlatf, timeDate, startmonth);
    [trendm_dRtoa_taOnly2, trends_dRtoa_taOnly2, trendyr_dRtoa_taOnly2, ~, ~] = autoCalTrend(taOnlyEffect2, nlonf, nlatf, timeDate, startmonth);
    [trendm_dRtoa_tas2, trends_dRtoa_tas2, trendyr_dRtoa_tas2, ~, ~] = autoCalTrend(tasEffect2, nlonf, nlatf, timeDate, startmonth);
    [trendm_dRtoa_nonTs, trends_dRtoa_nonTs, trendyr_dRtoa_nonTs, ~, ~] = autoCalTrend(dR_nonTs_toa, nlonf, nlatf, timeDate, startmonth);
    [trendm_dRtoa_residual, trends_dRtoa_residual, trendyr_dRtoa_residual, ~, ~] = autoCalTrend(dR_residual_cld_toa, nlonf, nlatf, timeDate, startmonth);
    [trendm_dRtoa_cloud, trends_dRtoa_cloud, trendyr_dRtoa_cloud, ~, ~] = autoCalTrend(dR_cloud_toa, nlonf, nlatf, timeDate, startmonth);
    [trendm_dRtoa_q, trends_dRtoa_q, trendyr_dRtoa_q, ~, ~] = autoCalTrend(husEffect, nlonf, nlatf, timeDate, startmonth);
    [trendm_dRtoa_alb, trends_dRtoa_alb, trendyr_dRtoa_alb, ~, ~] = autoCalTrend(albEffect, nlonf, nlatf, timeDate, startmonth);
    [trendm_dRtoa_ts, trends_dRtoa_ts, trendyr_dRtoa_ts, ~, ~] = autoCalTrend(tsEffect, nlonf, nlatf, timeDate, startmonth);
    clear dCRF_toa 
    clear -regexp ^alb ^ts ^ta ^hus

    readme.aboutTas2 = 'tas2 mean consider 2 levels near sfc';
    save([trend_radEffectPath, 'global_vars.mat'], 'lon_f', 'lat_f', 'time', 'plev_k', 'readme', 'timeseries', 'modelname')
    save([trend_radEffectPath, 'trend_dradEffect_toa_cld.mat'], ...
        'trendm_dRtoa_CRF', 'trends_dRtoa_CRF', 'trendyr_dRtoa_CRF', ...
        'trendm_dRtoa_ta', 'trends_dRtoa_ta', 'trendyr_dRtoa_ta', ...
        'trendm_dRtoa_taOnly2', 'trends_dRtoa_taOnly2', 'trendyr_dRtoa_taOnly2', ...
        'trendm_dRtoa_tas2', 'trends_dRtoa_tas2', 'trendyr_dRtoa_tas2', ...
        'trendm_dRtoa_nonTs', 'trends_dRtoa_nonTs', 'trendyr_dRtoa_nonTs', ...
        'trendm_dRtoa_residual', 'trends_dRtoa_residual', 'trendyr_dRtoa_residual', ...
        'trendm_dRtoa_cloud', 'trends_dRtoa_cloud', 'trendyr_dRtoa_cloud', ...
        'trendm_dRtoa_q', 'trends_dRtoa_q', 'trendyr_dRtoa_q', ...
        'trendm_dRtoa_alb', 'trends_dRtoa_alb', 'trendyr_dRtoa_alb', ...
        'trendm_dRtoa_ts', 'trends_dRtoa_ts', 'trendyr_dRtoa_ts')
    clear -regexp ^trendm ^trends ^trendyr
    %%%%Part2: cal the model output net radEffect
    load([radEffectPath, 'global_vars.mat'])% lat_f lon_f time plev_k readme
    load([radEffectPath, 'model_dradEffect.mat'])%'lw_net', 'sw_net', 'dR_net_cld', 'dR_net_clr', 'readme_realradEfect'
    modelRadEffect = dR_net_cld(:, :, :, 2);
    nlonf = length(lon_f); nlatf = length(lat_f);
    % cal the trend(10 vars)
    [trendm_dRtoa_modelRad, trends_dRtoa_modelRad, trendyr_dRtoa_modelRad, ~, ~] = autoCalTrend(modelRadEffect, nlonf, nlatf, timeDate, startmonth);
    save([trend_radEffectPath, 'trend_dmodelRadEfect_toa_cld.mat'], ...
        'trendm_dRtoa_modelRad', 'trends_dRtoa_modelRad', 'trendyr_dRtoa_modelRad')

    disp([level.model2{mdlNum}, ', ', level.time1{exmNum}, ', ', esmName{esmNum, 1}, ' ensemble is done!'])
    disp(' ')

end
