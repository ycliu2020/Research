%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-07-06 13:27:23
% LastEditTime : 2021-04-08 17:21:40
% LastEditors  : Please set LastEditors
% Description  :
% FilePath     : /code/p1_processObserveData/ERA5/ERA5_s4_trendVars_dp.m
%
%%---------------------------------------------------------
[readme, level, tLin, vars] = obsParameters('ERA5');
for p_1 = 4:5
    varsPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{1}); %rawdata
    dvarsPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{2}); %anomaly
    kernelCalPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{4}); % kernelCal
    radEfectPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{5}); %radEffect
    trend_radEfectPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{6});; %/data1/liuyincheng/cmip6-proces/aimp_2000-2014/MRI-ESM2-0/ensemble/radEffect_trend/
    auto_mkdir(trend_radEfectPath)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part1: SFC vars
    %%%%%1.1 cal radEffect
    % load radEffect
    load([radEfectPath, 'dCRF.mat'])% dCRF
    load([radEfectPath, 'global_vars.mat'])% 'lon_f', 'lat_f', 'lon_k', 'lat_k', 'plev_k', 'time'
    load([radEfectPath, 'dradEffect_sfc_cld.mat'])%'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', 'taEffect', 'tasEffect', taOnlyEffect', tasEffect2', 'taOnlyEffect2', 'totalEffect'
    load([radEfectPath, 'dR_tsAtom_cld.mat'])% dR_tsAtom_cld(ts effect on atoms)
    load([radEfectPath, 'dR_mainEffect_sfc.mat'])%dR_mainEffect_sfc (ta+alb+q+clod effect)
    load([radEfectPath, 'dR_residual_cld_sfc.mat'])%dR_residual_cld_sfc(co2 effect)
    load([radEfectPath, 'dR_cloud_sfc.mat'])%dR_residual_cld_sfc(co2 effect)

    dCRF_sfc = squeeze(dCRF(:, :, :, 1));
    nlonf = length(lon_f); nlatf = length(lat_f);
    timeDate = time;
    startMonth = tLin.startMonth{p_1};

    % cal the trend(10 vars)
    [trendm_dRsfc_CRF, trends_dRsfc_CRF, trendyr_dRsfc_CRF, ~, ~] = autoCalTrend(dCRF_sfc, nlonf, nlatf, timeDate, startMonth);
    [trendm_dRsfc_ta, trends_dRsfc_ta, trendyr_dRsfc_ta, ~, ~] = autoCalTrend(taEffect, nlonf, nlatf, timeDate, startMonth);
    [trendm_dRsfc_taOnly, trends_dRsfc_taOnly, trendyr_dRsfc_taOnly, ~, ~] = autoCalTrend(taOnlyEffect, nlonf, nlatf, timeDate, startMonth);
    [trendm_dRsfc_taOnly2, trends_dRsfc_taOnly2, trendyr_dRsfc_taOnly2, ~, ~] = autoCalTrend(taOnlyEffect2, nlonf, nlatf, timeDate, startMonth);
    [trendm_dRsfc_tas, trends_dRsfc_tas, trendyr_dRsfc_tas, ~, ~] = autoCalTrend(tasEffect, nlonf, nlatf, timeDate, startMonth);
    [trendm_dRsfc_tas2, trends_dRsfc_tas2, trendyr_dRsfc_tas2, ~, ~] = autoCalTrend(tasEffect2, nlonf, nlatf, timeDate, startMonth);
    [trendm_dRatm_tsAtom, trends_dRatm_tsAtom, trendyr_dRatm_tsAtom, ~, ~] = autoCalTrend(dR_tsAtom_cld, nlonf, nlatf, timeDate, startMonth);
    [trendm_dRsfc_mainEffect, trends_dRsfc_mainEffect, trendyr_dRsfc_mainEffect, ~, ~] = autoCalTrend(dR_mainEffect_sfc, nlonf, nlatf, timeDate, startMonth);
    [trendm_dRsfc_residual, trends_dRsfc_residual, trendyr_dRsfc_residual, ~, ~] = autoCalTrend(dR_residual_cld_sfc, nlonf, nlatf, timeDate, startMonth);
    [trendm_dRsfc_cloud, trends_dRsfc_cloud, trendyr_dRsfc_cloud, ~, ~] = autoCalTrend(dR_cloud_sfc, nlonf, nlatf, timeDate, startMonth);
    [trendm_dRsfc_q, trends_dRsfc_q, trendyr_dRsfc_q, ~, ~] = autoCalTrend(husEffect, nlonf, nlatf, timeDate, startMonth);
    [trendm_dRsfc_alb, trends_dRsfc_alb, trendyr_dRsfc_alb, ~, ~] = autoCalTrend(albEffect, nlonf, nlatf, timeDate, startMonth);
    [trendm_dRsfc_ts, trends_dRsfc_ts, trendyr_dRsfc_ts, ~, ~] = autoCalTrend(tsEffect, nlonf, nlatf, timeDate, startMonth);

    readme.aboutTas2 = 'tas2 mean consider 2 levels near sfc';
    save([trend_radEfectPath, 'global_vars.mat'], 'lon_f', 'lat_f', 'time', 'plev_k')
    save([trend_radEfectPath, 'trend_dradEffect_sfc_cld.mat'], ...
        'trendm_dRsfc_CRF', 'trends_dRsfc_CRF', 'trendyr_dRsfc_CRF', ...
        'trendm_dRsfc_ta', 'trends_dRsfc_ta', 'trendyr_dRsfc_ta', ...
        'trendm_dRsfc_taOnly', 'trends_dRsfc_taOnly', 'trendyr_dRsfc_taOnly', ...
        'trendm_dRsfc_taOnly2', 'trends_dRsfc_taOnly2', 'trendyr_dRsfc_taOnly2', ...
        'trendm_dRsfc_tas', 'trends_dRsfc_tas', 'trendyr_dRsfc_tas', ...
        'trendm_dRsfc_tas2', 'trends_dRsfc_tas2', 'trendyr_dRsfc_tas2', ...
        'trendm_dRatm_tsAtom', 'trends_dRatm_tsAtom', 'trendyr_dRatm_tsAtom', ...
        'trendm_dRsfc_mainEffect', 'trends_dRsfc_mainEffect', 'trendyr_dRsfc_mainEffect', ...
        'trendm_dRsfc_residual', 'trends_dRsfc_residual', 'trendyr_dRsfc_residual', ...
        'trendm_dRsfc_cloud', 'trends_dRsfc_cloud', 'trendyr_dRsfc_cloud', ...
        'trendm_dRsfc_q', 'trends_dRsfc_q', 'trendyr_dRsfc_q', ...
        'trendm_dRsfc_alb', 'trends_dRsfc_alb', 'trendyr_dRsfc_alb', ...
        'trendm_dRsfc_ts', 'trends_dRsfc_ts', 'trendyr_dRsfc_ts')

    %%%%% 1.2: cal the data output net radEffect
    load([radEfectPath, 'global_vars.mat'])% 'lon_f', 'lat_f', 'lon_k', 'lat_k', 'plev_k', 'time'
    load([radEfectPath, 'real_dradEffect.mat'])%'l_rad', 's_rad', 'dR_allsky', 'dR_clr', 'readme_realradEfect'
    realRadEffect = dR_allsky(:, :, :, 1);
    nlonf = length(lon_f); nlatf = length(lat_f);
    % cal the trend(10 vars)
    [trendm_dRsfc_realRad, trends_dRsfc_realRad, trendyr_dRsfc_realRad, ~, ~] = autoCalTrend(realRadEffect, nlonf, nlatf, timeDate, startMonth);
    save([trend_radEfectPath, 'trend_drealRadEfect_sfc_cld.mat'], ...
        'trendm_dRsfc_realRad', 'trends_dRsfc_realRad', 'trendyr_dRsfc_realRad')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part1: TOA vars
    %%%%%1.1 cal radEffect
    load([radEfectPath, 'global_vars.mat'])% 'lon_f', 'lat_f', 'lon_k', 'lat_k', 'plev_k', 'time'
    load([radEfectPath, 'dCRF.mat'])% dCRF
    load([radEfectPath, 'dradEffect_toa_cld.mat'])%'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', 'taEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect'
    load([radEfectPath, 'dR_mainEffect_toa.mat'])%dR_mainEffect_toa (ta+alb+q+clod effect)
    load([radEfectPath, 'dR_residual_cld_toa.mat'])%dR_residual_cld_toa(co2 effect)
    load([radEfectPath, 'dR_cloud_toa.mat'])%dR_residual_cld_toa(co2 effect)
    dCRF_toa = squeeze(dCRF(:, :, :, 2));
    % cal the trend(10 vars)
    [trendm_dRtoa_CRF, trends_dRtoa_CRF, trendyr_dRtoa_CRF, ~, ~] = autoCalTrend(dCRF_toa, nlonf, nlatf, timeDate, startMonth);
    [trendm_dRtoa_ta, trends_dRtoa_ta, trendyr_dRtoa_ta, ~, ~] = autoCalTrend(taEffect, nlonf, nlatf, timeDate, startMonth);
    [trendm_dRtoa_taOnly2, trends_dRtoa_taOnly2, trendyr_dRtoa_taOnly2, ~, ~] = autoCalTrend(taOnlyEffect2, nlonf, nlatf, timeDate, startMonth);
    [trendm_dRtoa_tas2, trends_dRtoa_tas2, trendyr_dRtoa_tas2, ~, ~] = autoCalTrend(tasEffect2, nlonf, nlatf, timeDate, startMonth);
    [trendm_dRtoa_mainEffect, trends_dRtoa_mainEffect, trendyr_dRtoa_mainEffect, ~, ~] = autoCalTrend(dR_mainEffect_toa, nlonf, nlatf, timeDate, startMonth);
    [trendm_dRtoa_residual, trends_dRtoa_residual, trendyr_dRtoa_residual, ~, ~] = autoCalTrend(dR_residual_cld_toa, nlonf, nlatf, timeDate, startMonth);
    [trendm_dRtoa_cloud, trends_dRtoa_cloud, trendyr_dRtoa_cloud, ~, ~] = autoCalTrend(dR_cloud_toa, nlonf, nlatf, timeDate, startMonth);
    [trendm_dRtoa_q, trends_dRtoa_q, trendyr_dRtoa_q, ~, ~] = autoCalTrend(husEffect, nlonf, nlatf, timeDate, startMonth);
    [trendm_dRtoa_alb, trends_dRtoa_alb, trendyr_dRtoa_alb, ~, ~] = autoCalTrend(albEffect, nlonf, nlatf, timeDate, startMonth);
    [trendm_dRtoa_ts, trends_dRtoa_ts, trendyr_dRtoa_ts, ~, ~] = autoCalTrend(tsEffect, nlonf, nlatf, timeDate, startMonth);
    readme.aboutTas2 = 'tas2 mean consider 2 levels near sfc';
    save([trend_radEfectPath, 'global_vars.mat'], 'lon_f', 'lat_f', 'time', 'plev_k')
    save([trend_radEfectPath, 'trend_dradEffect_toa_cld.mat'], ...
        'trendm_dRtoa_CRF', 'trends_dRtoa_CRF', 'trendyr_dRtoa_CRF', ...
        'trendm_dRtoa_ta', 'trends_dRtoa_ta', 'trendyr_dRtoa_ta', ...
        'trendm_dRtoa_taOnly2', 'trends_dRtoa_taOnly2', 'trendyr_dRtoa_taOnly2', ...
        'trendm_dRtoa_tas2', 'trends_dRtoa_tas2', 'trendyr_dRtoa_tas2', ...
        'trendm_dRtoa_mainEffect', 'trends_dRtoa_mainEffect', 'trendyr_dRtoa_mainEffect', ...
        'trendm_dRtoa_residual', 'trends_dRtoa_residual', 'trendyr_dRtoa_residual', ...
        'trendm_dRtoa_cloud', 'trends_dRtoa_cloud', 'trendyr_dRtoa_cloud', ...
        'trendm_dRtoa_q', 'trends_dRtoa_q', 'trendyr_dRtoa_q', ...
        'trendm_dRtoa_alb', 'trends_dRtoa_alb', 'trendyr_dRtoa_alb', ...
        'trendm_dRtoa_ts', 'trends_dRtoa_ts', 'trendyr_dRtoa_ts')

    %%%%Part2: cal the data output net radEffect
    load([radEfectPath, 'global_vars.mat'])% 'lon_f', 'lat_f', 'lon_k', 'lat_k', 'plev_k', 'time'
    load([radEfectPath, 'real_dradEffect.mat'])%'l_rad', 's_rad', 'dR_allsky', 'dR_clr', 'readme_realradEfect'
    realRadEffect = dR_allsky(:, :, :, 2);
    % cal the trend(10 vars)
    [trendm_dRtoa_realRad, trends_dRtoa_realRad, trendyr_dRtoa_realRad, ~, ~] = autoCalTrend(realRadEffect, nlonf, nlatf, timeDate, startMonth);
    save([trend_radEfectPath, 'trend_drealRadEfect_toa_cld.mat'], ...
        'trendm_dRtoa_realRad', 'trends_dRtoa_realRad', 'trendyr_dRtoa_realRad')


end
