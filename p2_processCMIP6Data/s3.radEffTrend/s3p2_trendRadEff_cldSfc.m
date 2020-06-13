%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-09 15:52:00
% LastEditTime : 2020-06-10 19:36:41
% LastEditors  : LYC
% Description  : 
% FilePath     : /Research/p2_processCMIP6Data/s3.radEffTrend/s4_trendRadEff_cldSfc.m
%  
%%---------------------------------------------------------
%need to be edit
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% mothly data
% cal cal 1.radEfect trend,2.year trend
% this script mainly include:
% PS: m mean month, s mean season, yr mean year
%
% hope to use an unite scrip to read and process all the varies we need
% Compare the pre rad and ts data from 2000-03 to 2018-02(18 years)
%
% time:2000.01-2014.12(interval:15*12);1980.01-2014.12(interval:35*12); 2015.01-2099.12(interval:85*12)
% initial time in hist(1740 total): 1,561 of 1740(2000.03);1,321 of 1740(1980.01)
% initial time in futrue(1032 total): 1 of 1032(2015.01);
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear; clc; tic;
nowpath = pwd;
% model settings
for p_1 = [2  4]%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    %CAMS-CSM1-0 didn't have sfc clear sky radiation, delete it
    [readme, Experiment, level, tLin, mPlev, vars] = modelParameters(p_1);
    % inputPath
    inputPath = ['/data1/liuyincheng/cmip6-process/', level.time1{p_1}]; %/data1/liuyincheng/cmip6-process/amip_1980-2014/
    startmonth = 1;

    % model loop
    for level1 = 1:length(level.model2)
        % load radEffect
        radEffectPath = [inputPath, level.model2{level1}, '/', level.process3{6}]; %/data1/liuyincheng/cmip6-process/amip_1980-2014/MRI-ESM2-0/radEffect/
        load([radEffectPath, 'dCRF.mat'])% dCRF
        load([radEffectPath, 'global_vars.mat'])% lat lon time plevf readme
        load([radEffectPath, 'dradEfect_sfc_cld.mat'])%'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', 'taEffect', 'tasEffect', taOnlyEffect', tasEffect2', 'taOnlyEffect2', 'totalEffect'
        load([radEffectPath, 'dR_tsAtom_cld.mat'])% dR_tsAtom_cld(ts effect on atoms)
        load([radEffectPath, 'dR_mainEffect_sfc.mat'])%dR_mainEffect_sfc (ta+alb+q+clod effect)
        load([radEffectPath, 'dR_residual_cld_sfc.mat'])%dR_residual_cld_sfc(co2 effect)
        load([radEffectPath, 'dR_cloud_sfc.mat'])%dR_residual_cld_sfc(co2 effect)
        
        dCRF_sfc=squeeze(dCRF(:,:,:,1)); 
        nlon = length(lon); nlat = length(lat);
        % cal the trend(10 vars)
        [trendm_dRsfc_CRF, trends_dRsfc_CRF, trendyr_dRsfc_CRF, ~, ~] = autoCalTrend(dCRF_sfc, nlon, nlat, time, startmonth);
        [trendm_dRsfc_ta, trends_dRsfc_ta, trendyr_dRsfc_ta, ~, ~] = autoCalTrend(taEffect, nlon, nlat, time, startmonth);
        [trendm_dRsfc_taOnly, trends_dRsfc_taOnly, trendyr_dRsfc_taOnly, ~, ~] = autoCalTrend(taOnlyEffect, nlon, nlat, time, startmonth);
        [trendm_dRsfc_taOnly2, trends_dRsfc_taOnly2, trendyr_dRsfc_taOnly2, ~, ~] = autoCalTrend(taOnlyEffect2, nlon, nlat, time, startmonth);
        [trendm_dRsfc_tas, trends_dRsfc_tas, trendyr_dRsfc_tas, ~, ~] = autoCalTrend(tasEffect, nlon, nlat, time, startmonth);
        [trendm_dRsfc_tas2, trends_dRsfc_tas2, trendyr_dRsfc_tas2, ~, ~] = autoCalTrend(tasEffect2, nlon, nlat, time, startmonth);
        [trendm_dRatm_tsAtom, trends_dRatm_tsAtom, trendyr_dRatm_tsAtom, ~, ~] = autoCalTrend(dR_tsAtom_cld, nlon, nlat, time, startmonth);
        [trendm_dRsfc_mainEffect, trends_dRsfc_mainEffect, trendyr_dRsfc_mainEffect, ~, ~] = autoCalTrend(dR_mainEffect_sfc, nlon, nlat, time, startmonth);
        [trendm_dRsfc_residual, trends_dRsfc_residual, trendyr_dRsfc_residual, ~, ~] = autoCalTrend(dR_residual_cld_sfc, nlon, nlat, time, startmonth);
        [trendm_dRsfc_cloud, trends_dRsfc_cloud, trendyr_dRsfc_cloud, ~, ~] = autoCalTrend(dR_cloud_sfc, nlon, nlat, time, startmonth);
        [trendm_dRsfc_q, trends_dRsfc_q, trendyr_dRsfc_q, ~, ~] = autoCalTrend(husEffect, nlon, nlat, time, startmonth);
        [trendm_dRsfc_alb, trends_dRsfc_alb, trendyr_dRsfc_alb, ~, ~] = autoCalTrend(albEffect, nlon, nlat, time, startmonth);
        [trendm_dRsfc_ts, trends_dRsfc_ts, trendyr_dRsfc_ts, ~, ~] = autoCalTrend(tsEffect, nlon, nlat, time, startmonth);

        trend_radefectOutPath = [inputPath, level.model2{level1}, '/', level.process3{7}]; %/data1/liuyincheng/cmip6-proces/aimp_2000-2014/MRI-ESM2-0/radEffect_trend/
        auto_mkdir(trend_radefectOutPath)
        readme.aboutTas2 = 'tas2 mean consider 2 levels near sfc';
        save([trend_radefectOutPath, 'global_vars.mat'], 'lon', 'lat', 'time', 'plevf', 'readme', 'timeseries', 'modelname')
        save([trend_radefectOutPath, 'trend_dradEfect_sfc_cld.mat'], ...
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
        disp([level.model2{level1}, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{p_1}, ' era is done!'])
end

t = toc; disp(t)
