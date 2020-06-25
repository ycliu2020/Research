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
for p_1 = 1:4%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    %CAMS-CSM1-0 didn't have sfc clear sky radiation, delete it
    [readme, Experiment, level, tLin, mPlev, vars] = modelParameters(p_1);
    % inputPath
    inputPath = ['/data1/liuyincheng/cmip6-process/', level.time1{p_1}]; %/data1/liuyincheng/cmip6-process/amip_1980-2014/
    startmonth = 1;

    % model loop
    for level1 = 1:length(level.model2)
        % load radEffect
        radEffectPath = [inputPath, level.model2{level1}, '/', level.process3{6}]; %/data1/liuyincheng/cmip6-process/amip_1980-2014/MRI-ESM2-0/radEffect/
        load([radEffectPath, 'global_vars.mat'])% lat lon time plevf readme
        load([radEffectPath, 'dCRF.mat'])% dCRF
        load([radEffectPath, 'dradEfect_toa_cld.mat'])%'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', 'taEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect'
        load([radEffectPath, 'dR_mainEffect_toa.mat'])%dR_mainEffect_toa (ta+alb+q+clod effect)
        load([radEffectPath, 'dR_residual_cld_toa.mat'])%dR_residual_cld_toa(co2 effect)
        load([radEffectPath, 'dR_cloud_toa.mat'])%dR_residual_cld_toa(co2 effect)
        
        dCRF_toa=squeeze(dCRF(:,:,:,2)); 
        nlon = length(lon); nlat = length(lat);
        % cal the trend(10 vars)
        [trendm_dRtoa_CRF, trends_dRtoa_CRF, trendyr_dRtoa_CRF, ~, ~] = autoCalTrend(dCRF_toa, nlon, nlat, time, startmonth);
        [trendm_dRtoa_ta, trends_dRtoa_ta, trendyr_dRtoa_ta, ~, ~] = autoCalTrend(taEffect, nlon, nlat, time, startmonth);
        [trendm_dRtoa_taOnly2, trends_dRtoa_taOnly2, trendyr_dRtoa_taOnly2, ~, ~] = autoCalTrend(taOnlyEffect2, nlon, nlat, time, startmonth);
        [trendm_dRtoa_tas2, trends_dRtoa_tas2, trendyr_dRtoa_tas2, ~, ~] = autoCalTrend(tasEffect2, nlon, nlat, time, startmonth);
        [trendm_dRtoa_mainEffect, trends_dRtoa_mainEffect, trendyr_dRtoa_mainEffect, ~, ~] = autoCalTrend(dR_mainEffect_toa, nlon, nlat, time, startmonth);
        [trendm_dRtoa_residual, trends_dRtoa_residual, trendyr_dRtoa_residual, ~, ~] = autoCalTrend(dR_residual_cld_toa, nlon, nlat, time, startmonth);
        [trendm_dRtoa_cloud, trends_dRtoa_cloud, trendyr_dRtoa_cloud, ~, ~] = autoCalTrend(dR_cloud_toa, nlon, nlat, time, startmonth);
        [trendm_dRtoa_q, trends_dRtoa_q, trendyr_dRtoa_q, ~, ~] = autoCalTrend(husEffect, nlon, nlat, time, startmonth);
        [trendm_dRtoa_alb, trends_dRtoa_alb, trendyr_dRtoa_alb, ~, ~] = autoCalTrend(albEffect, nlon, nlat, time, startmonth);
        [trendm_dRtoa_ts, trends_dRtoa_ts, trendyr_dRtoa_ts, ~, ~] = autoCalTrend(tsEffect, nlon, nlat, time, startmonth);

        trend_radefectOutPath = [inputPath, level.model2{level1}, '/', level.process3{7}]; %/data1/liuyincheng/cmip6-proces/aimp_2000-2014/MRI-ESM2-0/radEffect_trend/
        auto_mkdir(trend_radefectOutPath)
        readme.aboutTas2 = 'tas2 mean consider 2 levels near sfc';
        save([trend_radefectOutPath, 'global_vars.mat'], 'lon', 'lat', 'time', 'plevf', 'readme', 'timeseries', 'modelname')
        save([trend_radefectOutPath, 'trend_dradEfect_toa_cld.mat'], ...
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
        disp([level.model2{level1}, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{p_1}, ' era is done!'])
end

t = toc; disp(t)
