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
for p_1 = 1:5%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
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
        % load([radEffectPath, 'dradEfect_toa_cld.mat'])%'wvlwEffect', 'wvswEffect', 'tsEffect',  'albEffect', 'husEffect', 'taEffect', 'tasEffect2', 'taOnlyEffect2','realRadEffect'
        % load([radEffectPath, 'dR_mainEffect_toa.mat'])%dR_mainEffect_toa (ta+alb+q+clod effect)
        load([radEffectPath, 'real_dradEfect.mat'])%'l_rad', 's_rad', 'dR_allsky', 'dR_clr', 'readme_realradEfect'

        % load([radEffectPath, 'dR_residual_cld_toa.mat'])%dR_residual_cld_toa(co2 effect)
        % load([radEffectPath, 'dR_cloud_toa.mat'])%dR_residual_cld_toa(co2 effect)
        realRadEffect = dR_allsky(:, :, :, 2);
        nlon = length(lon); nlat = length(lat);
        % cal the trend(10 vars)
        [trendm_dRtoa_realRad, trends_dRtoa_realRad, trendyr_dRtoa_realRad, ~, ~] = autoCalTrend(realRadEffect, nlon, nlat, time, startmonth);

        trend_radefectOutPath = [inputPath, level.model2{level1}, '/', level.process3{7}]; %/data1/liuyincheng/cmip6-proces/aimp_2000-2014/MRI-ESM2-0/radEffect_trend/
        auto_mkdir(trend_radefectOutPath)
        save([trend_radefectOutPath, 'trend_drealRadEfect_toa_cld.mat'], ...
            'trendm_dRtoa_realRad', 'trends_dRtoa_realRad', 'trendyr_dRtoa_realRad')
        disp([level.model2{level1}, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{p_1}, ' era is done!'])
end

t = toc; disp(t)
