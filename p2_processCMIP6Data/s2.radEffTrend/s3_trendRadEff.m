%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-09 15:52:00
% LastEditTime : 2020-07-11 10:20:21
% LastEditors  : LYC
% Description  : cal mainly include 1.regrid vars, 2.vars anomly
%                CMIP6 mothly data
%                time:2000.01-2014.12(interval:15*12);1980.01-2014.12(interval:35*12); 2015.01-2099.12(interval:85*12)
%                initial time in amip(432 total): 253 of 432(2000.03);13 of 432(1980.01);
%                initial time in futrue(1032 total): 1 of 1032(2015.01);
%                initial time in amip-hist(1740 total): 1,561 of 1740(2000.03);1,321 of 1740(1980.01);
%                PS: m mean month, s mean season, yr mean year
% FilePath     : /code/p2_processCMIP6Data/s2.radEffTrend/s3_trendRadEff.m
%
%%---------------------------------------------------------
clear; clc; tic;
nowpath = pwd;
lon_k = 0:2.5:357.5; nlonk = length(lon_k);
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment
for p_1 = 4:4%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    %CAMS-CSM1-0 didn't have sfc clear sky radiation, delete it
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(p_1);
    % exmPath
    exmPath = ['/data1/liuyincheng/cmip6-process/', level.time1{p_1}]; %/data1/liuyincheng/cmip6-process/amip_1980-2014/
    startmonth = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model
    for level1 = 1:length(level.model2)
        % model path
        mdlPath = fullfile(exmPath, level.model2{level1});
        eval(['cd ', mdlPath]);
        disp(' ')
        disp([level.model2{level1}, ' model start!'])
        % ensemble member path
        esmName = getPath_fileName(mdlPath, '.');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensemble member
        for esmNum = 1:length(esmName)
            esmPath = fullfile(mdlPath, esmName{esmNum, 1});
            eval(['cd ', esmPath]);
            % input and outputpath
            radEfectPath = fullfile(esmPath, level.process3{6}); %/data1/liuyincheng/cmip6-process/amip_1980-2014/MRI-ESM2-0/ensemble/radEffect/
            trend_radEfectPath = fullfile(esmPath, level.process3{7}); %/data1/liuyincheng/cmip6-proces/aimp_2000-2014/MRI-ESM2-0/ensemble/radEffect_trend/
            auto_mkdir(trend_radEfectPath)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % cloud SFC rad trend
            %%%%Part1: cal radEffect
            % load radEffect
            load([radEfectPath, 'dCRF.mat'])% dCRF
            load([radEfectPath, 'global_vars.mat'])% lat_f lon_f time plev_k readme
            load([radEfectPath, 'dradEfect_sfc_cld.mat'])%'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', 'taEffect', 'tasEffect', taOnlyEffect', tasEffect2', 'taOnlyEffect2', 'totalEffect'
            load([radEfectPath, 'dR_tsAtom_cld.mat'])% dR_tsAtom_cld(ts effect on atoms)
            load([radEfectPath, 'dR_mainEffect_sfc.mat'])%dR_mainEffect_sfc (ta+alb+q+clod effect)
            load([radEfectPath, 'dR_residual_cld_sfc.mat'])%dR_residual_cld_sfc(co2 effect)
            load([radEfectPath, 'dR_cloud_sfc.mat'])%dR_residual_cld_sfc(co2 effect)

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
            [trendm_dRsfc_mainEffect, trends_dRsfc_mainEffect, trendyr_dRsfc_mainEffect, ~, ~] = autoCalTrend(dR_mainEffect_sfc, nlonf, nlatf, timeDate, startmonth);
            [trendm_dRsfc_residual, trends_dRsfc_residual, trendyr_dRsfc_residual, ~, ~] = autoCalTrend(dR_residual_cld_sfc, nlonf, nlatf, timeDate, startmonth);
            [trendm_dRsfc_cloud, trends_dRsfc_cloud, trendyr_dRsfc_cloud, ~, ~] = autoCalTrend(dR_cloud_sfc, nlonf, nlatf, timeDate, startmonth);
            [trendm_dRsfc_q, trends_dRsfc_q, trendyr_dRsfc_q, ~, ~] = autoCalTrend(husEffect, nlonf, nlatf, timeDate, startmonth);
            [trendm_dRsfc_alb, trends_dRsfc_alb, trendyr_dRsfc_alb, ~, ~] = autoCalTrend(albEffect, nlonf, nlatf, timeDate, startmonth);
            [trendm_dRsfc_ts, trends_dRsfc_ts, trendyr_dRsfc_ts, ~, ~] = autoCalTrend(tsEffect, nlonf, nlatf, timeDate, startmonth);

            readme.aboutTas2 = 'tas2 mean consider 2 levels near sfc';
            save([trend_radEfectPath, 'global_vars.mat'], 'lon_f', 'lat_f', 'time', 'plev_k', 'readme', 'timeseries', 'modelname')
            save([trend_radEfectPath, 'trend_dradEfect_sfc_cld.mat'], ...
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

            %%%%Part2: cal the model output net radEffect
            load([radEfectPath, 'global_vars.mat'])% lat_f lon_f time plev_k readme
            load([radEfectPath, 'real_dradEfect.mat'])%'l_rad', 's_rad', 'dR_allsky', 'dR_clr', 'readme_realradEfect'
            realRadEffect = dR_allsky(:, :, :, 1);
            nlonf = length(lon_f); nlatf = length(lat_f);
            % cal the trend(10 vars)
            [trendm_dRsfc_realRad, trends_dRsfc_realRad, trendyr_dRsfc_realRad, ~, ~] = autoCalTrend(realRadEffect, nlonf, nlatf, timeDate, startmonth);
            save([trend_radEfectPath, 'trend_drealRadEfect_sfc_cld.mat'], ...
                'trendm_dRsfc_realRad', 'trends_dRsfc_realRad', 'trendyr_dRsfc_realRad')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % cloud TOA rad trend
            %%%%Part1: cal radEffect
            load([radEfectPath, 'global_vars.mat'])% lat_f lon_f time plev_k readme
            load([radEfectPath, 'dCRF.mat'])% dCRF
            load([radEfectPath, 'dradEfect_toa_cld.mat'])%'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', 'taEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect'
            load([radEfectPath, 'dR_mainEffect_toa.mat'])%dR_mainEffect_toa (ta+alb+q+clod effect)
            load([radEfectPath, 'dR_residual_cld_toa.mat'])%dR_residual_cld_toa(co2 effect)
            load([radEfectPath, 'dR_cloud_toa.mat'])%dR_residual_cld_toa(co2 effect)
            dCRF_toa = squeeze(dCRF(:, :, :, 2));
            nlonf = length(lon_f); nlatf = length(lat_f);
            timeDate = time.date;
            % cal the trend(10 vars)
            [trendm_dRtoa_CRF, trends_dRtoa_CRF, trendyr_dRtoa_CRF, ~, ~] = autoCalTrend(dCRF_toa, nlonf, nlatf, timeDate, startmonth);
            [trendm_dRtoa_ta, trends_dRtoa_ta, trendyr_dRtoa_ta, ~, ~] = autoCalTrend(taEffect, nlonf, nlatf, timeDate, startmonth);
            [trendm_dRtoa_taOnly2, trends_dRtoa_taOnly2, trendyr_dRtoa_taOnly2, ~, ~] = autoCalTrend(taOnlyEffect2, nlonf, nlatf, timeDate, startmonth);
            [trendm_dRtoa_tas2, trends_dRtoa_tas2, trendyr_dRtoa_tas2, ~, ~] = autoCalTrend(tasEffect2, nlonf, nlatf, timeDate, startmonth);
            [trendm_dRtoa_mainEffect, trends_dRtoa_mainEffect, trendyr_dRtoa_mainEffect, ~, ~] = autoCalTrend(dR_mainEffect_toa, nlonf, nlatf, timeDate, startmonth);
            [trendm_dRtoa_residual, trends_dRtoa_residual, trendyr_dRtoa_residual, ~, ~] = autoCalTrend(dR_residual_cld_toa, nlonf, nlatf, timeDate, startmonth);
            [trendm_dRtoa_cloud, trends_dRtoa_cloud, trendyr_dRtoa_cloud, ~, ~] = autoCalTrend(dR_cloud_toa, nlonf, nlatf, timeDate, startmonth);
            [trendm_dRtoa_q, trends_dRtoa_q, trendyr_dRtoa_q, ~, ~] = autoCalTrend(husEffect, nlonf, nlatf, timeDate, startmonth);
            [trendm_dRtoa_alb, trends_dRtoa_alb, trendyr_dRtoa_alb, ~, ~] = autoCalTrend(albEffect, nlonf, nlatf, timeDate, startmonth);
            [trendm_dRtoa_ts, trends_dRtoa_ts, trendyr_dRtoa_ts, ~, ~] = autoCalTrend(tsEffect, nlonf, nlatf, timeDate, startmonth);
            readme.aboutTas2 = 'tas2 mean consider 2 levels near sfc';
            save([trend_radEfectPath, 'global_vars.mat'], 'lon_f', 'lat_f', 'time', 'plev_k', 'readme', 'timeseries', 'modelname')
            save([trend_radEfectPath, 'trend_dradEfect_toa_cld.mat'], ...
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

            %%%%Part2: cal the model output net radEffect
            load([radEfectPath, 'global_vars.mat'])% lat_f lon_f time plev_k readme
            load([radEfectPath, 'real_dradEfect.mat'])%'l_rad', 's_rad', 'dR_allsky', 'dR_clr', 'readme_realradEfect'
            realRadEffect = dR_allsky(:, :, :, 2);
            nlonf = length(lon_f); nlatf = length(lat_f);
            % cal the trend(10 vars)
            [trendm_dRtoa_realRad, trends_dRtoa_realRad, trendyr_dRtoa_realRad, ~, ~] = autoCalTrend(realRadEffect, nlonf, nlatf, timeDate, startmonth);
            save([trend_radEfectPath, 'trend_drealRadEfect_toa_cld.mat'], ...
                'trendm_dRtoa_realRad', 'trends_dRtoa_realRad', 'trendyr_dRtoa_realRad')

            disp([esmName{esmNum, 1}, ' ensemble is done!'])
        end

        disp([level.model2{level1}, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{p_1}, ' era is done!'])
end

t = toc; disp(t)
