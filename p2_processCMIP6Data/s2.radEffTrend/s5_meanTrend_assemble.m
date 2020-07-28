%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-07-28 10:52:10
% LastEditTime : 2020-07-28 10:57:13
% LastEditors  : LYC
% Description  : 
% FilePath     : /code/p2_processCMIP6Data/s2.radEffTrend/s5_meanTrend_assemble.m
%  
%%---------------------------------------------------------
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% mothly data
% plot the Ts trend and rhs,Ts, Ta, wv, albedo trend(TOA)
% note: the trendyr already be masked and caled cc
%
% experiment information:
% time:2000.01-2014.12(interval:15*12);1980.01-2014.12(interval:35*12); 2015.01-2099.12(interval:85*12)
% initial time in hist(1740 total): 1,561 of 1740(2000.03);1,321 of 1740(1980.01)
% initial time in futrue(1032 total): 1 of 1032(2015.01);
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clc; clear;
nowpath = pwd;
% load mask map
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')
[mlabels, areaNum] = cmipPlotParameters('sfc', 'land', 'radEffect'); % plot parameters
esm = 'r1i1p1f1';
%% global set
% Latitude range
latRange = 60;
lon1 = [2.5 357.5]; lat1 = [-latRange + 1 latRange - 1]; % world area
set(0, 'defaultfigurecolor', 'w')
% Experent ID
exm_left = 1; exm_right = 4;

for exmName = exm_left:exm_right
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmName);
    % mPath.input:E:/data/cmip6-process/2000-2014/
    mPath.input = fullfile('/data1/liuyincheng/cmip6-process/', level.time1{exmName});
    % mPath.output:a_research/P02.Ts_change_research/figure/04.cmip6Result/2000-2014/
    mPath.uniOutput = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/02.cmip6Result/1.6/', lower(mlabels.level), level.time1{exmName});
    mPath.Output = fullfile(mPath.uniOutput);
    auto_mkdir(mPath.Output)

    % model loop
    mdl_left = 1; mdl_right = length(level.model2);% differnt models%length(level.model2)

    for mdlName = mdl_left:mdl_right
        mdlPath = fullfile(mPath.input, level.model2{mdlName});

        % load data
        dradTrendPath = fullfile(mdlPath, esm, level.process3{7}); %/data1/liuyincheng/cmip6-process/amip_1980-2014/CESM2/Effect_trend
        danomTrendPath = fullfile(mdlPath, esm, level.process3{3}); %/data1/liuyincheng/cmip6-process/amip_1980-2014/CESM2/anomaly_trend

        if ~exist(dradTrendPath, 'dir')
            disp(['the ', esm, ' ensemble of ', level.model2{mdlName}, ' didnt exist']);
            continue
        end

        load([dradTrendPath, 'global_vars.mat'])% lat_f lon_f time plevf readme
        load([dradTrendPath, 'trend_dradEfect_toa_cld.mat'])% 10 vars:'trendyr_dRtoa_ta','trendyr_dRtoa_taOnly2', 'trendyr_dRtoa_tas2., 'trendyr_dRtoa_tsAtom', 'trendyr_dRtoa_mainEffect', 'trendyr_dRtoa_residual', 'trendyr_dRtoa_cloud', 'trendyr_dRtoa_q', 'trendyr_dRtoa_alb', 'trendyr_dRtoa_ts'
        load([dradTrendPath, 'trend_dradEfect_sfc_cld.mat'])% 10 vars:'trendyr_dRsfc_ta','trendyr_dRsfc_taOnly2', 'trendyr_dRsfc_tas2., 'trendyr_dRsfc_tsAtom', 'trendyr_dRsfc_mainEffect', 'trendyr_dRsfc_residual', 'trendyr_dRsfc_cloud', 'trendyr_dRsfc_q', 'trendyr_dRsfc_alb', 'trendyr_dRsfc_ts'
        load([danomTrendPath, 'trend_dnetTOA.mat'])% trendyr_dnetTOA
        load([danomTrendPath, 'trend_drhs.mat'])% trendyr_drhs
        load([danomTrendPath, 'trend_dhFlux.mat'])% trendyr_dhFlux
        load([danomTrendPath, 'trend_drhsPlus.mat'])% trendyr_dhFlux
        load([danomTrendPath, 'trend_dts.mat'])% trendyr_dts
        nlonf = length(lon_f); nlatf = length(lat_f);
        trendyr_dRatm_ta=trendyr_dRtoa_ta-trendyr_dRsfc_ta;
        trendyr_dRatm_tas2=trendyr_dRtoa_tas2-trendyr_dRsfc_tas2;
        trendyr_dRatm_taOnly2=trendyr_dRtoa_taOnly2-trendyr_dRsfc_taOnly2;
        trendyr_dRatm_cloud=trendyr_dRtoa_cloud-trendyr_dRsfc_cloud;
        trendyr_dRatm_q=trendyr_dRtoa_q-trendyr_dRsfc_q;
        trendyr_dRatm_alb=trendyr_dRtoa_alb-trendyr_dRsfc_alb;
        trendyr_dRatm_mainEffect=trendyr_dRtoa_mainEffect-trendyr_dRsfc_mainEffect;
        trendyr_dRatm_residual=trendyr_dRtoa_residual-trendyr_dRsfc_residual;

        % trendyr_dRatm_tsAtom=-trendyr_dRatm_tsAtom; % dim downward rad as positive
        % use one var to plot
        trendyr = zeros(nlonf, nlatf, 12);

        for varNum = 1:12
            eval(['trendyr(:,:,varNum)=', mlabels.vars{varNum}, ';'])
        end

        trendyr = trendyr * 365 * 10;
        trendyr(:, :, 1) = trendyr(:, :, 1);
        trendyr(:, :, 10) = trendyr(:, :, 10);
        % mask and cal the cc
        [trendyr, yr_cc, yr_pp] = maskArea(trendyr, lat_f, latRange, -latRange, areaNum);

 


    %model
    end
    %exm
end


