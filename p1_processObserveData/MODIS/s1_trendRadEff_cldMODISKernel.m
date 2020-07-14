%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-07-13 19:04:55
% LastEditTime : 2020-07-14 14:57:30
% LastEditors  : LYC
% Description  :
% FilePath     : /code/p1_processObserveData/MODIS/s1_trendRadEff_cldMODISKernel.m
%
%%---------------------------------------------------------
% cal 200207_201706 cloud Effect trend by MODIS data using kernel method
% raw data time series: 200207_201706
%
clc; clear; tic;
lon_f = 0:2.5:357.5; nlonf = length(lon_f);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f);

p_1 = 2; % different time series, 1mean 2000 - 03 to 2018 - 02(18 * 12). 2 mean 200207 - 201706(15 * 12)
[readme, level, tLin, vars] = obsParameters('ERA5');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part1: cal MODIS cloud effect and trend
% read data
rawDataPath = '/data1/liuyincheng/Observe-rawdata/';
anomPath_mod = '/data1/liuyincheng/Observe-process/200207-201706/MODIS/radEffect/';
anomTrendPath_mod = '/data1/liuyincheng/Observe-process/200207-201706/MODIS/radEffect_trend/';
auto_mkdir(anomPath_mod); auto_mkdir(anomTrendPath_mod)

load([rawDataPath, 'surface_flux_MODIS_kernel_estimation_200207_201908.mat'])% lat_m, lon_m, lwcf, swcf,, lat1,lon1,
% coordinate
lat_m = double(lat); lon_m = double(lon);
lon_m(2:end + 1) = lon_m; lon_m(1) = -0.5;
t0 = datetime(2002, 7, 1);
t = t0 + calmonths(0:179);
time = datenum(t);

Rsfc_cloud = swcf + lwcf;
Rsfc_cloud(:, 2:end + 1, :) = Rsfc_cloud; Rsfc_cloud(:, 1, :) = Rsfc_cloud(:, end, :);
Rsfc_cloud = permute(Rsfc_cloud, [2 1 3]);
% choose 200207_201706 time
Rsfc_cloud = Rsfc_cloud(:, :, 1:tLin.inter{p_1});

% Regrid
Rsfc_cloud = autoRegrid3(lon_m, lat_m, time, Rsfc_cloud, lon_f, lat_f, time);
% Deseasonalize
startMonth = tLin.startMonth{2};
[dRMOsfc_cloud, Clim_dRMOsfc_cloud] = monthlyAnomaly3D(nlonf, nlatf, time, Rsfc_cloud, startMonth);
% cal trend
[trendm_dRMOsfc_cloud, trends_dRMOsfc_cloud, trendyr_dRMOsfc_cloud, ~, ~] = autoCalTrend(dRMOsfc_cloud, nlonf, nlatf, time, startMonth);
save([anomPath_mod, 'global_vars.mat'], 'lon_f', 'lat_f', 'time');
save([anomPath_mod, 'dRMO_cloud.mat'], 'dRMOsfc_cloud', 'Clim_dRMOsfc_cloud');
save([anomTrendPath_mod, 'trend_dRMO_cloud.mat'], 'trendm_dRMOsfc_cloud', 'trends_dRMOsfc_cloud', 'trendyr_dRMOsfc_cloud');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part2: cal MODIS and ERA rad data
%%  era data ------------------------
sfcToaName = {'sfc', 'toa'};
allClrName = {'all', 'clr'};
eraName = {'ERA5', 'ERAi'};

%% Read date
for eraNum = 1:2% 1 mean ERA5, 2 mean ERAi
    anomPath_era = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, eraName{eraNum}, level.standVarPath{2}); % /anomaly
    radEfectPath_era = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, eraName{eraNum}, level.standVarPath{5}); %radEffect
    load([radEfectPath_era, 'dradEfect_union.mat']); %%'dR_hus', 'dR_alb', 'dR_ts', 'dR_ta', 'dR_total', 'readme_radEfect'
    load([radEfectPath_era, 'real_dradEfect.mat'])%'l_rad', 's_rad', 'dR_allsky', 'dR_clr', 'readme_realradEfect'
    %( 1sfc cld,2sfc clr,3toa cld,4toa clr)
    dR_main = dR_hus + dR_alb + dR_ta; % dR exclude ts
    dRMO_mainEffect_sfc = dRMOsfc_cloud + squeeze(dR_main(:, :, :, 1));
    dRMO_res_sfc = squeeze(dR_allsky(:, :, :, 1)) - dRMO_mainEffect_sfc - squeeze(dR_ts(:, :, :, 1));
    dRMO_mainEffect_sfc_addRes = dRMO_mainEffect_sfc + dRMO_res_sfc;

    % cal trend
    [trendm_dRMO_mainEffect_sfc, trends_dRMO_mainEffect_sfc, trendyr_dRMO_mainEffect_sfc, ~, ~] = autoCalTrend(dRMO_mainEffect_sfc, nlonf, nlatf, time, startMonth);
    [trendm_dRMO_mainEffect_sfc_addRes, trends_dRMO_mainEffect_sfc_addRes, trendyr_dRMO_mainEffect_sfc_addRes, ~, ~] = autoCalTrend(dRMO_mainEffect_sfc_addRes, nlonf, nlatf, time, startMonth);
    [trendm_dRMO_res_sfc, trends_dRMO_res_sfc, trendyr_dRMO_res_sfc, ~, ~] = autoCalTrend(dRMO_res_sfc, nlonf, nlatf, time, startMonth);

    if eraNum == 1
        save([anomPath_mod, 'dRMO_mainEffect_sfc_era5.mat'], 'dRMO_mainEffect_sfc', 'dRMO_mainEffect_sfc_addRes');
        save([anomPath_mod, 'dRMO_res_sfc_era5.mat'], 'dRMO_res_sfc');

        save([anomTrendPath_mod, 'global_vars.mat'], 'lon_f', 'lat_f', 'time');
        save([anomTrendPath_mod, 'trend_dRMOD_cloud_era5.mat'], 'trendm_dRMO_mainEffect_sfc', 'trends_dRMO_mainEffect_sfc', 'trendyr_dRMO_mainEffect_sfc', ...
            'trendm_dRMO_mainEffect_sfc_addRes', 'trends_dRMO_mainEffect_sfc_addRes', 'trendyr_dRMO_mainEffect_sfc_addRes', ...
            'trendm_dRMO_res_sfc', 'trends_dRMO_res_sfc', 'trendyr_dRMO_res_sfc')
    elseif eraNum == 2
        save([anomPath_mod, 'dRMO_mainEffect_sfc_erai.mat'], 'dRMO_mainEffect_sfc', 'dRMO_mainEffect_sfc_addRes');
        save([anomPath_mod, 'dRMO_res_sfc_erai.mat'], 'dRMO_res_sfc');

        save([anomTrendPath_mod, 'global_vars.mat'], 'lon_f', 'lat_f', 'time');
        save([anomTrendPath_mod, 'trend_dRMOD_cloud_erai.mat'], 'trendm_dRMO_mainEffect_sfc', 'trends_dRMO_mainEffect_sfc', 'trendyr_dRMO_mainEffect_sfc', ...
            'trendm_dRMO_mainEffect_sfc_addRes', 'trends_dRMO_mainEffect_sfc_addRes', 'trendyr_dRMO_mainEffect_sfc_addRes', ...
            'trendm_dRMO_res_sfc', 'trends_dRMO_res_sfc', 'trendyr_dRMO_res_sfc')
    end

end

t = toc; disp(t)
