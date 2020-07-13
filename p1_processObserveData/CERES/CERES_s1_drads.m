%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-07-08 10:22:15
% LastEditTime : 2020-07-13 10:10:15
% LastEditors  : LYC
% Description  :
% FilePath     : /code/p1_processObserveData/CERES/CERES_s1_drads.m
%
%%---------------------------------------------------------
% This program is to calculate the radiation anomalies from CERES/EBAF4.0
clc; clear; tic;
% kernel Coordinate
lon_k = 0:2.5:357.5; nlonk = length(lon_k);
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
% plot figure Coordinate
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f);
lon_f = lon_k; nlonf = length(lon_f);
% ceres Coordinate
lon_c = 0.5:1:359.5; nlonc = length(lon_c); lon_c = double(lon_c);
lat_c = -89.5:1:89.5; nlatc = length(lat_c); lat_c = double(lat_c);

vector_var_str{1} = 'sfc';
vector_var_str{2} = 'toa';
%modify path first
ceresDataPath = '/data1/liuyincheng/Observe-rawdata/CERES_EBAF/4.1/';
sfcDataPath = strcat(ceresDataPath, 'CERES_EBAF-SFC_Ed4.1_Subset_200003-201911.nc');
toaDataPath = strcat(ceresDataPath, 'CERES_EBAF-TOA_Ed4.1_Subset_200003-201911.nc');

% transfor mat time
time = double(ncread(sfcDatePath, 'time')); % days since 2000-03-01 00:00:00
timeUnits = ncreadatt(sfcDatePath, 'time', 'units');
timeCalendar = 'gregorian';
time = cdfdate2num(timeUnits, timeCalendar, time);
ntime = length(time);

%% read raw data
%% SFC
% All-sky
dR_tot_sfc_all = ncread(sfcDatePath, 'sfc_net_tot_all_mon'); % = sfc_net_lw_all_mon + sfc_net_sw_all_mon
dR_cre_sfc_all_t = ncread(sfcDatePath, 'sfc_cre_net_tot_mon'); % Cloud radiative effect(CRF: sfc_net_tot_all_mon-sfc_net_tot_clr_mon_t)
dR_lwDown_sfc_all = ncread(sfcDatePath, 'sfc_lw_down_all_mon'); % LW downwards
dR_swNet_sfc_all = ncread(sfcDatePath, 'sfc_net_sw_all_mon'); % surface_net_downward_shortwave_flux
% Clear-sky
dR_tot_sfc_clr_t = ncread(sfcDatePath, 'sfc_net_tot_clr_mon_t'); % transfer model result
dR_tot_sfc_clr_c = ncread(sfcDatePath, 'sfc_net_tot_clr_mon_c'); % obs for cloud-free areas result

%% TOA
% allsky
dR_tot_toa_all = ncread(toaDataPath, 'toa_net_all_mon'); % = - toa_sw_all_mon - toa_lw_all_mon + solar_mon
dR_cre_toa_all_t = ncread(toaDataPath, 'toa_cre_net_mon'); % Cloud radiative effect
% clear sky
dR_tot_toa_clr_t = ncread(toaDataPath, 'toa_net_clr_mon_t');
dR_tot_toa_clr_c = ncread(toaDataPath, 'toa_net_clr_mon_c');

month = 1:12;
%% Regrid (Because of unite plot)
% For interpolate(??????????lonf???lon_c?????(lon_c????????), latf???????lat_c????????????????)
lon_c(2:end + 1) = lon_c;
lon_c(1) = -0.5;

dR_tot_sfc_all = transfCordi(lon_c, lat_c, time, dR_tot_sfc_all, lon_f, lat_f, time);
dR_cre_sfc_all_t = transfCordi(lon_c, lat_c, time, dR_cre_sfc_all_t, lon_f, lat_f, time);
dR_lwDown_sfc_all = transfCordi(lon_c, lat_c, time, dR_lwDown_sfc_all, lon_f, lat_f, time);
dR_swNet_sfc_all = transfCordi(lon_c, lat_c, time, dR_swNet_sfc_all, lon_f, lat_f, time);
dR_tot_sfc_clr_t = transfCordi(lon_c, lat_c, time, dR_tot_sfc_clr_t, lon_f, lat_f, time);
dR_tot_sfc_clr_c = transfCordi(lon_c, lat_c, time, dR_tot_sfc_clr_c, lon_f, lat_f, time);

dR_tot_toa_all = transfCordi(lon_c, lat_c, time, dR_tot_toa_all, lon_f, lat_f, time);
dR_cre_toa_all_t = transfCordi(lon_c, lat_c, time, dR_cre_toa_all_t, lon_f, lat_f, time);
dR_tot_toa_clr_t = transfCordi(lon_c, lat_c, time, dR_tot_toa_clr_t, lon_f, lat_f, time);
dR_tot_toa_clr_c = transfCordi(lon_c, lat_c, time, dR_tot_toa_clr_c, lon_f, lat_f, time);

%% different time series, 1mean 2000-03 to 2018-02(18*12). 2 mean 200207-201706(15*12)
[readme, level, tLin, vars] = obsParameters('CERES');

for p_1 = 1:2
    anomPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'CERES', level.standVarPath{2}); % /anomaly
    anomPath_era5 = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{2}); % /anomaly
    anomPath_erai = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi', level.standVarPath{2}); % /anomaly
    anomTrendPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'CERES', level.standVarPath{3}); % /anomaly_trend
    radEfectPath_era5 = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{5}); %radEffect
    radEfectPath_erai = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi', level.standVarPath{5}); %radEffect
    auto_mkdir(anomPath); auto_mkdir(anomTrendPath);

    % find start and end month
    timeSpec = tLin.start{p_1};
    formatOut = 'yyyy-mm';
    timeStr = string(datestr(time, formatOut));
    startT = find(timeStr == timeSpec);
    endT = startT + tLin.inter{p_1} - 1;

    % read specific time series
    dR_tot_sfc_all = dR_tot_sfc_all(:, :, startT:endT);
    dR_cre_sfc_all_t = dR_cre_sfc_all_t(:, :, startT:endT);
    dR_lwDown_sfc_all = dR_lwDown_sfc_all(:, :, startT:endT);
    dR_swNet_sfc_all = dR_swNet_sfc_all(:, :, startT:endT);
    dR_tot_sfc_clr_t = dR_tot_sfc_clr_t(:, :, startT:endT);
    dR_tot_sfc_clr_c = dR_tot_sfc_clr_c(:, :, startT:endT);

    dR_tot_toa_all = dR_tot_toa_all(:, :, startT:endT);
    dR_cre_toa_all_t = dR_cre_toa_all_t(:, :, startT:endT);
    dR_tot_toa_clr_t = dR_tot_toa_clr_t(:, :, startT:endT);
    dR_tot_toa_clr_c = dR_tot_toa_clr_c(:, :, startT:endT);

    dR_cre_sfc_all_c = dR_tot_sfc_all - dR_tot_sfc_clr_c;
    dR_cre_toa_all_c = dR_tot_toa_all - dR_tot_toa_clr_c;
    time = time(startT:endT, 1);
    ntime = length(time);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Part1: deseasonlize vars
    startMonth = tLin.startMonth{p_1};
    [dR_tot_sfc_all, Clim_dR_tot_sfc_all] = monthlyAnomaly3D(nlonf, nlatf, time, dR_tot_sfc_all, startMonth);
    [dR_cre_sfc_all_t, Clim_dR_cre_sfc_all_t] = monthlyAnomaly3D(nlonf, nlatf, time, dR_cre_sfc_all_t, startMonth);
    [dR_cre_sfc_all_c, Clim_dR_cre_sfc_all_c] = monthlyAnomaly3D(nlonf, nlatf, time, dR_cre_sfc_all_c, startMonth);
    [dR_lwDown_sfc_all, Clim_dR_lwDown_sfc_all] = monthlyAnomaly3D(nlonf, nlatf, time, dR_lwDown_sfc_all, startMonth);
    [dR_swNet_sfc_all, Clim_dR_swNet_sfc_all] = monthlyAnomaly3D(nlonf, nlatf, time, dR_swNet_sfc_all, startMonth);
    [dR_tot_sfc_clr_t, Clim_dR_tot_sfc_clr_t] = monthlyAnomaly3D(nlonf, nlatf, time, dR_tot_sfc_clr_t, startMonth);
    [dR_tot_sfc_clr_c, Clim_dR_tot_sfc_clr_c] = monthlyAnomaly3D(nlonf, nlatf, time, dR_tot_sfc_clr_c, startMonth);

    [dR_tot_toa_all, Clim_dR_tot_toa_all] = monthlyAnomaly3D(nlonf, nlatf, time, dR_tot_toa_all, startMonth);
    [dR_cre_toa_all_t, Clim_dR_cre_toa_all_t] = monthlyAnomaly3D(nlonf, nlatf, time, dR_cre_toa_all_t, startMonth);
    [dR_cre_toa_all_c, Clim_dR_cre_toa_all_c] = monthlyAnomaly3D(nlonf, nlatf, time, dR_cre_toa_all_c, startMonth);
    [dR_tot_toa_clr_t, Clim_dR_tot_toa_clr_t] = monthlyAnomaly3D(nlonf, nlatf, time, dR_tot_toa_clr_t, startMonth);
    [dR_tot_toa_clr_c, Clim_dR_tot_toa_clr_c] = monthlyAnomaly3D(nlonf, nlatf, time, dR_tot_toa_clr_c, startMonth);

    %% save mat files
    save([anomPath, 'global_vars.mat'], 'lon_f', 'lat_f', 'time', 'readme')
    save([anomPath, 'dCErad.mat'], ...
        'dR_tot_sfc_all', 'dR_cre_sfc_all_t', 'dR_cre_sfc_all_c', 'dR_lwDown_sfc_all', 'dR_swNet_sfc_all', 'dR_tot_sfc_clr_t', 'dR_tot_sfc_clr_c', ...
        'Clim_dR_tot_sfc_all', 'Clim_dR_cre_sfc_all_t', 'Clim_dR_cre_sfc_all_c', 'Clim_dR_lwDown_sfc_all', 'Clim_dR_swNet_sfc_all', 'Clim_dR_tot_sfc_clr_t', 'Clim_dR_tot_sfc_clr_c', ...
        'dR_tot_toa_all', 'dR_cre_toa_all_t', 'dR_cre_toa_all_c', 'dR_tot_toa_clr_t', 'dR_tot_toa_clr_c', ...
        'Clim_dR_tot_toa_all', 'Clim_dR_cre_toa_all_t', 'Clim_dR_cre_toa_all_c', 'Clim_dR_tot_toa_clr_t', 'Clim_dR_tot_toa_clr_c');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Part2: cal specific rad componet (rhs, rhsPlus)
    % SFC
    % cal dCErhs
    dCErhs = dR_lwDown_sfc_all + dR_swNet_sfc_all;
    [trendm_dCErhs, trends_dCErhs, trendyr_dCErhs, ~, ~] = autoCalTrend(dCErhs, nlonf, nlatf, time, startMonth);
    save([anomPath, 'dCErhs.mat'], 'dCErhs');
    save([anomTrendPath, 'global_vars.mat'], 'lon_f', 'lat_f', 'time')
    save([anomTrendPath, 'trend_dCErhs.mat'], 'trendm_dCErhs', 'trends_dCErhs', 'trendyr_dCErhs');
    % cal the drhsPlus
    load([anomPath_era5, 'drhs.mat'], 'dhFlux')% ERA5 dhFlux
    dCErhsPlus_era5 = dCErhs + dhFlux;
    [trendm_dCErhsPlus_era5, trends_dCErhsPlus_era5, trendyr_dCErhsPlus_era5, ~, ~] = autoCalTrend(dCErhsPlus_era5, nlonf, nlatf, time, startMonth);
    save([anomPath, 'dCErhsPlus_era5.mat'], 'dCErhsPlus_era5');
    save([anomTrendPath, 'trend_dCErhsPlus_era5.mat'], 'trendm_dCErhsPlus_era5', 'trends_dCErhsPlus_era5', 'trendyr_dCErhsPlus_era5');

    load([anomPath_erai, 'drhs.mat'], 'dhFlux')% ERAi dhFlux
    dCErhsPlus_erai = dCErhs + dhFlux;
    [trendm_dCErhsPlus_erai, trends_dCErhsPlus_erai, trendyr_dCErhsPlus_erai, ~, ~] = autoCalTrend(dCErhsPlus_erai, nlonf, nlatf, time, startMonth);
    save([anomPath, 'dCErhsPlus_erai.mat'], 'dCErhsPlus_erai');
    save([anomTrendPath, 'trend_dCErhsPlus_erai.mat'], 'trendm_dCErhsPlus_erai', 'trends_dCErhsPlus_erai', 'trendyr_dCErhsPlus_erai');
    % TOA
    dCEnetTOA = dR_tot_toa_all;
    [trendm_dCEnetTOA, trends_dCEnetTOA, trendyr_dCEnetTOA, p_dCEnetTOA, cons_dCEnetTOA] = autoCalTrend(dCEnetTOA, nlonf, nlatf, time, startMonth);
    save([anomTrendPath, 'trend_dCEnetTOA.mat'], 'trendm_dCEnetTOA', 'cons_dCEnetTOA', 'p_dCEnetTOA', 'trends_dCEnetTOA', 'trendyr_dCEnetTOA');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Part3: cal cloud Effect and trend
    %% ERA5 ,  _t
    % load  ( 1sfc cld,2sfc clr,3toa cld,4toa clr)
    load([radEfectPath_era5, 'dradEfect_union.mat']); %'dR_hus', 'dR_alb', 'dR_ts', 'dR_ta', 'dR_total', 'readme_radEfect'
    dCE_CRF_t(:, :, :, 1) = dR_cre_sfc_all_t;
    dCE_CRF_t(:, :, :, 2) = dR_cre_toa_all_t;
    dvarsFeedback = dR_total(:, :, :, [2 4]) - dR_total(:, :, :, [1 3]);
    dRCE_cloud = dCE_CRF_t + dvarsFeedback;
    dRCE_clr(:, :, :, 1) = dR_tot_sfc_clr_t;
    dRCE_clr(:, :, :, 2) = dR_tot_toa_clr_t;
    dRCE_residual = dRCE_clr - dR_total(:, :, :, [2 4]);
    %save as one file
    dCE_era5_t.Rcloud_sfc = squeeze(dRCE_cloud(:, :, :, 1));
    dCE_era5_t.Rcloud_toa = squeeze(dRCE_cloud(:, :, :, 2));
    dCE_era5_t.residual_sfc = squeeze(dRCE_residual(:, :, :, 1));
    dCE_era5_t.residual_toa = squeeze(dRCE_residual(:, :, :, 2));
    % _c
    dCE_CRF_c(:, :, :, 1) = dR_cre_sfc_all_c;
    dCE_CRF_c(:, :, :, 2) = dR_cre_toa_all_c;
    dvarsFeedback = dR_total(:, :, :, [2 4]) - dR_total(:, :, :, [1 3]);
    dRCE_cloud = dCE_CRF_c + dvarsFeedback;
    dRCE_clr(:, :, :, 1) = dR_tot_sfc_clr_c;
    dRCE_clr(:, :, :, 2) = dR_tot_toa_clr_c;
    dRCE_residual = dRCE_clr - dR_total(:, :, :, [2 4]);
    %save as one file
    dCE_era5_c.Rcloud_sfc = squeeze(dRCE_cloud(:, :, :, 1));
    dCE_era5_c.Rcloud_toa = squeeze(dRCE_cloud(:, :, :, 2));
    dCE_era5_c.residual_sfc = squeeze(dRCE_residual(:, :, :, 1));
    dCE_era5_c.residual_toa = squeeze(dRCE_residual(:, :, :, 2));
    [trendm_dCEcloud_sfc_t, trends_dCEcloud_sfc_t, trendyr_dCEcloud_sfc_t, ~, ~] = autoCalTrend(dCE_era5_t.Rcloud_sfc, nlonf, nlatf, time, startMonth);
    [trendm_dCEcloud_toa_t, trends_dCEcloud_toa_t, trendyr_dCEcloud_toa_t, ~, ~] = autoCalTrend(dCE_era5_t.Rcloud_toa, nlonf, nlatf, time, startMonth);
    [trendm_dCEcloud_sfc_c, trends_dCEcloud_sfc_c, trendyr_dCEcloud_sfc_c, ~, ~] = autoCalTrend(dCE_era5_c.Rcloud_sfc, nlonf, nlatf, time, startMonth);
    [trendm_dCEcloud_toa_c, trends_dCEcloud_toa_c, trendyr_dCEcloud_toa_c, ~, ~] = autoCalTrend(dCE_era5_c.Rcloud_toa, nlonf, nlatf, time, startMonth);
    % save
    save([anomPath, 'dCEcloud_era5.mat'], 'dCR_ear5_t', 'dCR_ear5_c')
    save([anomTrendPath, 'trend_dCEcloud_era5.mat'], 'trendm_dCEcloud_sfc_t', 'trends_dCEcloud_sfc_t', 'trendyr_dCEcloud_sfc_t', ...
        'trendm_dCEcloud_toa_t', 'trends_dCEcloud_toa_t', 'trendyr_dCEcloud_toa_t', ...
        'trendm_dCEcloud_sfc_c', 'trends_dCEcloud_sfc_c', 'trendyr_dCEcloud_sfc_c', ...
        'trendm_dCEcloud_toa_c', 'trends_dCEcloud_toa_c', 'trendyr_dCEcloud_toa_c');

    %% ERAi ,  _t
    % load  ( 1sfc cld,2sfc clr,3toa cld,4toa clr)
    load([radEfectPath_erai, 'dradEfect_union.mat']); %'dR_hus', 'dR_alb', 'dR_ts', 'dR_ta', 'dR_total', 'readme_radEfect'
    dCE_CRF_t(:, :, :, 1) = dR_cre_sfc_all_t;
    dCE_CRF_t(:, :, :, 2) = dR_cre_toa_all_t;
    dvarsFeedback = dR_total(:, :, :, [2 4]) - dR_total(:, :, :, [1 3]);
    dRCE_cloud = dCE_CRF_t + dvarsFeedback;
    dRCE_clr(:, :, :, 1) = dR_tot_sfc_clr_t;
    dRCE_clr(:, :, :, 2) = dR_tot_toa_clr_t;
    dRCE_residual = dRCE_clr - dR_total(:, :, :, [2 4]);
    %save as one file
    dCE_erai_t.Rcloud_sfc = squeeze(dRCE_cloud(:, :, :, 1));
    dCE_erai_t.Rcloud_toa = squeeze(dRCE_cloud(:, :, :, 2));
    dCE_erai_t.residual_sfc = squeeze(dRCE_residual(:, :, :, 1));
    dCE_erai_t.residual_toa = squeeze(dRCE_residual(:, :, :, 2));
    % _c
    dCE_CRF_c(:, :, :, 1) = dR_cre_sfc_all_c;
    dCE_CRF_c(:, :, :, 2) = dR_cre_toa_all_c;
    dvarsFeedback = dR_total(:, :, :, [2 4]) - dR_total(:, :, :, [1 3]);
    dRCE_cloud = dCE_CRF_c + dvarsFeedback;
    dRCE_clr(:, :, :, 1) = dR_tot_sfc_clr_c;
    dRCE_clr(:, :, :, 2) = dR_tot_toa_clr_c;
    dRCE_residual = dRCE_clr - dR_total(:, :, :, [2 4]);
    %save as one file
    dCE_erai_c.Rcloud_sfc = squeeze(dRCE_cloud(:, :, :, 1));
    dCE_erai_c.Rcloud_toa = squeeze(dRCE_cloud(:, :, :, 2));
    dCE_erai_c.residual_sfc = squeeze(dRCE_residual(:, :, :, 1));
    dCE_erai_c.residual_toa = squeeze(dRCE_residual(:, :, :, 2));
    [trendm_dCEcloud_sfc_t, trends_dCEcloud_sfc_t, trendyr_dCEcloud_sfc_t, ~, ~] = autoCalTrend(dCE_erai_t.Rcloud_sfc, nlonf, nlatf, time, startMonth);
    [trendm_dCEcloud_toa_t, trends_dCEcloud_toa_t, trendyr_dCEcloud_toa_t, ~, ~] = autoCalTrend(dCE_erai_t.Rcloud_toa, nlonf, nlatf, time, startMonth);
    [trendm_dCEcloud_sfc_c, trends_dCEcloud_sfc_c, trendyr_dCEcloud_sfc_c, ~, ~] = autoCalTrend(dCE_erai_c.Rcloud_sfc, nlonf, nlatf, time, startMonth);
    [trendm_dCEcloud_toa_c, trends_dCEcloud_toa_c, trendyr_dCEcloud_toa_c, ~, ~] = autoCalTrend(dCE_erai_c.Rcloud_toa, nlonf, nlatf, time, startMonth);
    % save
    save([anomPath, 'dCEcloud_erai.mat'], 'dCR_eari_t', 'dCR_eari_c')
    save([anomTrendPath, 'trend_dCEcloud_erai.mat'], 'trendm_dCEcloud_sfc_t', 'trends_dCEcloud_sfc_t', 'trendyr_dCEcloud_sfc_t', ...
    'trendm_dCEcloud_toa_t', 'trends_dCEcloud_toa_t', 'trendyr_dCEcloud_toa_t', ...
    'trendm_dCEcloud_sfc_c', 'trends_dCEcloud_sfc_c', 'trendyr_dCEcloud_sfc_c', ...
    'trendm_dCEcloud_toa_c', 'trends_dCEcloud_toa_c', 'trendyr_dCEcloud_toa_c');


end

t = toc; disp(t)

function inputData = transfCordi(lon_ori, lat_ori, time_ori, inputData, lon_aft, lat_aft, time_aft)
    inputData(2:end + 1, :, :) = inputData;
    inputData(1, :, :) = inputData(end, :, :);
    inputData = autoRegrid3(lon_ori, lat_ori, time_ori, inputData, lon_aft, lat_aft, time_aft);
end
