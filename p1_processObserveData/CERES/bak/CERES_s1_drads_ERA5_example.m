%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-07-08 10:22:15
% LastEditTime : 2021-05-02 15:12:19
% LastEditors  : Please set LastEditors
% Description  :
% FilePath     : /code/p1_processObserveData/CERES/CERES_s1_drads_ERA5_example.m
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

% modify path first
ceresDataPath = '/data2/liuyincheng/Observe-rawdata/CERES_EBAF/4.1/';
sfcDataPath = strcat(ceresDataPath, 'CERES_EBAF-SFC_Ed4.1_Subset_200003-202011.nc');
toaDataPath = strcat(ceresDataPath, 'CERES_EBAF-TOA_Ed4.1_Subset_200003-202011.nc');

% transfor mat time
time = double(ncread(sfcDataPath, 'time')); % days since 2000-03-01 00:00:00
timeUnits = ncreadatt(sfcDataPath, 'time', 'units');
timeCalendar = 'gregorian';
time = cdfdate2num(timeUnits, timeCalendar, time);
ntime = length(time);

%% read CERES raw data
%% SFC
% All-sky
dR_tot_sfc_all = ncread(sfcDataPath, 'sfc_net_tot_all_mon'); % = sfc_net_lw_all_mon + sfc_net_sw_all_mon
dR_cre_sfc_all_t = ncread(sfcDataPath, 'sfc_cre_net_tot_mon'); % Cloud radiative effect(CRF: sfc_net_tot_all_mon-sfc_net_tot_clr_mon_t)
dR_lwDown_sfc_all = ncread(sfcDataPath, 'sfc_lw_down_all_mon'); % LW downwards
dR_swNet_sfc_all = ncread(sfcDataPath, 'sfc_net_sw_all_mon'); % surface_net_downward_shortwave_flux
% Clear-sky
dR_tot_sfc_clr_t = ncread(sfcDataPath, 'sfc_net_tot_clr_t_mon'); % transfer model result
dR_tot_sfc_clr_c = ncread(sfcDataPath, 'sfc_net_tot_clr_c_mon'); % obs for cloud-free areas result

%% TOA
% allsky
dR_tot_toa_all = ncread(toaDataPath, 'toa_net_all_mon'); % = - toa_sw_all_mon - toa_lw_all_mon + solar_mon
dR_cre_toa_all_t = ncread(toaDataPath, 'toa_cre_net_mon'); % Cloud radiative effect
% clear sky
dR_tot_toa_clr_t = ncread(toaDataPath, 'toa_net_clr_t_mon');
dR_tot_toa_clr_c = ncread(toaDataPath, 'toa_net_clr_c_mon');

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

% necessary parameters
[readme, level, tLin, vars] = obsParameters('CERES');
obsPath='/data2/liuyincheng/Observe-process';
%% different time series
for p_1 = 6:6 % p_1 mean different time series, 1 mean 2000-03 to 2018-02(18*12). 2 mean 200207-201706(15*12). 6 mean '200003-202011'
    anomPath = fullfile(obsPath, tLin.time{p_1}, 'CERES', level.standVarPath{2}); % /anomaly
    anomPath_era5 = fullfile(obsPath, tLin.time{p_1}, 'ERA5', level.standVarPath{2}); % /anomaly
    anomPath_erai = fullfile(obsPath, tLin.time{p_1}, 'ERAi', level.standVarPath{2}); % /anomaly
    anomTrendPath = fullfile(obsPath, tLin.time{p_1}, 'CERES', level.standVarPath{3}); % /anomaly_trend
    radEffectPath_era5 = fullfile(obsPath, tLin.time{p_1}, 'ERA5', level.standVarPath{5}); %radEffect
    radEffectPath_erai = fullfile(obsPath, tLin.time{p_1}, 'ERAi', level.standVarPath{5}); %radEffect
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
    %% Part2: cal cloud Effect and trend
    %% ERA5 ,  _t
    % load  ERA5 ('final column: 1sfc cld,2sfc clr,3toa cld,4toa clr') 
    load([radEffectPath_era5, 'dradEffect_union.mat']); %'dR_hus', 'dR_alb', 'dR_ts', 'dR_ta', 'dR_nonCloud', 'readme_radEfect'
    dCE_CRF_t(:, :, :, 1) = dR_cre_sfc_all_t;
    dCE_CRF_t(:, :, :, 2) = dR_cre_toa_all_t;
    dvarsFeedback = dR_nonCloud(:, :, :, [2 4]) - dR_nonCloud(:, :, :, [1 3]);
    dRCE_cloud = dCE_CRF_t + dvarsFeedback;
    dRCE_clr(:, :, :, 1) = dR_tot_sfc_clr_t;
    dRCE_clr(:, :, :, 2) = dR_tot_toa_clr_t;
    dRCE_residual = dRCE_clr - dR_nonCloud(:, :, :, [2 4]);
    dR_main = dR_hus + dR_alb + dR_ta; % dR exclude ts
    dR_mainEffect_sfc = squeeze(dRCE_cloud(:, :, :, 1)) + squeeze(dR_main(:, :, :, 1));
    dR_mainEffect_toa = squeeze(dRCE_cloud(:, :, :, 2)) + squeeze(dR_main(:, :, :, 3));
    % save as one file
    dCE_era5_t.Rcloud_sfc = squeeze(dRCE_cloud(:, :, :, 1));
    dCE_era5_t.Rcloud_toa = squeeze(dRCE_cloud(:, :, :, 2));
    dCE_era5_t.residual_sfc = squeeze(dRCE_residual(:, :, :, 1));
    dCE_era5_t.residual_toa = squeeze(dRCE_residual(:, :, :, 2));
    dCE_era5_t.mainEffect_sfc = dR_mainEffect_sfc;
    dCE_era5_t.mainEffect_toa = dR_mainEffect_toa;
    % _c
    dCE_CRF_c(:, :, :, 1) = dR_cre_sfc_all_c;
    dCE_CRF_c(:, :, :, 2) = dR_cre_toa_all_c;
    dvarsFeedback = dR_nonCloud(:, :, :, [2 4]) - dR_nonCloud(:, :, :, [1 3]);
    dRCE_cloud = dCE_CRF_c + dvarsFeedback;
    dRCE_clr(:, :, :, 1) = dR_tot_sfc_clr_c;
    dRCE_clr(:, :, :, 2) = dR_tot_toa_clr_c;
    dRCE_residual = dRCE_clr - dR_nonCloud(:, :, :, [2 4]);
    dR_main = dR_hus + dR_alb + dR_ta; % dR exclude ts
    dR_mainEffect_sfc = squeeze(dRCE_cloud(:, :, :, 1)) + squeeze(dR_main(:, :, :, 1));
    dR_mainEffect_toa = squeeze(dRCE_cloud(:, :, :, 2)) + squeeze(dR_main(:, :, :, 3));
    % save as one file
    dCE_era5_c.Rcloud_sfc = squeeze(dRCE_cloud(:, :, :, 1));
    dCE_era5_c.Rcloud_toa = squeeze(dRCE_cloud(:, :, :, 2));
    dCE_era5_c.residual_sfc = squeeze(dRCE_residual(:, :, :, 1));
    dCE_era5_c.residual_toa = squeeze(dRCE_residual(:, :, :, 2));
    dCE_era5_c.mainEffect_sfc = dR_mainEffect_sfc;
    dCE_era5_c.mainEffect_toa = dR_mainEffect_toa;
    save([anomPath, 'dCEcloud_era5.mat'], 'dCE_era5_t', 'dCE_era5_c')

    % %% ERAi ,  _t
    % % load  ( 1sfc cld,2sfc clr,3toa cld,4toa clr)
    % load([radEffectPath_erai, 'dradEffect_union.mat']); %'dR_hus', 'dR_alb', 'dR_ts', 'dR_ta', 'dR_nonCloud', 'readme_radEfect'
    % dCE_CRF_t(:, :, :, 1) = dR_cre_sfc_all_t;
    % dCE_CRF_t(:, :, :, 2) = dR_cre_toa_all_t;
    % dvarsFeedback = dR_nonCloud(:, :, :, [2 4]) - dR_nonCloud(:, :, :, [1 3]);
    % dRCE_cloud = dCE_CRF_t + dvarsFeedback;
    % dRCE_clr(:, :, :, 1) = dR_tot_sfc_clr_t;
    % dRCE_clr(:, :, :, 2) = dR_tot_toa_clr_t;
    % dRCE_residual = dRCE_clr - dR_nonCloud(:, :, :, [2 4]);
    % dR_main = dR_hus + dR_alb + dR_ta; % dR exclude ts
    % dR_mainEffect_sfc = squeeze(dRCE_cloud(:, :, :, 1)) + squeeze(dR_main(:, :, :, 1));
    % dR_mainEffect_toa = squeeze(dRCE_cloud(:, :, :, 2)) + squeeze(dR_main(:, :, :, 3));
    % %save as one file
    % dCE_erai_t.Rcloud_sfc = squeeze(dRCE_cloud(:, :, :, 1));
    % dCE_erai_t.Rcloud_toa = squeeze(dRCE_cloud(:, :, :, 2));
    % dCE_erai_t.residual_sfc = squeeze(dRCE_residual(:, :, :, 1));
    % dCE_erai_t.residual_toa = squeeze(dRCE_residual(:, :, :, 2));
    % dCE_erai_t.mainEffect_sfc = dR_mainEffect_sfc;
    % dCE_erai_t.mainEffect_toa = dR_mainEffect_toa;
    % % _c
    % dCE_CRF_c(:, :, :, 1) = dR_cre_sfc_all_c;
    % dCE_CRF_c(:, :, :, 2) = dR_cre_toa_all_c;
    % dvarsFeedback = dR_nonCloud(:, :, :, [2 4]) - dR_nonCloud(:, :, :, [1 3]);
    % dRCE_cloud = dCE_CRF_c + dvarsFeedback;
    % dRCE_clr(:, :, :, 1) = dR_tot_sfc_clr_c;
    % dRCE_clr(:, :, :, 2) = dR_tot_toa_clr_c;
    % dRCE_residual = dRCE_clr - dR_nonCloud(:, :, :, [2 4]);
    % dR_main = dR_hus + dR_alb + dR_ta; % dR exclude ts
    % dR_mainEffect_sfc = squeeze(dRCE_cloud(:, :, :, 1)) + squeeze(dR_main(:, :, :, 1));
    % dR_mainEffect_toa = squeeze(dRCE_cloud(:, :, :, 2)) + squeeze(dR_main(:, :, :, 3));
    % %save as one file
    % dCE_erai_c.Rcloud_sfc = squeeze(dRCE_cloud(:, :, :, 1));
    % dCE_erai_c.Rcloud_toa = squeeze(dRCE_cloud(:, :, :, 2));
    % dCE_erai_c.residual_sfc = squeeze(dRCE_residual(:, :, :, 1));
    % dCE_erai_c.residual_toa = squeeze(dRCE_residual(:, :, :, 2));
    % dCE_erai_c.mainEffect_sfc = dR_mainEffect_sfc;
    % dCE_erai_c.mainEffect_toa = dR_mainEffect_toa;
    % save([anomPath, 'dCEcloud_erai.mat'], 'dCE_erai_t', 'dCE_erai_c')

    
    clear dCE_CRF_t dCE_CRF_c dRCE_clr dRCE_residual dRCE_cloud

end

t = toc; disp(t)

function inputData = transfCordi(lon_ori, lat_ori, time_ori, inputData, lon_aft, lat_aft, time_aft)
    inputData(2:end + 1, :, :) = inputData;
    inputData(1, :, :) = inputData(end, :, :);
    inputData = autoRegrid3(lon_ori, lat_ori, time_ori, inputData, lon_aft, lat_aft, time_aft);
end

