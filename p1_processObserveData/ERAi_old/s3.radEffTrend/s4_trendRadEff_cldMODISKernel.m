% cal 200207_201706 cloud Effect trend by MODIS data using kernel method
% raw data time series: 200207_201706
%
clc; clear; tic;
%modify path first
inputpath = '/data1/liuyincheng/Observe-rawdata/';
outputpath1 = '/data1/liuyincheng/Observe-process/200207-201706/MODIS/radEffect/';
outputpath2 = '/data1/liuyincheng/Observe-process/200207-201706/MODIS/radEffect_trend/';
auto_mkdir(outputpath1)
auto_mkdir(outputpath2)
lon_f = 0:2.5:357.5; nlonf = length(lon_f);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f);

p_1 = 2; % different time series, 1mean 2000 - 03 to 2018 - 02(18 * 12). 2 mean 200207 - 201706(15 * 12)
[readme, level, tLin, vars] = obsParameters('ERA5');
%------------------------ MODIS data ------------------------
load([inputpath, 'surface_flux_MODIS_kernel_estimation_200207_201908.mat'])% lat_m, lon_m, lwcf, swcf,, lat1,lon1,
% lat_m=89.5:-1:-89.5;lon_m=.5:1:359.5;
lat_m = double(lat); lon_m = double(lon);
lon_m(2:end + 1) = lon_m; lon_m(1) = -0.5;
Rsfc_cloud = swcf + lwcf;
Rsfc_cloud(:, 2:end + 1, :) = Rsfc_cloud; Rsfc_cloud(:, 1, :) = Rsfc_cloud(:, end, :);
Rsfc_cloud = permute(Rsfc_cloud, [2 1 3]);
% choose 200207_201706 time
Rsfc_cloud = Rsfc_cloud(:, :, 1:tLin.inter{p_1});

t0 = datetime(2002, 7, 1);
t = t0 + calmonths(0:179);
time = datenum(t);

% Regrid
Rsfc_cloud = autoRegrid3(lon_m, lat_m, time, Rsfc_cloud, lon_f, lat_f, time);

%% Deseasonalize
startMonth = tLin.startMonth{2};
[dRsfc_cloudMOD, Clim_Rsfc_cloudMOD] = monthlyAnomaly3D(nlonf, nlatf, time, Rsfc_cloud, startMonth);


%% cal trend
[trendm_dRsfc_cloudMOD, trends_dRsfc_cloudMOD, trendyr_dRsfc_cloudMOD, ~, ~] = autoCalTrend(dRsfc_cloudMOD, nlonf, nlatf, time, startMonth);


%%  ------------------------ erai data ------------------------
varLevel1{1} = 'sfc'; varLevel1{2} = 'toa';
varLevel2{1} = 'all'; varLevel2{2} = 'clr';
varLevel3{1} = 'ERAi'; varLevel3{2} = 'CERES';
realRad_source = 'ERAi';

%modify path first
inputpath_erai = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, '/');

%% ------------------------Read date------------------------
% radiative effect(144,72,tLin.inter{p_1},2,2)sfc/toa*cld/clr
s_albEff = zeros(144, 72, tLin.inter{p_1}, 2, 2); l_taEff = s_albEff; l_tsEff = s_albEff; l_wvEff = s_albEff; s_wvEff = s_albEff; totalEff = s_albEff;
l_rad = s_albEff; s_rad = s_albEff;

for ii = 1:2% 1.sfc, 2.toa
    var1 = char(varLevel1{ii});

    for jj = 1:2% 1.all 2.clr
        var2 = char(varLevel2{jj});
        eff_path = strcat(inputpath_erai, varLevel3{1}, '/radEffect/dradEffect_', var1, '_', var2, '.nc');
        realRad_path = strcat(inputpath_erai, realRad_source, '/anomaly/drad.mat');
        % grobal vars
        lon_f = ncread(eff_path, 'longitude'); % 144x72
        lat_f = ncread(eff_path, 'latitude');
        nlon = length(lon_f); nlat = length(lat_f);
        time = ncread(eff_path, 'time'); ntime = length(time);
        % read effects
        s_albEff(:, :, :, ii, jj) = ncread(eff_path, 'albRadEff');
        l_taEff(:, :, :, ii, jj) = ncread(eff_path, 'taRadEff');
        l_tsEff(:, :, :, ii, jj) = ncread(eff_path, 'tsRadEff');
        l_wvEff(:, :, :, ii, jj) = ncread(eff_path, 'wvlwRadEff');
        s_wvEff(:, :, :, ii, jj) = ncread(eff_path, 'wvswRadEff');
        totalEff(:, :, :, ii, jj) = ncread(eff_path, 'totalRadEff'); % alb + ta + ts + wv
        % read real rad Anomly
        l_rad(:, :, :, ii, jj) = cell2mat(struct2cell(load(realRad_path, ['dR_lw_', var1, '_', var2]))); %net thermal radiation
        s_rad(:, :, :, ii, jj) = cell2mat(struct2cell(load(realRad_path, ['dR_sw_', var1, '_', var2]))); %net thermal radiation
    end

end

% cal allsky vars radEffect(only need allsky)
% vars:  q, alb, ts, ta
husEff = l_wvEff + s_wvEff;
dR_rad = l_rad + s_rad;
dR_husSfcAll = squeeze(squeeze(husEff(:, :, :, 1, 1))); dR_husToaAll = squeeze(squeeze(husEff(:, :, :, 2, 1))); % q
dR_albSfcAll = squeeze(squeeze(s_albEff(:, :, :, 1, 1))); dR_albToaAll = squeeze(squeeze(s_albEff(:, :, :, 2, 1))); %albedo
dR_tsSfcAll = squeeze(squeeze(l_tsEff(:, :, :, 1, 1))); dR_tsToaAll = squeeze(squeeze(l_tsEff(:, :, :, 2, 1))); % ts
dR_taSfcAll = squeeze(squeeze(l_taEff(:, :, :, 1, 1))); dR_taToaAll = squeeze(squeeze(l_taEff(:, :, :, 2, 1))); % ta
dR_totalSfcAll = squeeze(squeeze(totalEff(:, :, :, 1, 1))); dR_totalToaAll = squeeze(squeeze(totalEff(:, :, :, 2, 1))); % ta
dR_realRadSfcAll = squeeze(squeeze(dR_rad(:, :, :, 1, 1))); dR_realRadToaAll = squeeze(squeeze(dR_rad(:, :, :, 2, 1))); % ta
dR_all = squeeze(dR_rad(:, :, :, :, 1));
dR_clr = squeeze(dR_rad(:, :, :, :, 2));
dRtotalEff_all = squeeze(totalEff(:, :, :, :, 1)); % dRfeedback_all
dRtotalEff_clr = squeeze(totalEff(:, :, :, :, 2)); % dRfeedback_clr

% cal residual effect(use modis result) % 1.sfc, 2.toa
% only sfc can be caled
dR_residualMOD_sfc = squeeze(dR_all(:, :, :, 1) - dRtotalEff_all(:, :, :, 1)) - dRsfc_cloudMOD;
dRheat_varsCal_sfc_Resi = dR_husSfcAll + dR_albSfcAll + dR_taSfcAll + dRsfc_cloudMOD + dR_residualMOD_sfc;
dRheat_varsCal_sfc_noResi = dR_husSfcAll + dR_albSfcAll + dR_taSfcAll + dRsfc_cloudMOD;
testVar_equal_Resi=squeeze(dR_all(:, :, :, 1))-dR_tsSfcAll;% 验证这一项和上一项是否相等

%% Save to mat
save([outputpath1, 'dcloud_MOD.mat'], 'lon_f', 'lat_f', 'time', 'readme'...
    , 'dRsfc_cloudMOD', 'Clim_Rsfc_cloudMOD'...
    ,'dR_residualMOD_sfc','dRheat_varsCal_sfc_Resi','dRheat_varsCal_sfc_noResi','testVar_equal_Resi')

load('/data1/liuyincheng/Observe-process/200207-201706/ERAi/radEffect/dcldRadEffect_ERAi.mat')% 'lon_f', 'lat_f', 'time', ...'dR_cloud_sfc', 'dR_cloud_toa','dR_residual_sfc','dR_residual_toa'
dRheat_varsCal_sfc_ResiERAi = dR_husSfcAll + dR_albSfcAll + dR_taSfcAll + dRsfc_cloudMOD + dR_residual_sfc;

%% cal trend
[trendm_dRsfc_residualMOD, trends_dRsfc_residualMOD, trendyr_dRsfc_residualMOD, ~, ~] = autoCalTrend(dR_residualMOD_sfc, nlon, nlat, time, startMonth);
[trendm_dRsfcHeat_varsCal_Resi, trends_dRsfcHeat_varsCal_Resi, trendyr_dRsfcHeat_varsCal_Resi, ~, ~] = autoCalTrend(dRheat_varsCal_sfc_Resi, nlon, nlat, time, startMonth);
[trendm_dRsfcHeat_varsCal_ResiERAi, trends_dRsfcHeat_varsCal_ResiERAi, trendyr_dRsfcHeat_varsCal_ResiERAi, ~, ~] = autoCalTrend(dRheat_varsCal_sfc_ResiERAi, nlon, nlat, time, startMonth);
[trendm_dRsfcHeat_varsCal_noResi, trends_dRsfcHeat_varsCal_noResi, trendyr_dRsfcHeat_varsCal_noResi, ~, ~] = autoCalTrend(dRheat_varsCal_sfc_noResi, nlon, nlat, time, startMonth);
[trendm_testVar_equal_Resi, trends_testVar_equal_Resi, trendyr_testVar_equal_Resi, ~, ~] = autoCalTrend(testVar_equal_Resi, nlon, nlat, time, startMonth);

% save as mat file
save([outputpath2, 'trend_dcloud_MOD.mat'], 'lon_f', 'lat_f', 'time', 'readme', ...
    'trendm_dRsfc_cloudMOD', 'trends_dRsfc_cloudMOD', 'trendyr_dRsfc_cloudMOD', ...
    'trendm_dRsfc_residualMOD', 'trends_dRsfc_residualMOD', 'trendyr_dRsfc_residualMOD', ...
    'trendm_dRsfcHeat_varsCal_Resi', 'trends_dRsfcHeat_varsCal_Resi', 'trendyr_dRsfcHeat_varsCal_Resi', ...
    'trendm_dRsfcHeat_varsCal_ResiERAi', 'trends_dRsfcHeat_varsCal_ResiERAi', 'trendyr_dRsfcHeat_varsCal_ResiERAi', ...
    'trendm_dRsfcHeat_varsCal_noResi', 'trends_dRsfcHeat_varsCal_noResi', 'trendyr_dRsfcHeat_varsCal_noResi', ...
    'trendm_testVar_equal_Resi', 'trends_testVar_equal_Resi', 'trendyr_testVar_equal_Resi')
t = toc; disp(t)
