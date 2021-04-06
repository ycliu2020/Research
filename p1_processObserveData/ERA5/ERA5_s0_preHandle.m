%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-07-02 15:23:02
% LastEditTime : 2021-04-06 11:13:21
% LastEditors  : Please set LastEditors
% Description  : 预先将所有的ERA5文件读取进一个大文件插值并保存, 这样下次就不用每次都读取了
%                只需执行一次    
% FilePath     : /code/p1_processObserveData/ERA5/ERA5_s0_preHandle.m
%
%%---------------------------------------------------------
clc; clear; tic;
% constant
sigma = 5.67e-8; % Stefan-Boltzmann constant: Wm-2K-4
kernels_path = '/data1/liuyincheng/y_kernels/YiH/kernels_YiH/toa/dp.nc';
lon_k = 0:2.5:357.5; nlonk = length(lon_k); % kernel lat lon
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);
var_state = {'d', 'clim_', 'trendm_d', 'trends_d', 'trendyr_d'};
%% modify path first
% ouput path
outputPath_rawdata= '/data1/liuyincheng/Observe-process/rawdata/ERA5/';
outputPath_rawdataRegrid= '/data1/liuyincheng/Observe-process/rawdata_regrid/ERA5/';
auto_mkdir(outputPath_rawdata);
auto_mkdir(outputPath_rawdataRegrid);
% input path
obsPath = '/data1/liuyincheng/Observe-rawdata/ERA/ERA5';
metoVarsPath = fullfile(obsPath, 'meto_vars/');
sfc_radVarsPath = fullfile(obsPath, 'rad_vars/SFC/');
toa_radVarsPath = fullfile(obsPath, 'rad_vars/TOA/');
[readme, level, tLin, vars] = obsParameters('ERA5');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read raw data from 2000-2018 for futher process to 200003-201802,,,,
% 1.1 连续时间的变量
alb_Path = fullfile(metoVarsPath, 'alb_1979-2019.nc');
ps_Path = fullfile(metoVarsPath, 'ps_1979-2019.nc');
ts_Path = fullfile(metoVarsPath, 'ts_1979-2019.nc');
% read time
% locate first year
formatOut = 'yyyy-mm';
layerSfc_dir = dir(alb_Path);
startT = cdftime2loc(layerSfc_dir, formatOut, '1979-01');
% read time and transfor mat time
yearTotal=2019-1979+1;
startLoc = startT; count = yearTotal * 12; stride = 1;
time_raw = cmipTimeRead(alb_Path, startLoc, count, stride); %time.date, Units, Calendar, length
ntime_raw = length(time_raw.date);
time1 = 1:ntime_raw;
% read rawginal lat and lon
lat_raw = double(ncread(alb_Path, 'latitude'));
lon_raw = double(ncread(alb_Path, 'longitude'));
% read vars
startLoc = [1, 1, 1, startT]; count = [inf, inf, 1, yearTotal * 12]; stride = [1, 1, 1, 1];
alb_raw = squeeze(ncread(alb_Path, 'fal', startLoc, count, stride));
ps_raw = squeeze(ncread(ps_Path, 'sp', startLoc, count, stride));
ts_raw = squeeze(ncread(ts_Path, 'skt', startLoc, count, stride));
% Regrid to 2.5*2.5
% 3D meto vars
ps_regrid = autoRegrid3(lon_raw, lat_raw, time1, ps_raw, lon_k, lat_k, time1);
alb_regrid = autoRegrid3(lon_raw, lat_raw, time1, alb_raw, lon_k, lat_k, time1);
ts_regrid = autoRegrid3(lon_raw, lat_raw, time1, ts_raw, lon_k, lat_k, time1);
% save raw data 
save([outputPath_rawdata, 'global_vars.mat'], 'lon_raw', 'lat_raw', 'time_raw')
save([outputPath_rawdata, 'ps_raw.mat'], 'ps_raw')
save([outputPath_rawdata, 'alb_raw.mat'], 'alb_raw')
save([outputPath_rawdata, 'ts_raw.mat'], 'ts_raw')
clear alb_raw ps_raw ts_raw
% save regrided data
save([outputPath_rawdataRegrid, 'global_vars.mat'], 'lon_f', 'lat_f', 'time_raw')
save([outputPath_rawdataRegrid, 'ps_regrid.mat'], 'ps_regrid')
save([outputPath_rawdataRegrid, 'alb_regrid.mat'], 'alb_regrid')
save([outputPath_rawdataRegrid, 'ts_regrid.mat'], 'ts_regrid')
clear alb_regrid ps_regrid ts_regrid

% 1.2 按时间划分文件的变量
for yearNum = 0:yearTotal-1% (read raw range: 2000-2018)
    % meto data
    TQ_Path = strcat(metoVarsPath, 'TQ', num2str((yearNum + 1979), '%04i'), '.nc');
    ta_raw(:, :, :, yearNum * 12 + 1:(yearNum + 1) * 12) = ncread(TQ_Path, 't');
    hus_raw(:, :, :, yearNum * 12 + 1:(yearNum + 1) * 12) = ncread(TQ_Path, 'q');

    % rad data({'tsr', 'ttr', 'tsrcs', 'ttrcs', 'slhf', 'sshf', 'ssr', 'str', 'strd', 'ssrd', 'ssrc', 'strc', 'strdc', 'ssrdc'};)
    sfcRad_Path = strcat(sfc_radVarsPath, 'rad_sfc', num2str((yearNum + 2000), '%04i'), '.nc');
    toaRad_Path = strcat(toa_radVarsPath, 'rad_toa', num2str((yearNum + 2000), '%04i'), '.nc');

    for radVarNum = 1:length(vars.sfcRad)
        rad_sfc_raw(:, :, yearNum * 12 + 1:(yearNum + 1) * 12, radVarNum) = ncread(sfcRad_Path, vars.sfcRad{radVarNum});
    end

    for radVarNum = 1:length(vars.toaRad)
        rad_toa_raw(:, :, yearNum * 12 + 1:(yearNum + 1) * 12, radVarNum) = ncread(toaRad_Path, vars.toaRad{radVarNum});
    end

end

% read level
plev_raw = ncread(TQ_Path, 'level');
plev_raw = double(plev_raw);
plev_raw = plev_raw(end:-1:1);
plev_k = ncread('/data1/liuyincheng/y_kernels/YiH/kernels_YiH/toa/dp.nc', 'player'); % player 最大1000符合条件pay attention to plev_k's range must smaller than plev_raw's
nplevk = length(plev_k);

% transform plev direction
hus_raw = hus_raw(:, :, end:-1:1, :);
ta_raw = ta_raw(:, :, end:-1:1, :);

% Regrid to 2.5*2.5
% 4D meto vars
hus_regrid = autoRegrid4(lon_raw, lat_raw, plev_raw, time1, hus_raw, lon_k, lat_k, plev_k, time1);
ta_regrid = autoRegrid4(lon_raw, lat_raw, plev_raw, time1, ta_raw, lon_k, lat_k, plev_k, time1);

% rad vars
rad_sfc_regrid = zeros(nlonk, nlatk, ntime_raw, length(vars.sfcRad));
rad_toa_regrid = zeros(nlonk, nlatk, ntime_raw, length(vars.toaRad));

for radVarNum = 1:length(vars.toaRad)
    rad_toa_regrid(:, :, :, radVarNum) = autoRegrid3(lon_raw, lat_raw, time1, rad_toa_raw(:, :, :, radVarNum), lon_k, lat_k, time1);
end

for radVarNum = 1:length(vars.sfcRad)
    rad_sfc_regrid(:, :, :, radVarNum) = autoRegrid3(lon_raw, lat_raw, time1, rad_sfc_raw(:, :, :, radVarNum), lon_k, lat_k, time1);
end

% save raw data 
save([outputPath_rawdata, 'global_vars.mat'], 'lon_raw', 'lat_raw', 'time_raw', 'plev_raw')
save([outputPath_rawdata, 'hus_raw.mat'], 'hus_raw')
save([outputPath_rawdata, 'ta_raw.mat'], 'ta_raw')
save([outputPath_rawdata, 'rad_sfc_raw.mat'], 'rad_sfc_raw')
save([outputPath_rawdata, 'rad_toa_raw.mat'], 'rad_toa_raw')
clear rad_sfc_raw rad_toa_raw hus_raw ta_raw

% save regrided data
save([outputPath_rawdataRegrid, 'global_vars.mat'], 'lon_f', 'lat_f', 'time_raw', 'plev_k')
save([outputPath_rawdataRegrid, 'hus_regrid.mat'], 'hus_regrid')
save([outputPath_rawdataRegrid, 'ta_regrid.mat'], 'ta_regrid')
save([outputPath_rawdataRegrid, 'rad_sfc_regrid.mat'], 'rad_sfc_regrid')
save([outputPath_rawdataRegrid, 'rad_toa_regrid.mat'], 'rad_toa_regrid')
clear rad_sfc_regrid rad_toa_regrid hus_regrid ta_regrid

t = toc; disp(t)
