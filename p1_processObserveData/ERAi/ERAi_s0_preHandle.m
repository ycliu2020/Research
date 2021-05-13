%%---------------------------------------------------------
% Author       : LYC
% Date         : 2021-04-22 14:43:30
% LastEditors  : Please set LastEditors
% Description  : 
% FilePath     : /code/p1_processObserveData/ERAi/ERAi_s0_preHandle.m
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
outputPath_rawdata= '/data2/liuyincheng/Observe-process/rawdata/ERAi/';
outputPath_rawdataRegrid= '/data2/liuyincheng/Observe-process/rawdata_regrid/ERAi/';
auto_mkdir(outputPath_rawdata);
auto_mkdir(outputPath_rawdataRegrid);
% input path
obsPath = '/data2/liuyincheng/Observe-rawdata/ERA/ERAi';
metoVarsPath = fullfile(obsPath, 'meto_vars/');
sfc_radVarsPath = fullfile(obsPath, 'rad_vars/SFC/');
toa_radVarsPath = fullfile(obsPath, 'rad_vars/TOA/');
[readme, level, tLin, vars] = obsParameters('ERAi');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read raw data from 2000-2018 for futher process to 200003-201802,,,,
% 1.1 连续时间的变量
alb_Path = fullfile(metoVarsPath, 'alb_1979-2018.nc');
ps_Path = fullfile(metoVarsPath, 'ps_1979-2018.nc');
ts_Path = fullfile(metoVarsPath, 'ts_1979-2018.nc');
% read time
% locate first year
formatOut = 'yyyy-mm';
layerSfc_dir = dir(alb_Path);
startT = cdftime2loc(layerSfc_dir, formatOut, '1979-01');
% read time and transfor mat time
yearTotal=2018-1979+1;
startLoc = startT; count = yearTotal * 12; stride = 1;
time_raw = cmipTimeRead(alb_Path, startLoc, count, stride); %time.date, Units, Calendar, length
ntime_raw = length(time_raw.date);
time1 = 1:ntime_raw;
% read rawginal lat and lon
lat_raw = double(ncread(alb_Path, 'latitude'));
lon_raw = double(ncread(alb_Path, 'longitude'));
% read level
TQ_Path = strcat(metoVarsPath, 'tq', num2str((0 + 1979), '%04i'), '.nc');
plev_raw = ncread(TQ_Path, 'level');
plev_raw = double(plev_raw);
plev_raw = plev_raw(end:-1:1);
plev_k = ncread('/data1/liuyincheng/y_kernels/YiH/kernels_YiH/toa/dp.nc', 'player'); % player 最大1000符合条件pay attention to plev_k's range must smaller than plev_raw's
nplevk = length(plev_k);

% read vars
startLoc = [1, 1, startT]; count = [inf, inf, yearTotal * 12]; stride = [1, 1, 1];
alb_raw = squeeze(ncread(alb_Path, 'fal', startLoc, count, stride));
ps_raw = squeeze(ncread(ps_Path, 'sp', startLoc, count, stride));
ts_raw = squeeze(ncread(ts_Path, 'skt', startLoc, count, stride));
% Regrid to 2.5*2.5
% 3D meto vars
ps_regrid = autoRegrid3(lon_raw, lat_raw, time1, ps_raw, lon_k, lat_k, time1);
alb_regrid = autoRegrid3(lon_raw, lat_raw, time1, alb_raw, lon_k, lat_k, time1);
ts_regrid = autoRegrid3(lon_raw, lat_raw, time1, ts_raw, lon_k, lat_k, time1);
% save raw data 
save([outputPath_rawdata, 'global_vars.mat'], 'lon_raw', 'lat_raw', 'time_raw', 'plev_raw')
save([outputPath_rawdata, 'ps_raw.mat'], 'ps_raw','-v7.3')
save([outputPath_rawdata, 'alb_raw.mat'], 'alb_raw','-v7.3')
save([outputPath_rawdata, 'ts_raw.mat'], 'ts_raw','-v7.3')
clear alb_raw ps_raw ts_raw
% save regrided data
save([outputPath_rawdataRegrid, 'global_vars.mat'], 'lon_k', 'lat_k', 'time_raw', 'plev_k')
save([outputPath_rawdataRegrid, 'ps_regrid.mat'], 'ps_regrid')
save([outputPath_rawdataRegrid, 'alb_regrid.mat'], 'alb_regrid')
save([outputPath_rawdataRegrid, 'ts_regrid.mat'], 'ts_regrid')
clear alb_regrid ps_regrid ts_regrid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1.2 按时间划分文件的变量

%%%%%%%%%%%%%%%%%%%%%%%%%%% meto data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read data
ta_raw=zeros(length(lon_raw),length(lat_raw),length(plev_raw),yearTotal*12);
hus_raw=ta_raw;
for yearNum = 0:yearTotal-1% (read raw range: 2000-2018)
    TQ_Path = strcat(metoVarsPath, 'tq', num2str((yearNum + 1979), '%04i'), '.nc');
    ta_raw(:, :, :, yearNum * 12 + 1:(yearNum + 1) * 12) = ncread(TQ_Path, 't');
    hus_raw(:, :, :, yearNum * 12 + 1:(yearNum + 1) * 12) = ncread(TQ_Path, 'q');

end
% transform plev direction
hus_raw = hus_raw(:, :, end:-1:1, :);
ta_raw = ta_raw(:, :, end:-1:1, :);
% save raw data 
save([outputPath_rawdata, 'global_vars.mat'], 'lon_raw', 'lat_raw', 'time_raw', 'plev_raw')
save([outputPath_rawdata, 'hus_raw.mat'], 'hus_raw','-v7.3')
save([outputPath_rawdata, 'ta_raw.mat'], 'ta_raw','-v7.3')

%  Regrid to 2.5*2.5, method 1
% 4D meto vars
time2=1:12;
for yearNum = 0:yearTotal-1 % (read raw range: 2000-2018)
    hus_regrid(:, :, :, yearNum * 12 + 1:(yearNum + 1) * 12) = autoRegrid4(lon_raw, lat_raw, plev_raw, time2, hus_raw(:, :, :, yearNum * 12 + 1:(yearNum + 1) * 12), lon_k, lat_k, plev_k, time2);
    ta_regrid(:, :, :, yearNum * 12 + 1:(yearNum + 1) * 12) = autoRegrid4(lon_raw, lat_raw, plev_raw, time2, ta_raw(:, :, :, yearNum * 12 + 1:(yearNum + 1) * 12), lon_k, lat_k, plev_k, time2);
end
clear hus_raw ta_raw

% save regrided data
save([outputPath_rawdataRegrid, 'global_vars.mat'], 'lon_k', 'lat_k', 'time_raw', 'plev_k')
save([outputPath_rawdataRegrid, 'hus_regrid.mat'], 'hus_regrid')
save([outputPath_rawdataRegrid, 'ta_regrid.mat'], 'ta_regrid')
clear hus_regrid ta_regrid

%%%%%%%%%%%%%%%%%%%%%%%%%%% rad data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rad_sfc_raw=zeros(length(lon_raw),length(lat_raw),yearTotal*12,length(vars.sfcRad));
rad_toa_raw=rad_sfc_raw;
for yearNum = 0:yearTotal-1% (read raw range: 2000-2018)
    sfcRad_Path = strcat(sfc_radVarsPath, 'rad_sfc', num2str((yearNum + 1979), '%04i'), '.nc');
    toaRad_Path = strcat(toa_radVarsPath, 'rad_toa', num2str((yearNum + 1979), '%04i'), '.nc');
    for radVarNum = 1:length(vars.sfcRad)
        rad_sfc_raw(:, :, yearNum * 12 + 1:(yearNum + 1) * 12, radVarNum) = ncread(sfcRad_Path, vars.sfcRad{radVarNum});
    end

    for radVarNum = 1:length(vars.toaRad)
        rad_toa_raw(:, :, yearNum * 12 + 1:(yearNum + 1) * 12, radVarNum) = ncread(toaRad_Path, vars.toaRad{radVarNum});
    end

end

rad_sfc_regrid = zeros(nlonk, nlatk, ntime_raw, length(vars.sfcRad));
rad_toa_regrid = zeros(nlonk, nlatk, ntime_raw, length(vars.toaRad));

for radVarNum = 1:length(vars.toaRad)
    rad_toa_regrid(:, :, :, radVarNum) = autoRegrid3(lon_raw, lat_raw, time1, rad_toa_raw(:, :, :, radVarNum), lon_k, lat_k, time1);
end

for radVarNum = 1:length(vars.sfcRad)
    rad_sfc_regrid(:, :, :, radVarNum) = autoRegrid3(lon_raw, lat_raw, time1, rad_sfc_raw(:, :, :, radVarNum), lon_k, lat_k, time1);
end

% save raw data 
save([outputPath_rawdata, 'rad_sfc_raw.mat'], 'rad_sfc_raw','-v7.3')
save([outputPath_rawdata, 'rad_toa_raw.mat'], 'rad_toa_raw','-v7.3')
clear rad_sfc_raw rad_toa_raw 

% save regrided data
save([outputPath_rawdataRegrid, 'rad_sfc_regrid.mat'], 'rad_sfc_regrid')
save([outputPath_rawdataRegrid, 'rad_toa_regrid.mat'], 'rad_toa_regrid')
clear rad_sfc_regrid rad_toa_regrid 

t = toc; disp(t)
