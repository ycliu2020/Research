%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-07-02 15:23:02
% LastEditTime : 2020-07-02 20:31:03
% LastEditors  : LYC
% Description  : process ERA5 data to anomaly (includ meto vars and rad)
% FilePath     : /code/p1_processObserveData/s1.obersveDataProcess/ERA5/s1_dvars_ERA5.m
%
%%---------------------------------------------------------
clc; clear; tic;
% constant
sigma = 5.67e-8; % Stefan-Boltzmann constant: Wm-2K-4
kernels_path = '/data1/liuyincheng/y_kernels/kernels_YiH/toa/dp.nc';
lonf = 0:2.5:357.5; nlonf = length(lonf);
latf = 90:-2.5:-90; nlatf = length(latf);
%modify path first
obsPath = '/data1/liuyincheng/Observe-rawdata/ERA/ERA5';
metoVarsPath = fullfile(obsPath, 'meto_vars/');
sfc_radVarsPath = fullfile(obsPath, 'meto_vars/SFC/');
toa_radVarsPath = fullfile(obsPath, 'meto_vars/TOA/');

%%Part1 read raw data.
% 1.1 零碎时间的变量
for yearNum = 0:18% (read raw range: 2000-2018)
    % meto data
    TQ_Path = strcat(metoVarsPath, 'TQ', num2str((yearNum + 2000), '%04i'), '.nc');
    t0(:, :, :, yearNum * 12 + 1:(yearNum + 1) * 12) = ncread(TQ_Path, 't');
    q0(:, :, :, yearNum * 12 + 1:(yearNum + 1) * 12) = ncread(TQ_Path, 'q');
    time0(yearNum * 12 + 1:(yearNum + 1) * 12) = ncread(TQ_Path, 'time');

    % rad data
    sfcRad_Path = strcat(sfc_radVarsPath, 'rad_sfc', num2str((yearNum + 2000), '%04i'), '.nc');
    toaRad_Path = strcat(toa_radVarsPath, 'rad_toa', num2str((yearNum + 2000), '%04i'), '.nc');

    for radVarNum = 1:length(vars.sfcRad)
        rad_sfc0(:, :, :, yearNum * 12 + 1:(yearNum + 1) * 12, radVarNum) = ncread(sfcRad_Path, vars.sfcRad{radVarNum});
    end

    for radVarNum = 1:length(vars.toaRad)
        rad_toa0(:, :, :, yearNum * 12 + 1:(yearNum + 1) * 12, radVarNum) = ncread(toaRad_Path, vars.toaRad{radVarNum});
    end

end

% transfor mat time
temp.name = strcat(sfc_radVarsPath, 'TQ', num2str((4 + 2000), '%04i'), '.nc');
timeUnits = ncreadatt(temp.name, 'time', 'units');
% timeCalendar = ncreadatt(temp.name, 'time', 'calendar');
timeCalendar = 'gregorian';
time0 = cdfdate2num(timeUnits, timeCalendar, time0);
% read level
plev = ncread(temp.name, 'level');
plev = double(plev);
plevf1 = ncread('/data1/liuyincheng/y_kernels/kernels_YiH/toa/dp.nc', 'player'); % player 最大1000符合条件pay attention to plevf's range must smaller than plev's
plevf = plevf1(end:-1:1); % 调整顺序

% 1.2 连续时间的变量
alb_Path = fullfile(metoVarsPath, 'alb_1979-2019.nc');
ps_Path = fullfile(metoVarsPath, 'ps_1979-2019.nc');
ts_Path = fullfile(metoVarsPath, 'ts_1979-2019.nc');


for p_1 = 1:1
    [readme, level, tLin, vars] = observeParameters(p_1);

end

t = toc; disp(t)
