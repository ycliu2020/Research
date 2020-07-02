%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-09 15:52:00
% LastEditTime : 2020-07-02 20:14:43
% LastEditors  : LYC
% Description  : 
% FilePath     : /code/p1_processObserveData/s1.obersveDataProcess/ERAi/s1_dvars_ERAi.m
%  
%%---------------------------------------------------------
%% Combine the ERAi climatological data from 2000-03 to 2018-02
clc; clear; tic;
%modify path first
inputpath1 = '/data1/xieyan/ERAi/Surface_ori/';
inputpath2 = '/data1/xieyan/ERAi/TQ_ori/';

%% read raw data. 228 mean raw date 19yrs time line
a0 = zeros(360, 181, 228); % Forecast albedo
ts0 = zeros(360, 181, 228); % Skin temperature
q0 = zeros(360, 181, 37, 228); % Specific humidity
t0 = zeros(360, 181, 37, 228); % Temperature
time0 = zeros(228, 1);
for ii = 0:18%(read raw range: 2000-2018)
    filename1 = strcat(inputpath1, 'sfc', num2str((ii + 2000), '%04i'), '.nc');
    filename2 = strcat(inputpath2, 'tq', num2str((ii + 2000), '%04i'), '.nc');

    temp1 = ncread(filename1, 'fal');
    temp2 = ncread(filename1, 'skt');
    temp3 = ncread(filename2, 'q');
    temp4 = ncread(filename2, 't');
    temp5 = ncread(filename2, 'time');
    % only read  0-12 is ok
    a0(:, :, ii * 12 + 1:(ii + 1) * 12) = temp1;
    ts0(:, :, ii * 12 + 1:(ii + 1) * 12) = temp2;
    q0(:, :, :, ii * 12 + 1:(ii + 1) * 12) = temp3;
    t0(:, :, :, ii * 12 + 1:(ii + 1) * 12) = temp4;
    time0(ii * 12 + 1:(ii + 1) * 12, 1) = temp5;
end
% transfor mat time
temp.name = strcat(inputpath1, 'sfc', num2str((4 + 2000), '%04i'), '.nc');
timeUnits = ncreadatt(temp.name, 'time', 'units');
% timeCalendar = ncreadatt(temp.name, 'time', 'calendar');
timeCalendar = 'gregorian';
time0 = cdfdate2num(timeUnits, timeCalendar, time0);

%% different time series
for p_1 = 1:2% 1mean 2000-03 to 2018-02(18*12). 2 mean 200207-201706(15*12)
    [readme, tLin] = observeParameters(p_1); % fuction% readme,  tLin,
    outpath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi/anomaly/');
    outpath1 = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi/rawdata_regrid/');
    auto_mkdir(outpath)
    auto_mkdir(outpath1)

    % find start and end month
    timeSpec = tLin.start{p_1};
    formatOut = 'yyyy-mm';
    timeStr = string(datestr(time0, formatOut));
    startT = find(timeStr == timeSpec);
    endT = startT + tLin.inter{p_1} - 1;
    
    % read specific time series
    alb = a0(:, :, startT:endT);
    ts = ts0(:, :, startT:endT);
    q = q0(:, :, :, startT:endT);
    ta = t0(:, :, :, startT:endT);
    time = time0(startT:endT, 1);
    ntime = length(time);

    %% Regrid to 2.5*2.5
    lon = 0:359; lat = 90:-1:-90;
    plev = ncread('/data1/xieyan/ERAi/TQ_ori/tq2017.nc', 'level');
    plev = double(plev);
    time2 = 1:ntime; %
    [Xlon, Ylat, Zplev, Ttime] = ndgrid(lon, lat, plev, time2);
    lonf = 0:2.5:357.5; nlonf = length(lonf); %144
    latf = 90:-2.5:-90; nlatf = length(latf); %73 is kernels
    plevf1 = ncread('/data/pub/kernels_YiH/toa/dp.nc', 'player'); % pay attention to plevf's range must smaller than plev's

    % kernel question contact McGill University Prof.Yi Huang Email: yi.huang@mcgill.ca
    plevf = plevf1(end:-1:1);
    plevfnum = length(plevf);
    [Xlonf, Ylatf, Zplevf, Ttimef] = ndgrid(lonf, latf, plevf, time2);
    q = interpn(Xlon, Ylat, Zplev, Ttime, q, Xlonf, Ylatf, Zplevf, Ttimef);
    ta = interpn(Xlon, Ylat, Zplev, Ttime, ta, Xlonf, Ylatf, Zplevf, Ttimef);

    [Xlon, Ylat, Ttime] = meshgrid(lat, lon, time2);
    [Xlonf, Ylatf, Ttimef] = meshgrid(latf, lonf, time2);
    ts = interp3(Xlon, Ylat, Ttime, ts, Xlonf, Ylatf, Ttimef);
    alb = interp3(Xlon, Ylat, Ttime, alb, Xlonf, Ylatf, Ttimef);

    startmonth = tLin.startmonth{p_1};
    %% Deseasonalize
    [dalb, Clim_alb] = monthlyAnomaly3D(nlonf, nlatf, time, alb, startmonth);
    [dts, Clim_ts] = monthlyAnomaly3D(nlonf, nlatf, time, ts, startmonth);
    [dq, Clim_q] = monthlyAnomaly4D(nlonf, nlatf, plevfnum, time, q, startmonth);
    [dta, Clim_ta] = monthlyAnomaly4D(nlonf, nlatf, plevfnum, time, ta, startmonth);

    %% Save to mat
    save([outpath, 'dvars.mat'], 'dalb', 'dts', 'dq', 'dta'...
        , 'Clim_alb', 'Clim_ts', 'Clim_q', 'Clim_ta'...
        , 'lonf', 'latf', 'plevf', 'time', 'readme')

    save([outpath1, 'vars.mat'], 'alb', 'ts', 'q', 'ta'...
        , 'lonf', 'latf', 'plevf', 'time', 'readme')
end
t=toc;disp(t)