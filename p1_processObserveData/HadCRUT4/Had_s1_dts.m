%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-07-13 10:31:04
% LastEditTime : 2020-07-13 11:19:44
% LastEditors  : LYC
% Description  :
% FilePath     : /code/p1_processObserveData/HadCRUT4/Had_s1_dts.m
%
%%---------------------------------------------------------
clc; clear; tic
% HadCRUT4 data
hadDataPath = '/data1/liuyincheng/Observe-rawdata/HadCRUT4/HadCRUT.4.6.0.0.median.nc';
dts_ori = ncread(hadDataPath, 'temperature_anomaly');
lon_h = ncread(hadDataPath, 'longitude'); nlonh = length(lon_h); % -177.5:5:177.5
lat_h = ncread(hadDataPath, 'latitude'); nlath = length(lat_h);

% transfor mat time
time_h = ncread(hadDataPath, 'time'); % 1979/03 - 2018/02
timeUnits_h = ncreadatt(hadDataPath, 'time', 'units');
timeCalendar_h = 'gregorian';
time_h = double(cdfdate2num(timeUnits_h, timeCalendar_h, time_h));
ntimeh = length(time_h);

[readme, level, tLin, vars] = obsParameters('CERES');
%% different time series, 1mean 2000-03 to 2018-02(18*12). 2 mean 200207-201706(15*12)
for p_1 = 1:2
    anomPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'HadCRUT4', level.standVarPath{2}); % /anomaly
    auto_mkdir(anomPath);

    %% calculate the ts anomalies from HadCRUT4
    % find start and end month\
    timeSpec = tLin.start{p_1};
    formatOut = 'yyyy-mm';
    timeStr_h = string(datestr(time_h, formatOut));
    startT_h = find(timeStr_h == timeSpec);
    endT_h = startT_h + tLin.inter{p_1} - 1;

    dts0 = dts_ori(:, :, startT_h:endT_h); % 1803:2000.03, 2018:2018.02, 72*36*216 surface temperature anomaly, Global land
    time_h = time_h(startT_h:endT_h);
    ntimeh = length(time_h);

    % Missing too many data,so it will be eliminated from the time series after missing only one test
    test = dts0;
    test1 = sum(isnan(dts0), 3);

    for ii = 1:ntimeh
        temp1 = test(:, :, ii);
        temp1(test1 ~= 0) = 0;
        test(:, :, ii) = temp1;
    end

    dts0 = test;
    startMonth = tLin.startMonth{p_1};
    [trendm_dtsHad, trends_dtsHad, trendyr_dtsHad, ~, ~] = autoCalTrend(dts0, nlonh, nlath, time_h, startMonth);
    save([anomPath, 'trend_dtsHad4.0.mat'], 'lon_h', 'lat_h', 'time_h', ...
        'trendm_dtsHad', 'trends_dtsHad', 'trendyr_dtsHad')
end

t = toc; disp(t)
