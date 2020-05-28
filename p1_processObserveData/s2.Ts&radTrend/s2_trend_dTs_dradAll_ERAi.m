% x_trend_Ts_radcld_lyc.m
% cal the trend, if figure show the result is good, then cal the the coef
% note monthly data from 2000-03 to 2018-02
% trend1:(12,2) k and b, from 1-12 month
% cons: (12,1) means of month in 18 years
% pvalue: (12,1)
% time seris:
% 1.every month
% 2.MAM,JJA,SON,DJF
% 3.year

clc; clear; tic;


%% different time series, 1mean 2000-03 to 2018-02(18*12). 2 mean 200207-201706(15*12)
for p_1 = 1:2
    [readme, tLin] = observeParameters(p_1); % fuction% readme,  tLin,
    startmonth = tLin.startmonth{p_1};

    %modify path first
    inputpath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi/anomaly/dts_dradAll_ERAi.nc');

    outpath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi/anomaly_trend/');
    auto_mkdir(outpath)

    time = ncread(inputpath, 'time'); ntime = length(time);
    lon = ncread(inputpath, 'longitude'); nlon = length(lon);
    lat = ncread(inputpath, 'latitude'); nlat = length(lat);

    dts = ncread(inputpath, 'dts_sfc'); % skin temperature K
    lwdow_anom = ncread(inputpath, 'dRdownlw_sfcAll'); % Radiation Anomaly (downwards lw) at SFC W/m*2
    swnet_anom = ncread(inputpath, 'dRnetsw_sfcAll'); % Radiation Anomaly (net sw) at SFC
    drhs = lwdow_anom + swnet_anom; % R.H.S = lwdow_anom + swnet_anom

    trendm_dts = zeros(nlon, nlat, 12, 2); trendm_drhs = zeros(nlon, nlat, 12, 2); % trend, 12 months
    % im==1 is MAM, im==2 is JJA, im==3 is SON, im==4 is DJF

    p_ts = zeros(nlon, nlat, 12); p_rhs = zeros(nlon, nlat, 12); % f test
    cons_ts = zeros(nlon, nlat, 12); cons_rhs = zeros(nlon, nlat, 12); % mean
    % months mean trend(unite:per day)
    for i = 1:nlon

        for j = 1:nlat
            [~, trendm_dts(i, j, :, :), cons_ts(i, j, :), p_ts(i, j, :)] = detrend_yc(dts(i, j, :), time,startmonth);
            [~, trendm_drhs(i, j, :, :), cons_rhs(i, j, :), p_rhs(i, j, :)] = detrend_yc(drhs(i, j, :), time,startmonth);
        end

    end


    % years mean trend (unite:per day)
    trendyr_dts = squeeze(mean(trendm_dts(:, :, :, 1), 3));
    trendyr_drhs = squeeze(mean(trendm_drhs(:, :, :, 1), 3));

    %----------------------save varity----------------------------------
    outfilename = [outpath, 'trend_dtsdradERAi.mat'];
    save(outfilename, 'lon', 'lat', 'time', 'trendm_dts', 'trendm_drhs', 'cons_ts', 'cons_rhs', 'p_ts', 'p_rhs', ...
        'trendyr_dts', 'trendyr_drhs')
end

t = toc; disp(t)
