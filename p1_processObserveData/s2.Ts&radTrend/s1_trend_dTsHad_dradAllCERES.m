%% Combine the surface CERES cldrad and HadCRUT4 Ts data from 2000-03 to 2018-02 (18 years)
% attention:  matlab's calendar doesnt consist with the nc file, so result nc show's wrong time in ncview, but doest affect result
% only matter that data from 2000-03 to 2018-02 (18 years), not from 2000-01 to 2017-12
%
% this function only cal surface varity anomaly (including cloud and clc):
% strd(Surface thermal radiation downwards),
% ssr(surface_net_downward_shortwave_flux),
% ts
% cal time seris:
% every month
%

clc; clear; tic
% -------section 1: calculate the radiation anomalies from CERES/EBAF4.0-----------
vector_var_str{1} = 'sfc';
% vector_var_str{2} = 'toa';
var = char(vector_var_str{1});

%modify path first
inputpath = '/data1/liuyincheng/Observe-rawdata/CERES_EBAF/4.0/';
filename1 = 'CERES_EBAF-Surface_Ed4.0_Subset_200003-201802.nc';
filename = strcat(inputpath, filename1);

%% read raw data
% All-sky
lon = 0.5:1:359.5; lon = double(lon);
lat = -89.5:1:89.5; lat = double(lat);
dRnetCE = ncread(filename, 'sfc_net_tot_all_mon');
dRnetCE_lw_cld = ncread(filename, 'sfc_net_lw_all_mon');
dRdwnCE_lw_cld = ncread(filename, 'sfc_lw_down_all_mon');
dRnetCE_sw_cld = ncread(filename, 'sfc_net_sw_all_mon');
% transfor mat time
time = double(ncread(filename, 'time')); % days since 2000-03-01 00:00:00
timeUnits = ncreadatt(filename, 'time', 'units');
timeCalendar = 'gregorian';
time = cdfdate2num(timeUnits, timeCalendar, time);
ntime = length(time);
% For interpolate
for ll = 360:-1:1
    dRnetCE(ll + 1, :, :) = dRnetCE(ll, :, :);
    dRnetCE_lw_cld(ll + 1, :, :) = dRnetCE_lw_cld(ll, :, :);
    dRdwnCE_lw_cld(ll + 1, :, :) = dRdwnCE_lw_cld(ll, :, :);
    dRnetCE_sw_cld(ll + 1, :, :) = dRnetCE_sw_cld(ll, :, :);

end

dRnetCE(1, :, :) = dRnetCE(361, :, :);
dRnetCE_lw_cld(1, :, :) = dRnetCE_lw_cld(361, :, :);
dRdwnCE_lw_cld(1, :, :) = dRdwnCE_lw_cld(361, :, :);
dRnetCE_sw_cld(1, :, :) = dRnetCE_sw_cld(361, :, :);
lon(2:end + 1) = lon;
lon(1) = -0.5;
nlat = length(lat);
nlon = length(lon);
%% Regrid (unite plot)
% Regraded to 2.5*2.5
lonPlot = 0:2.5:357.5; nlonPlot = length(lonPlot);
latPlot = 88.75:-2.5:-88.75; nlatPlot = length(latPlot);
time2 = 1:ntime;
[Xlon, Ylat, Ttime] = meshgrid(lat, lon, time2);
[Xlonf, Ylatf, Ttimef] = meshgrid(latPlot, lonPlot, time2);

dRnetCE = interp3(Xlon, Ylat, Ttime, dRnetCE, Xlonf, Ylatf, Ttimef);
dRnetCE_lw_cld = interp3(Xlon, Ylat, Ttime, dRnetCE_lw_cld, Xlonf, Ylatf, Ttimef);
dRdwnCE_lw_cld = interp3(Xlon, Ylat, Ttime, dRdwnCE_lw_cld, Xlonf, Ylatf, Ttimef);
dRnetCE_sw_cld = interp3(Xlon, Ylat, Ttime, dRnetCE_sw_cld, Xlonf, Ylatf, Ttimef);

% HadCRUT4 data
filename2 = '/data1/liuyincheng/Observe-rawdata/HadCRUT4/HadCRUT.4.6.0.0.median.nc';
dts0 = ncread(filename2, 'temperature_anomaly');
lonh = ncread(filename2, 'longitude'); nlonh = length(lonh); % -177.5:5:177.5
lath = ncread(filename2, 'latitude'); nlath = length(lath);
% transfor mat time
timeh = ncread(filename2, 'time'); % 1979/03 - 2018/02
timeUnits_h = ncreadatt(filename, 'time', 'units');
timeCalendar_h = 'gregorian';
timeh = double(cdfdate2num(timeUnits, timeCalendar, timeh));
ntimeh = length(timeh);

%% different time series, 1mean 2000-03 to 2018-02(18*12). 2 mean 200207-201706(15*12)
for p_1 = 1:2
    [readme, tLin] = observeParameters(p_1); % fuction% readme,  tLin,
    outpath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'CERES/anomaly_trend/');
    auto_mkdir(outpath)

    % find start and end month
    timeSpec = tLin.start{p_1}; 
    formatOut = 'yyyy-mm';
    timeStr = string(datestr(time, formatOut));
    startT = find(timeStr == timeSpec);
    endT = startT + tLin.inter{p_1} - 1;

    % read specific time series
    dRnetCE = dRnetCE(:, :, startT:endT);
    dRnetCE_lw_cld = dRnetCE_lw_cld(:, :, startT:endT);
    dRdwnCE_lw_cld = dRdwnCE_lw_cld(:, :, startT:endT);
    dRnetCE_sw_cld = dRnetCE_sw_cld(:, :, startT:endT);
    time = time(startT:endT);
    ntime = length(time);
    month = 1:12;

    %% Deseasonalized Anomaly
    startmonth = tLin.startmonth{p_1};
    [dRnetCE, Clim_dRnet] = monthlyAnomaly3D(nlonPlot, nlatPlot, time, dRnetCE, startmonth);
    [dRnetCE_lw_cld, Clim_dRnet_lw] = monthlyAnomaly3D(nlonPlot, nlatPlot, time, dRnetCE_lw_cld, startmonth);
    [dRdwnCE_lw_cld, Clim_dRdwn_lw] = monthlyAnomaly3D(nlonPlot, nlatPlot, time, dRdwnCE_lw_cld, startmonth);
    [dRnetCE_sw_cld, Clim_dRnet_sw] = monthlyAnomaly3D(nlonPlot, nlatPlot, time, dRnetCE_sw_cld, startmonth);

    % -------section 2: calculate the ts anomalies from HadCRUT4-----------
    % find start and end month
    timeStr_h = string(datestr(timeh, formatOut));
    startT_h = find(timeStr_h == timeSpec);
    endT_h = startT_h + tLin.inter{p_1} - 1;

    dts0_had = dts0(:, :, startT_h:endT_h); % 1803:2000.03, 2018:2018.02, 72*36*216 surface temperature anomaly, Global land
    timeh = timeh(startT_h:endT_h);
    ntimeh = length(timeh);

    % Missing too many data,so it will be eliminated from the time series after missing only one test
    test = dts0_had;
    test1 = sum(isnan(dts0_had), 3);

    for ii = 1:ntimeh
        temp1 = test(:, :, ii);
        temp1(test1 ~= 0) = 0;
        test(:, :, ii) = temp1;
    end

    dts0_had = test;

    % -------section 3: calculate rhs(LW_down+SW) and ts anomalies' trend(mon,season,year)----------
    % cal the rhs(LW_down+SW)
    CErhs = dRdwnCE_lw_cld + dRnetCE_sw_cld;
    % cal the ts trend
    trendm_dTsHad = zeros(nlonh, nlath, 12, 2); trendm_CErhs = zeros(nlonPlot, nlatPlot, 12, 2); % trend, 12 months
    % im==1 is MAM, im==2 is JJA, im==3 is SON, im==4 is DJF
    trends_dTsHad = zeros(nlonh, nlath, 4); trends_CErhs = zeros(nlonPlot, nlatPlot, 4); % trend, 4 seasons

    p_dTsHad = zeros(nlonh, nlath, 12); p_CErhs = zeros(nlonPlot, nlatPlot, 12); % f test
    cons_dTsHad = zeros(nlonh, nlath, 12); cons_CErhs = zeros(nlonPlot, nlatPlot, 12); % mean
    % months mean trend(unite:per day)
    for i = 1:nlonPlot

        for j = 1:nlatPlot
            [~, trendm_CErhs(i, j, :, :), cons_CErhs(i, j, :), p_CErhs(i, j, :)] = detrend_yan(CErhs(i, j, :), time);
        end

    end

    for i = 1:nlonh

        for j = 1:nlath
            [~, trendm_dTsHad(i, j, :, :), cons_dTsHad(i, j, :), p_dTsHad(i, j, :)] = detrend_yan(dts0_had(i, j, :), timeh);
        end

    end

    % years mean trend (unite:per day)
    trendyr_dTsHad = squeeze(mean(trendm_dTsHad(:, :, :, 1), 3));
    trendyr_CErhs = squeeze(mean(trendm_CErhs(:, :, :, 1), 3));

    % ----save-----
    outfilename1 = [outpath, 'dradCERES_dTsHad4.0.mat'];
    save(outfilename1, 'lonPlot', 'latPlot', 'dRnetCE', 'dRnetCE_lw_cld', 'dRdwnCE_lw_cld', 'dRnetCE_sw_cld', 'Clim_dRnet', 'Clim_dRnet_lw', 'Clim_dRdwn_lw', 'Clim_dRnet_sw', 'time', ...
        'lonh', 'lath', 'dts0_had', 'timeh')

    outfilename2 = [outpath, 'trend_dradCERES_dTsHad4.0.mat'];
    save(outfilename2, 'lonPlot', 'latPlot', 'lonh', 'lath', 'time', 'timeh', ...
        'trendm_dTsHad', 'trendm_CErhs', 'trendyr_dTsHad', 'trendyr_CErhs', ...
        'p_dTsHad', 'p_CErhs', 'cons_dTsHad', 'cons_CErhs')
end

t = toc; disp(t)
