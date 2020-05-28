% This program is to calculate the radiation anomalies from CERES/EBAF4.0
clc; clear; tic;
vector_var_str{1} = 'sfc';
vector_var_str{2} = 'toa';
%modify path first
inputpath = '/data1/liuyincheng/Observe-rawdata/CERES_EBAF/4.0/';
filename1 = 'CERES_EBAF-Surface_Ed4.0_Subset_200003-201802.nc';
filename2 = 'CERES_EBAF-TOA_Ed4.0_Subset_200003-201802.nc';
filename = strcat(inputpath, filename1);
lon = 0.5:1:359.5; nlon = length(lon); lon = double(lon);
lat = -89.5:1:89.5; nlat = length(lat); lat = double(lat);

% transfor mat time
time = double(ncread(filename, 'time')); % days since 2000-03-01 00:00:00
timeUnits = ncreadatt(filename, 'time', 'units');
timeCalendar = 'gregorian';
time = cdfdate2num(timeUnits, timeCalendar, time);
ntime = length(time);

%% read raw data
%% SFC
% All-sky
dR_tot_sfc_all = ncread(filename, 'sfc_net_tot_all_mon');
dR_lw_sfc_all = ncread(filename, 'sfc_net_lw_all_mon');
dR_sw_sfc_all = ncread(filename, 'sfc_net_sw_all_mon');
dR_cre_sfc_all = ncread(filename, 'sfc_cre_net_tot_mon'); % Cloud radiative effect
dR_crelw_sfc_all = ncread(filename, 'sfc_cre_net_lw_mon');
dR_cresw_sfc_all = ncread(filename, 'sfc_cre_net_sw_mon');
% Clear-sky
dR_tot_sfc_clr = ncread(filename, 'sfc_net_tot_clr_mon');
dR_lw_sfc_clr = ncread(filename, 'sfc_net_lw_clr_mon');
dR_sw_sfc_clr = ncread(filename, 'sfc_net_sw_clr_mon');

%% TOA
filename = strcat(inputpath, filename2);
% allsky
dR_tot_toa_all = ncread(filename, 'toa_net_all_mon');
dR_lw_toa_all = -ncread(filename, 'toa_lw_all_mon'); % define positive as downward
dR_sw_toa_all = -ncread(filename, 'toa_sw_all_mon');
dR_cre_toa_all = ncread(filename, 'toa_cre_net_mon'); % Cloud radiative effect
dR_crelw_toa_all = ncread(filename, 'toa_cre_lw_mon');
dR_cresw_toa_all = ncread(filename, 'toa_cre_sw_mon');
% clear sky
dR_tot_toa_clr = ncread(filename, 'toa_net_clr_mon');
dR_lw_toa_clr = -ncread(filename, 'toa_lw_clr_mon'); % define positive as downward
dR_sw_toa_clr = -ncread(filename, 'toa_sw_clr_mon');

month = 1:12;
%% Regrid (Because of unite plot)
lonPlot = 0:2.5:357.5; nlonf = length(lonPlot);
latPlot = 88.75:-2.5:-88.75; nlatf = length(latPlot);
% For interpolate(??????????lonf???lon?????(lon????????), latf???????lat????????????????)
lon(2:end + 1) = lon;
lon(1) = -0.5;

dR_tot_sfc_all(2:end + 1, :, :) = dR_tot_sfc_all;
dR_tot_sfc_all(1, :, :) = dR_tot_sfc_all(end, :, :);
dR_lw_sfc_all(2:end + 1, :, :) = dR_lw_sfc_all;
dR_lw_sfc_all(1, :, :) = dR_lw_sfc_all(end, :, :);
dR_sw_sfc_all(2:end + 1, :, :) = dR_sw_sfc_all;
dR_sw_sfc_all(1, :, :) = dR_sw_sfc_all(end, :, :);
dR_cre_sfc_all(2:end + 1, :, :) = dR_cre_sfc_all;
dR_cre_sfc_all(1, :, :) = dR_cre_sfc_all(end, :, :);
dR_crelw_sfc_all(2:end + 1, :, :) = dR_crelw_sfc_all;
dR_crelw_sfc_all(1, :, :) = dR_crelw_sfc_all(end, :, :);
dR_cresw_sfc_all(2:end + 1, :, :) = dR_cresw_sfc_all;
dR_cresw_sfc_all(1, :, :) = dR_cresw_sfc_all(end, :, :);

dR_tot_sfc_clr(2:end + 1, :, :) = dR_tot_sfc_clr;
dR_tot_sfc_clr(1, :, :) = dR_tot_sfc_clr(end, :, :);
dR_lw_sfc_clr(2:end + 1, :, :) = dR_lw_sfc_clr;
dR_lw_sfc_clr(1, :, :) = dR_lw_sfc_clr(end, :, :);
dR_sw_sfc_clr(2:end + 1, :, :) = dR_sw_sfc_clr;
dR_sw_sfc_clr(1, :, :) = dR_sw_sfc_clr(end, :, :);

dR_tot_toa_all(2:end + 1, :, :) = dR_tot_toa_all;
dR_tot_toa_all(1, :, :) = dR_tot_toa_all(end, :, :);
dR_lw_toa_all(2:end + 1, :, :) = dR_lw_toa_all;
dR_lw_toa_all(1, :, :) = dR_lw_toa_all(end, :, :);
dR_sw_toa_all(2:end + 1, :, :) = dR_sw_toa_all;
dR_sw_toa_all(1, :, :) = dR_sw_toa_all(end, :, :);
dR_cre_toa_all(2:end + 1, :, :) = dR_cre_toa_all;
dR_cre_toa_all(1, :, :) = dR_cre_toa_all(end, :, :);
dR_crelw_toa_all(2:end + 1, :, :) = dR_crelw_toa_all;
dR_crelw_toa_all(1, :, :) = dR_crelw_toa_all(end, :, :);
dR_cresw_toa_all(2:end + 1, :, :) = dR_cresw_toa_all;
dR_cresw_toa_all(1, :, :) = dR_cresw_toa_all(end, :, :);

dR_tot_toa_clr(2:end + 1, :, :) = dR_tot_toa_clr;
dR_tot_toa_clr(1, :, :) = dR_tot_toa_clr(end, :, :);
dR_lw_toa_clr(2:end + 1, :, :) = dR_lw_toa_clr;
dR_lw_toa_clr(1, :, :) = dR_lw_toa_clr(end, :, :);
dR_sw_toa_clr(2:end + 1, :, :) = dR_sw_toa_clr;
dR_sw_toa_clr(1, :, :) = dR_sw_toa_clr(end, :, :);

% Regraded to 2.5*2.5
time2 = 1:ntime;
[Xlon, Ylat, Ttime] = meshgrid(lat, lon, time2);
[Xlonf, Ylatf, Ttimef] = meshgrid(latPlot, lonPlot, time2);

dR_tot_sfc_all = interp3(Xlon, Ylat, Ttime, dR_tot_sfc_all, Xlonf, Ylatf, Ttimef);
dR_tot_toa_all = interp3(Xlon, Ylat, Ttime, dR_tot_toa_all, Xlonf, Ylatf, Ttimef);
dR_tot_sfc_clr = interp3(Xlon, Ylat, Ttime, dR_tot_sfc_clr, Xlonf, Ylatf, Ttimef);
dR_tot_toa_clr = interp3(Xlon, Ylat, Ttime, dR_tot_toa_clr, Xlonf, Ylatf, Ttimef);

dR_lw_sfc_all = interp3(Xlon, Ylat, Ttime, dR_lw_sfc_all, Xlonf, Ylatf, Ttimef);
dR_lw_toa_all = interp3(Xlon, Ylat, Ttime, dR_lw_toa_all, Xlonf, Ylatf, Ttimef);
dR_lw_sfc_clr = interp3(Xlon, Ylat, Ttime, dR_lw_sfc_clr, Xlonf, Ylatf, Ttimef);
dR_lw_toa_clr = interp3(Xlon, Ylat, Ttime, dR_lw_toa_clr, Xlonf, Ylatf, Ttimef);

dR_sw_sfc_all = interp3(Xlon, Ylat, Ttime, dR_sw_sfc_all, Xlonf, Ylatf, Ttimef);
dR_sw_toa_all = interp3(Xlon, Ylat, Ttime, dR_sw_toa_all, Xlonf, Ylatf, Ttimef);
dR_sw_sfc_clr = interp3(Xlon, Ylat, Ttime, dR_sw_sfc_clr, Xlonf, Ylatf, Ttimef);
dR_sw_toa_clr = interp3(Xlon, Ylat, Ttime, dR_sw_toa_clr, Xlonf, Ylatf, Ttimef);

dR_cre_toa_all = interp3(Xlon, Ylat, Ttime, dR_cre_toa_all, Xlonf, Ylatf, Ttimef);
dR_crelw_toa_all = interp3(Xlon, Ylat, Ttime, dR_crelw_toa_all, Xlonf, Ylatf, Ttimef);
dR_cresw_toa_all = interp3(Xlon, Ylat, Ttime, dR_cresw_toa_all, Xlonf, Ylatf, Ttimef);

%% different time series, 1mean 2000-03 to 2018-02(18*12). 2 mean 200207-201706(15*12)
for p_1 = 1:2
    [readme, tLin] = observeParameters(p_1); % fuction% readme,  tLin,
    outpath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'CERES/anomaly/');
    auto_mkdir(outpath)

    % find start and end month
    timeSpec = tLin.start{p_1};
    formatOut = 'yyyy-mm';
    timeStr = string(datestr(time, formatOut));
    startT = find(timeStr == timeSpec);
    endT = startT + tLin.inter{p_1} - 1;

    % read specific time series
    dR_tot_sfc_all = dR_tot_sfc_all(:, :, startT:endT);
    dR_lw_sfc_all = dR_lw_sfc_all(:, :, startT:endT);
    dR_sw_sfc_all = dR_sw_sfc_all(:, :, startT:endT);
    dR_cre_sfc_all = dR_cre_sfc_all(:, :, startT:endT);
    dR_crelw_sfc_all = dR_crelw_sfc_all(:, :, startT:endT);
    dR_cresw_sfc_all = dR_cresw_sfc_all(:, :, startT:endT);

    dR_tot_sfc_clr = dR_tot_sfc_clr(:, :, startT:endT);
    dR_lw_sfc_clr = dR_lw_sfc_clr(:, :, startT:endT);
    dR_sw_sfc_clr = dR_sw_sfc_clr(:, :, startT:endT);

    dR_tot_toa_all = dR_tot_toa_all(:, :, startT:endT);
    dR_lw_toa_all = dR_lw_toa_all(:, :, startT:endT);
    dR_sw_toa_all = dR_sw_toa_all(:, :, startT:endT);
    dR_cre_toa_all = dR_cre_toa_all(:, :, startT:endT);
    dR_crelw_toa_all = dR_crelw_toa_all(:, :, startT:endT);
    dR_cresw_toa_all = dR_cresw_toa_all(:, :, startT:endT);

    dR_tot_toa_clr = dR_tot_toa_clr(:, :, startT:endT);
    dR_lw_toa_clr = dR_lw_toa_clr(:, :, startT:endT);
    dR_sw_toa_clr = dR_sw_toa_clr(:, :, startT:endT);

    time = time(startT:endT, 1);
    ntime = length(time);

    %% Deseasonalized Anomaly
    startmonth = tLin.startmonth{p_1};
    [dR_tot_sfc_all, Clim_dR_tot_sfc_all] = monthlyAnomaly3D(nlonf, nlatf, time, dR_tot_sfc_all, startmonth);
    [dR_tot_sfc_clr, Clim_dR_tot_sfc_clr] = monthlyAnomaly3D(nlonf, nlatf, time, dR_tot_sfc_clr, startmonth);
    [dR_lw_sfc_all, Clim_dR_lw_sfc_all] = monthlyAnomaly3D(nlonf, nlatf, time, dR_lw_sfc_all, startmonth);
    [dR_lw_sfc_clr, Clim_dR_lw_sfc_clr] = monthlyAnomaly3D(nlonf, nlatf, time, dR_lw_sfc_clr, startmonth);
    [dR_sw_sfc_all, Clim_dR_sw_sfc_all] = monthlyAnomaly3D(nlonf, nlatf, time, dR_sw_sfc_all, startmonth);
    [dR_sw_sfc_clr, Clim_dR_sw_sfc_clr] = monthlyAnomaly3D(nlonf, nlatf, time, dR_sw_sfc_clr, startmonth);

    [dR_tot_toa_all, Clim_dR_tot_toa_all] = monthlyAnomaly3D(nlonf, nlatf, time, dR_tot_toa_all, startmonth);
    [dR_tot_toa_clr, Clim_dR_tot_toa_clr] = monthlyAnomaly3D(nlonf, nlatf, time, dR_tot_toa_clr, startmonth);
    [dR_lw_toa_clr, Clim_dR_lw_toa_clr] = monthlyAnomaly3D(nlonf, nlatf, time, dR_lw_toa_clr, startmonth);
    [dR_lw_toa_all, Clim_dR_lw_toa_all] = monthlyAnomaly3D(nlonf, nlatf, time, dR_lw_toa_all, startmonth);
    [dR_sw_toa_clr, Clim_dR_sw_toa_clr] = monthlyAnomaly3D(nlonf, nlatf, time, dR_sw_toa_clr, startmonth);
    [dR_sw_toa_all, Clim_dR_sw_toa_all] = monthlyAnomaly3D(nlonf, nlatf, time, dR_sw_toa_all, startmonth);
    [dR_cre_toa_all, Clim_dR_cre_toa_all] = monthlyAnomaly3D(nlonf, nlatf, time, dR_cre_toa_all, startmonth);
    [dR_crelw_toa_all, Clim_dR_crelw_toa_all] = monthlyAnomaly3D(nlonf, nlatf, time, dR_crelw_toa_all, startmonth);
    [dR_cresw_toa_all, Clim_dR_cresw_toa_all] = monthlyAnomaly3D(nlonf, nlatf, time, dR_cresw_toa_all, startmonth);

    %% save mat files
    save([outpath, 'drad.mat'], 'dR_lw_toa_all', 'dR_sw_toa_all', 'dR_lw_sfc_all', 'dR_sw_sfc_all'...
        , 'dR_lw_toa_clr', 'dR_sw_toa_clr', 'dR_lw_sfc_clr', 'dR_sw_sfc_clr'...
        , 'Clim_dR_lw_toa_all', 'Clim_dR_sw_toa_all', 'Clim_dR_lw_sfc_all', 'Clim_dR_sw_sfc_all'...
        , 'Clim_dR_lw_toa_clr', 'Clim_dR_sw_toa_clr', 'Clim_dR_lw_sfc_clr', 'Clim_dR_sw_sfc_clr'...
        , 'lonPlot', 'latPlot', 'time', 'readme'...
        , 'dR_tot_sfc_all', 'dR_tot_sfc_clr', 'dR_tot_toa_all', 'dR_tot_toa_clr', 'dR_cre_toa_all', 'dR_crelw_toa_all', 'dR_cresw_toa_all'...
        , 'Clim_dR_tot_sfc_all', 'Clim_dR_tot_sfc_clr', 'Clim_dR_tot_toa_all', 'Clim_dR_tot_toa_clr', 'Clim_dR_cre_toa_all', 'Clim_dR_crelw_toa_all', 'Clim_dR_cresw_toa_all')
end

t = toc; disp(t)
