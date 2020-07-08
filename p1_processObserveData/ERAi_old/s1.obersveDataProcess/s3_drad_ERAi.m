%% Combine the ERAi rad data from 2000-03 to 2018-02(18 years)
% contain allsky and clear sky.
clc; clear; tic;
%modify path first
inputpath1 = '/data1/liuyincheng/Observe-rawdata/ERAi/Rad_ori/';
inputpath2 = '/data1/liuyincheng/Observe-rawdata/ERAi/Heat_Flux/Heat_Flux_Accumlate/';

%% read raw data. 228 mean raw date 19yrs time line
dR_lw_toa_all0 = zeros(360, 181, 228);
dR_sw_toa_all0 = zeros(360, 181, 228);
dR_lw_sfc_all0 = zeros(360, 181, 228);
dR_sw_sfc_all0 = zeros(360, 181, 228);
dR_lw_toa_clr0 = zeros(360, 181, 228);
dR_sw_toa_clr0 = zeros(360, 181, 228);
dR_lw_sfc_clr0 = zeros(360, 181, 228);
dR_sw_sfc_clr0 = zeros(360, 181, 228);
lahf0 = zeros(360, 181, 228);
sehf0 = zeros(360, 181, 228);
slwr0 = zeros(360, 181, 228);
slwrd0 = zeros(360, 181, 228);
time0 = zeros(228, 1);
for ii = 0:18
    filename = strcat(inputpath1, 'rad', num2str((ii + 2000), '%04i'), '.nc');
    filename1 = strcat(inputpath2, 'hf', num2str((ii + 2000), '%04i'), '.nc');

    temp1 = ncread(filename, 'ttr'); % Top net thermal radiation,
    temp2 = ncread(filename, 'tsr'); % Top net solar radiation,
    temp3 = ncread(filename, 'str'); % Surface net thermal radiation,
    temp4 = ncread(filename, 'ssr'); % Surface net solar radiation,
    temp5 = ncread(filename, 'time'); % hours since 1900-01-01 00:00:00.0
    temp6 = ncread(filename, 'ttrc'); % Top net thermal radiation, clear sky
    temp7 = ncread(filename, 'tsrc'); % Top net solar radiation, clear sky
    temp8 = ncread(filename, 'strc'); % Surface net thermal radiation, clear sky
    temp9 = ncread(filename, 'ssrc'); % Surface net solar radiation, clear sky
    % heat flux
    temp10 = ncread(filename1, 'slhf'); % surface_upward_latent_heat_flux
    temp11 = ncread(filename1, 'sshf'); % surface_upward_sensible_heat_flux
    temp12 = ncread(filename1, 'str'); % surface_net_upward_longwave_flux(input rad + earth output rad)
    temp13 = ncread(filename1, 'strd'); % Surface thermal radiation downwards

    % only read  0-12 is ok
    dR_lw_toa_all0(:, :, ii * 12 + 1:(ii + 1) * 12) = temp1(:, :, 1:2:23);
    dR_sw_toa_all0(:, :, ii * 12 + 1:(ii + 1) * 12) = temp2(:, :, 1:2:23);
    dR_lw_sfc_all0(:, :, ii * 12 + 1:(ii + 1) * 12) = temp3(:, :, 1:2:23);
    dR_sw_sfc_all0(:, :, ii * 12 + 1:(ii + 1) * 12) = temp4(:, :, 1:2:23);
    time0(ii * 12 + 1:(ii + 1) * 12, 1) = temp5(1:2:23);
    dR_lw_toa_clr0(:, :, ii * 12 + 1:(ii + 1) * 12) = temp6(:, :, 1:2:23);
    dR_sw_toa_clr0(:, :, ii * 12 + 1:(ii + 1) * 12) = temp7(:, :, 1:2:23);
    dR_lw_sfc_clr0(:, :, ii * 12 + 1:(ii + 1) * 12) = temp8(:, :, 1:2:23);
    dR_sw_sfc_clr0(:, :, ii * 12 + 1:(ii + 1) * 12) = temp9(:, :, 1:2:23);
    lahf0(:, :, ii * 12 + 1:(ii + 1) * 12) = temp10(:, :, 1:2:23);
    sehf0(:, :, ii * 12 + 1:(ii + 1) * 12) = temp11(:, :, 1:2:23);
    slwr0(:, :, ii * 12 + 1:(ii + 1) * 12) = temp12(:, :, 1:2:23);
    slwrd0(:, :, ii * 12 + 1:(ii + 1) * 12) = temp13(:, :, 1:2:23);

end

% transfor mat time
temp.name = strcat(inputpath1, 'rad', num2str((4 + 2000), '%04i'), '.nc');
timeUnits = ncreadatt(temp.name, 'time', 'units');
% timeCalendar = ncreadatt(temp.name, 'time', 'calendar');
timeCalendar = 'gregorian';
time0 = cdfdate2num(timeUnits, timeCalendar, time0);


%% different time series, 1mean 2000-03 to 2018-02(18*12). 2 mean 200207-201706(15*12)
for p_1 = 1:2
    [readme, tLin] = observeParameters(p_1); % fuction% readme,  tLin,
    outpath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi/anomaly/');
    auto_mkdir(outpath)

    % find start and end month
    timeSpec = tLin.start{p_1};
    formatOut = 'yyyy-mm';
    timeStr = string(datestr(time0, formatOut));
    startT = find(timeStr == timeSpec);
    endT = startT + tLin.inter{p_1} - 1;

    % read specific time series
    dR_sw_toa_all = dR_sw_toa_all0(:, :, startT:endT);
    dR_lw_sfc_all = dR_lw_sfc_all0(:, :, startT:endT);
    dR_lw_toa_all = dR_lw_toa_all0(:, :, startT:endT);
    dR_sw_sfc_all = dR_sw_sfc_all0(:, :, startT:endT);
    dR_lw_toa_clr = dR_lw_toa_clr0(:, :, startT:endT);
    dR_sw_toa_clr = dR_sw_toa_clr0(:, :, startT:endT);
    dR_lw_sfc_clr = dR_lw_sfc_clr0(:, :, startT:endT);
    dR_sw_sfc_clr = dR_sw_sfc_clr0(:, :, startT:endT);
    slwru0 = slwr0 + slwrd0; % surface upward longwave flux
    dlahf = -lahf0(:, :, startT:endT); % define downward as positive
    dsehf = -sehf0(:, :, startT:endT); % define downward as positive
    dslwrd = slwrd0(:, :, startT:endT);
    dslwr = -slwr0(:, :, startT:endT); % define downward as positive
    dslwru = -slwru0(:, :, startT:endT); % define downward as positive

    time = time0(startT:endT, 1);
    ntime = length(time);

    % transform unite: W/m2
    dR_lw_toa_all = dR_lw_toa_all ./ (3600 * 24);
    dR_sw_toa_all = dR_sw_toa_all ./ (3600 * 24);
    dR_lw_sfc_all = dR_lw_sfc_all ./ (3600 * 24);
    dR_sw_sfc_all = dR_sw_sfc_all ./ (3600 * 24);
    dR_lw_toa_clr = dR_lw_toa_clr ./ (3600 * 24);
    dR_sw_toa_clr = dR_sw_toa_clr ./ (3600 * 24);
    dR_lw_sfc_clr = dR_lw_sfc_clr ./ (3600 * 24);
    dR_sw_sfc_clr = dR_sw_sfc_clr ./ (3600 * 24);
    dlahf = dlahf ./ (3600 * 24);
    dsehf = dsehf ./ (3600 * 24);
    dslwr = dslwr ./ (3600 * 24);
    dslwrd = dslwrd ./ (3600 * 24);
    dslwru = dslwru ./ (3600 * 24);

    %% Regrid to 2.5*2.5
    lon = 0:359;
    lat = 90:-1:-90;
    lonPlot = 0:2.5:357.5; nlonk = length(lonPlot); %144
    latPlot = 88.75:-2.5:-88.75; nlatk = length(latPlot); %72 is plot, 73 is kernel
    time2 = 1:ntime;

    [Xlon, Ylat, Ttime] = meshgrid(lat, lon, time2);
    [Xlonf, Ylatf, Ttimef] = meshgrid(latPlot, lonPlot, time2);

    dR_lw_toa_all = interp3(Xlon, Ylat, Ttime, dR_lw_toa_all, Xlonf, Ylatf, Ttimef);
    dR_sw_toa_all = interp3(Xlon, Ylat, Ttime, dR_sw_toa_all, Xlonf, Ylatf, Ttimef);
    dR_lw_sfc_all = interp3(Xlon, Ylat, Ttime, dR_lw_sfc_all, Xlonf, Ylatf, Ttimef);
    dR_sw_sfc_all = interp3(Xlon, Ylat, Ttime, dR_sw_sfc_all, Xlonf, Ylatf, Ttimef);
    dR_lw_toa_clr = interp3(Xlon, Ylat, Ttime, dR_lw_toa_clr, Xlonf, Ylatf, Ttimef);
    dR_sw_toa_clr = interp3(Xlon, Ylat, Ttime, dR_sw_toa_clr, Xlonf, Ylatf, Ttimef);
    dR_lw_sfc_clr = interp3(Xlon, Ylat, Ttime, dR_lw_sfc_clr, Xlonf, Ylatf, Ttimef);
    dR_sw_sfc_clr = interp3(Xlon, Ylat, Ttime, dR_sw_sfc_clr, Xlonf, Ylatf, Ttimef);

    dlahf = interp3(Xlon, Ylat, Ttime, dlahf, Xlonf, Ylatf, Ttimef);
    dsehf = interp3(Xlon, Ylat, Ttime, dsehf, Xlonf, Ylatf, Ttimef);
    dslwr = interp3(Xlon, Ylat, Ttime, dslwr, Xlonf, Ylatf, Ttimef);
    dslwrd = interp3(Xlon, Ylat, Ttime, dslwrd, Xlonf, Ylatf, Ttimef);
    dslwru = interp3(Xlon, Ylat, Ttime, dslwru, Xlonf, Ylatf, Ttimef);

    % Deseasonalized(var-mean)
    startmonth = tLin.startmonth{p_1};
    [dR_lw_toa_all, climdR_lw_toa_all] = monthlyAnomaly3D(nlonk, nlatk, time, dR_lw_toa_all, startmonth);
    [dR_sw_toa_all, climdR_sw_toa_all] = monthlyAnomaly3D(nlonk, nlatk, time, dR_sw_toa_all, startmonth);
    [dR_lw_sfc_all, climdR_lw_sfc_all] = monthlyAnomaly3D(nlonk, nlatk, time, dR_lw_sfc_all, startmonth);
    [dR_sw_sfc_all, climdR_sw_sfc_all] = monthlyAnomaly3D(nlonk, nlatk, time, dR_sw_sfc_all, startmonth);
    [dR_lw_toa_clr, climdR_lw_toa_clr] = monthlyAnomaly3D(nlonk, nlatk, time, dR_lw_toa_clr, startmonth);
    [dR_sw_toa_clr, climdR_sw_toa_clr] = monthlyAnomaly3D(nlonk, nlatk, time, dR_sw_toa_clr, startmonth);
    [dR_lw_sfc_clr, climdR_lw_sfc_clr] = monthlyAnomaly3D(nlonk, nlatk, time, dR_lw_sfc_clr, startmonth);
    [dR_sw_sfc_clr, climdR_sw_sfc_clr] = monthlyAnomaly3D(nlonk, nlatk, time, dR_sw_sfc_clr, startmonth);

    [dlahf, Clim_dlahf] = monthlyAnomaly3D(nlonk, nlatk, time, dlahf, startmonth);
    [dsehf, Clim_dsehf] = monthlyAnomaly3D(nlonk, nlatk, time, dsehf, startmonth);
    [dslwr, Clim_dslwr] = monthlyAnomaly3D(nlonk, nlatk, time, dslwr, startmonth);
    [dslwrd, Clim_dslwrd] = monthlyAnomaly3D(nlonk, nlatk, time, dslwrd, startmonth);
    [dslwru, Clim_dslwru] = monthlyAnomaly3D(nlonk, nlatk, time, dslwru, startmonth);
    %% Save to mat
    % outpath = '/data1/liuyincheng/Observe-process/ERAi-process/';
    save([outpath, 'drad.mat'], 'dR_lw_toa_all', 'dR_sw_toa_all', 'dR_lw_sfc_all', 'dR_sw_sfc_all'...
        , 'dR_lw_toa_clr', 'dR_sw_toa_clr', 'dR_lw_sfc_clr', 'dR_sw_sfc_clr'...
        , 'climdR_lw_toa_all', 'climdR_sw_toa_all', 'climdR_lw_sfc_all', 'climdR_sw_sfc_all'...
        , 'climdR_lw_toa_clr', 'climdR_sw_toa_clr', 'climdR_lw_sfc_clr', 'climdR_sw_sfc_clr'...
        , 'lonPlot', 'latPlot', 'time', 'readme')
    save([outpath, 'dheatFlux.mat'], 'dlahf', 'dsehf', 'dslwr', 'dslwrd', 'dslwru'...
        , 'Clim_dlahf', 'Clim_dsehf', 'Clim_dslwr', 'Clim_dslwrd', 'Clim_dslwru'...
        , 'lonPlot', 'latPlot', 'time', 'readme')

end

t = toc; disp(t)
