%% Combine the ERAi rad and Ts data from 2000-03 to 2018-02 (18 years)
% this function only cal surface varity anomaly (including cloud and clc):
% strd(Surface thermal radiation downwards),
% ssr(surface_net_downward_shortwave_flux),
% ts
% cal time seris:
% every month
%

clc; clear; tic;
% constant
sigma = 5.67e-8; % Stefan-Boltzmann constant: Wm-2K-4

%modify path first
inputpath1 = '/data1/liuyincheng/Observe-rawdata/ERAi/Rad_ori/';
inputpath2 = '/data1/liuyincheng/Observe-rawdata/ERAi/Heat_Flux/Heat_Flux_Accumlate/';
inputpath3 = '/data1/liuyincheng/Observe-rawdata/ERAi/Ts/';

%% read raw data. 228 mean raw date 19yrs time line
dR_downlw_sfc_cld0 = zeros(360, 181, 228);
dts_sfc0 = zeros(360, 181, 228);
dR_netsw_sfc_cld0 = zeros(360, 181, 228);
time0 = zeros(228, 1);
for ii = 0:18
    % inputpath1 = strcat('E:\xieyan\ERAi\Rad_ori\rad2004.nc');
    filename1 = strcat(inputpath1, 'rad', num2str((ii + 2000), '%04i'), '.nc');
    filename2 = strcat(inputpath2, 'hf', num2str((ii + 2000), '%04i'), '.nc');
    filename3 = strcat(inputpath3, 't', num2str((ii + 2000), '%04i'), '.nc');
    temp1 = ncread(filename2, 'strd'); % Surface thermal radiation downwards
    temp2 = ncread(filename3, 'skt'); % ts
    temp3 = ncread(filename1, 'ssr'); % surface_net_downward_shortwave_flux
    temp4 = ncread(filename1, 'time'); % hours since 1900-01-01 00:00:00.0
    % sum up the cumulative radiation of time step 0-12 and 12-24
    dR_downlw_sfc_cld0(:, :, ii * 12 + 1:(ii + 1) * 12) = temp1(:, :, 1:2:23);
    dts_sfc0(:, :, ii * 12 + 1:(ii + 1) * 12) = temp2(:, :, :);
    dR_netsw_sfc_cld0(:, :, ii * 12 + 1:(ii + 1) * 12) = temp3(:, :, 1:2:23);
    time0(ii * 12 + 1:(ii + 1) * 12, 1) = temp4(1:2:23);
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
    dRdownlw_sfcAll = dR_downlw_sfc_cld0(:, :, startT:endT);
    dts_sfc = dts_sfc0(:, :, startT:endT);
    dRnetsw_sfcAll = dR_netsw_sfc_cld0(:, :, startT:endT);
    time = time0(startT:endT, 1);
    ntime = length(time);
    % transform unite
    dRdownlw_sfcAll = dRdownlw_sfcAll ./ (3600 * 24);
    dRnetsw_sfcAll = dRnetsw_sfcAll ./ (3600 * 24);
    %% Regrid to 2.5*2.5
    lon = 0:359;
    lat = 90:-1:-90;
    lonPlot = 0:2.5:357.5; nlonf = length(lonPlot);
    latPlot = 88.75:-2.5:-88.75; nlatf = length(latPlot);
    time2 = 1:ntime;

    [Xlon, Ylat, Ttime] = meshgrid(lat, lon, time2);
    [Xlonf, Ylatf, Ttimef] = meshgrid(latPlot, lonPlot, time2);

    dRdownlw_sfcAll = interp3(Xlon, Ylat, Ttime, dRdownlw_sfcAll, Xlonf, Ylatf, Ttimef);
    dts_sfc = interp3(Xlon, Ylat, Ttime, dts_sfc, Xlonf, Ylatf, Ttimef);
    dRnetsw_sfcAll = interp3(Xlon, Ylat, Ttime, dRnetsw_sfcAll, Xlonf, Ylatf, Ttimef);

    % Deseasonalized
    startmonth = tLin.startmonth{p_1};
    [dRdownlw_sfcAll, ClimdRdownlw_sfcAll] = monthlyAnomaly3D(144, 72, time, dRdownlw_sfcAll, startmonth);
    [dts_sfc, Climdts_sfc] = monthlyAnomaly3D(144, 72, time, dts_sfc, startmonth);
    [dRnetsw_sfcAll, ClimdRnetsw_sfcAll] = monthlyAnomaly3D(144, 72, time, dRnetsw_sfcAll, startmonth);
    month = 1:12;
    % cal the deltaE(4*epis*T^3^delatT)
    epis=1;% emissivity
    pnum=startmonth-2;
    dDeltaE=zeros(nlonf,nlatf,ntime);
    for ii = 1:length(dts_sfc)
        dDeltaE(:,:,ii)=4*epis*sigma*squeeze(Climdts_sfc(:,:,mod(ii+pnum,12)+1)).^3.*dts_sfc(:,:,ii);
    end
    outpath1 = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi/radEffect/');
    save([outpath1, 'dDeltaE.mat'], 'dDeltaE', 'lonPlot', 'latPlot', 'time', 'readme')
    
    [trendm_dDeltaE, trends_dDeltaE, trendyr_dDeltaE, ~, ~] = autoCalTrend(dDeltaE, nlonf, nlatf, time, startmonth);
    outpath2 = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi/radEffect_trend/');

    save([outpath2, 'trend_dDeltaE.mat'], 'lonPlot', 'latPlot', 'time', 'readme',...
    'trendm_dDeltaE', 'trends_dDeltaE', 'trendyr_dDeltaE')


    %% Save to nc
    saveto = strcat(outpath, 'dts_dradAll_ERAi.nc');
    ncid = netcdf.create(saveto, 'NC_WRITE');

    %Define the dimensions
    dimidlon = netcdf.defDim(ncid, 'longitude', nlonf);
    dimidlat = netcdf.defDim(ncid, 'latitude', nlatf);
    dimidtime = netcdf.defDim(ncid, 'time', ntime);
    dimidmonth = netcdf.defDim(ncid, 'month', 12);
    %Define IDs for the dimension variables
    longitude_ID = netcdf.defVar(ncid, 'longitude', 'double', dimidlon);
    latitude_ID = netcdf.defVar(ncid, 'latitude', 'double', dimidlat);
    time_ID = netcdf.defVar(ncid, 'time', 'double', dimidtime);
    month_ID = netcdf.defVar(ncid, 'month', 'double', dimidmonth);
    %Define units for the dimension variables
    netcdf.putAtt(ncid, longitude_ID, 'standard_name', 'longitude');
    netcdf.putAtt(ncid, latitude_ID, 'standard_name', 'latitude');
    netcdf.putAtt(ncid, time_ID, 'long_name', 'time');
    netcdf.putAtt(ncid, time_ID, 'units', 'days since 0000-00-00 00:00:0.0');
    netcdf.putAtt(ncid, time_ID, 'calendar', 'gregorian');
    netcdf.putAtt(ncid, month_ID, 'long_name', 'month');
    netcdf.putAtt(ncid, month_ID, 'units', 'Month (JFMAMJJASOND)');

    %Define the main variable ()
    dR_downlw_sfc_cld_ID = netcdf.defVar(ncid, 'dRdownlw_sfcAll', 'double', [dimidlon dimidlat dimidtime]);
    dTs_sfc_ID = netcdf.defVar(ncid, 'dts_sfc', 'double', [dimidlon dimidlat dimidtime]);
    dR_netsw_sfc_cld_ID = netcdf.defVar(ncid, 'dRnetsw_sfcAll', 'double', [dimidlon dimidlat dimidtime]);

    ClimdR_downlw_sfc_cld_ID = netcdf.defVar(ncid, 'ClimdRdownlw_sfcAll', 'double', [dimidlon dimidlat dimidmonth]);
    Climdts_sfc_ID = netcdf.defVar(ncid, 'Climdts_sfc', 'double', [dimidlon dimidlat dimidmonth]);
    ClimdRnetsw_sfcAll_ID = netcdf.defVar(ncid, 'ClimdRnetsw_sfcAll', 'double', [dimidlon dimidlat dimidmonth]);
    %Define units for the main variables
    netcdf.putAtt(ncid, dR_downlw_sfc_cld_ID, 'long_name', 'Radiation Anomaly (downwards lw) at SFC');
    netcdf.putAtt(ncid, dR_downlw_sfc_cld_ID, 'units', 'W m**-2');
    netcdf.putAtt(ncid, dTs_sfc_ID, 'long_name', 'Skin temperature');
    netcdf.putAtt(ncid, dTs_sfc_ID, 'units', 'K');
    netcdf.putAtt(ncid, dR_netsw_sfc_cld_ID, 'long_name', 'Radiation Anomaly (net sw) at SFC');
    netcdf.putAtt(ncid, dR_netsw_sfc_cld_ID, 'units', 'W m**-2');

    netcdf.putAtt(ncid, ClimdR_downlw_sfc_cld_ID, 'long_name', 'Radiation Climatology (downwards lw) at SFC');
    netcdf.putAtt(ncid, ClimdR_downlw_sfc_cld_ID, 'units', 'W m**-2');
    netcdf.putAtt(ncid, Climdts_sfc_ID, 'long_name', 'Ts Climatology at SFC');
    netcdf.putAtt(ncid, Climdts_sfc_ID, 'units', 'K');
    netcdf.putAtt(ncid, ClimdRnetsw_sfcAll_ID, 'long_name', 'Radiation Climatology (net sw) at SFC');
    netcdf.putAtt(ncid, ClimdRnetsw_sfcAll_ID, 'units', 'W m**-2');
    %We are done defining the NetCdf
    netcdf.endDef(ncid);

    %Then store the dimension variables in
    netcdf.putVar(ncid, longitude_ID, lonPlot);
    netcdf.putVar(ncid, latitude_ID, latPlot);
    netcdf.putVar(ncid, time_ID, time);
    netcdf.putVar(ncid, month_ID, month);

    %Then store my main variable
    netcdf.putVar(ncid, dR_downlw_sfc_cld_ID, dRdownlw_sfcAll);
    netcdf.putVar(ncid, dTs_sfc_ID, dts_sfc);
    netcdf.putVar(ncid, dR_netsw_sfc_cld_ID, dRnetsw_sfcAll);

    netcdf.putVar(ncid, ClimdR_downlw_sfc_cld_ID, ClimdRdownlw_sfcAll);
    netcdf.putVar(ncid, Climdts_sfc_ID, Climdts_sfc);
    netcdf.putVar(ncid, ClimdRnetsw_sfcAll_ID, ClimdRnetsw_sfcAll);
    %We're done, close the netcdf
    netcdf.close(ncid)

end

t = toc; disp(t)
