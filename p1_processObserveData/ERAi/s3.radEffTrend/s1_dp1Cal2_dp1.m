%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% mothly data
% cal ps_m and dp1, dp
% month mean ps: ps_m(144x73x12)
% name is consist with HuanYi
% dp1 denotes the thickness (in hPa) of surface layer(month mean): dp1((144x73x12))
% dps denotes the thickness (in hPa) of each layer(month mean): dps((144x73x24x12))
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clc; clear;
% calculate the surface thickness dp1
%modify path first
inputpath = '/data1/liuyincheng/y_kernels/kernels_YiH/surface/';

%% different time series, 1mean 2000-03 to 2018-02(18*12). 2 mean 200207-201706(15*12)
for p_1 = 1:2
    [readme, tLin] = observeParameters(p_1); % fuction% readme,  tLin,
    outpath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi/kernelsCal/');
    auto_mkdir(outpath)

    % lonf = ncread('E:\Repeat_the_experiment\testdata\kernels\kernels_YiH\toa\RRTMG_wv_lw_toa_cld_highR.nc', 'lon'); nlonf = length(lonf);
    % latf = ncread('E:\Repeat_the_experiment\testdata\kernels\kernels_YiH\toa\RRTMG_wv_lw_toa_cld_highR.nc', 'lat'); nlatf = length(latf);
    lonf = ncread('/data1/liuyincheng/y_kernels/kernels_YiH/toa/RRTMG_wv_lw_toa_cld_highR.nc', 'lon'); nlonf = length(lonf);
    latf = ncread('/data1/liuyincheng/y_kernels/kernels_YiH/toa/RRTMG_wv_lw_toa_cld_highR.nc', 'lat'); nlatf = length(latf);

    plevel = ncread([inputpath, 'dp.nc'], 'plevel');
    player = ncread([inputpath, 'dp.nc'], 'player');
    dp_raw = ncread([inputpath, 'dp.nc'], 'dp');
    dps = zeros(nlonf, nlatf, 24, 12); %revise dps in all layers, under the lowest layer are zeros

    p1_path = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi/kernelsCal/');
    ps_m = ncread([p1_path, 'ERAi_ps_month.nc'], 'Ps'); % read surface pressure
    % dp_bottom = ncread([p1_path, 'dp_bottom.nc'], 'dp_bottom');
    dp1 = zeros(nlonf, nlatf, 12);

    for i = 1:nlonf
        for j = 1:nlatf
            for nt = 1:12
                temp = find(player < ps_m(i, j, nt), 1, 'first');
                dp1(i, j, nt) = ps_m(i, j, nt) - plevel(temp + 1);
                % dp1(i, j, nt) = dp_bottom(i,j)123;
                dps(i, j, temp, nt) = dp1(i, j, nt);
                dps(i, j, temp + 1:24, nt) = dp_raw(temp + 1:24);
            end
        end
    end

    %% create a .nc file from the regridded data (reference: p-martineau.com/saving-netcdf-in-matlab/)
    saveto = strcat(outpath, 'ERAi_dp1_month.nc');
    ncid = netcdf.create(saveto, 'NC_WRITE');

    %Define the dimensions
    dimidlon = netcdf.defDim(ncid, 'longitude', nlonf);
    dimidlat = netcdf.defDim(ncid, 'latitude', nlatf);
    dimidlevel = netcdf.defDim(ncid, 'level', 24);
    dimidlayer = netcdf.defDim(ncid, 'layer', 25);
    dimidmonth = netcdf.defDim(ncid, 'month', 12);

    %Define IDs for the dimension variables
    longitude_ID = netcdf.defVar(ncid, 'longitude', 'double', dimidlon);
    latitude_ID = netcdf.defVar(ncid, 'latitude', 'double', dimidlat);
    level_ID = netcdf.defVar(ncid, 'level', 'double', dimidlevel);
    month_ID = netcdf.defVar(ncid, 'month', 'double', dimidmonth);

    %Define units for the dimension variables
    netcdf.putAtt(ncid, longitude_ID, 'standard_name', 'longitude');
    netcdf.putAtt(ncid, latitude_ID, 'standard_name', 'latitude');

    %Define the main variable
    dps_ID = netcdf.defVar(ncid, 'dp1', 'double', [dimidlon dimidlat dimidmonth]);
    dp_ID = netcdf.defVar(ncid, 'dps', 'double', [dimidlon dimidlat dimidlevel dimidmonth]);

    % netcdf.putAtt(ncid,pressure_ID,'standard_name','pressure');
    netcdf.putAtt(ncid, dps_ID, 'long_name', 'surface thickness dp1');
    netcdf.putAtt(ncid, dps_ID, 'units', 'hPa');

    netcdf.putAtt(ncid, dp_ID, 'long_name', 'revise dp in all layers, under the lowest layer are zeros');
    netcdf.putAtt(ncid, dp_ID, 'units', 'hPa');

    %We are done defining the NetCdf
    netcdf.endDef(ncid);

    %Then store the dimension variables in
    netcdf.putVar(ncid, longitude_ID, lonf);
    netcdf.putVar(ncid, latitude_ID, latf);

    %Then store my main variable
    netcdf.putVar(ncid, dps_ID, dp1);
    netcdf.putVar(ncid, dp_ID, dps);

    %We're done, close the netcdf
    netcdf.close(ncid)
end
