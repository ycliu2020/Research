%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% mothly data
% cal ps_m and dps, dp
% month mean ps: ps_m(144x73x12)
% dps denotes the thickness (in hPa) of surface layer(month mean): dps((144x73x12))
% dp denotes the thickness (in hPa) of each layer(month mean): dp((144x73x12))
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clc; clear;
% calculate the surface thickness dps
inputpath = '/data1/liuyincheng/y_kernels/kernels_YiH/surface/';
outpath = '/data1/liuyincheng/Observe-process/ERAi/kernelsCal/';
% outpath = 'E:\Repeat_the_experiment\testdata\';

% lonf = ncread('E:\Repeat_the_experiment\testdata\kernels\kernels_YiH\toa\RRTMG_wv_lw_toa_cld_highR.nc', 'lon'); nlonf = length(lonf);
% latf = ncread('E:\Repeat_the_experiment\testdata\kernels\kernels_YiH\toa\RRTMG_wv_lw_toa_cld_highR.nc', 'lat'); nlatf = length(latf);
lonf = ncread('/data1/liuyincheng/y_kernels/kernels_YiH/toa/RRTMG_wv_lw_toa_cld_highR.nc', 'lon'); nlonf = length(lonf);
latf = ncread('/data1/liuyincheng/y_kernels/kernels_YiH/toa/RRTMG_wv_lw_toa_cld_highR.nc', 'lat'); nlatf = length(latf);

plevel = ncread([inputpath, 'dp.nc'], 'plevel');
player = ncread([inputpath, 'dp.nc'], 'player');
dp_raw = ncread([inputpath, 'dp.nc'], 'dp');
dp = zeros(nlonf,nlatf,24,12); %revise dp in all layers, under the lowest layer are zeros

% ps_m_path = 'E:\Repeat_the_experiment\testdata\kernels\Ps\';
ps_m_path = '/data1/liuyincheng/Observe-process/ERAi/kernelsCal/';
ps_m = ncread([ps_m_path, 'ERAi_ps_month.nc'], 'Ps'); % read surface pressure
% dp_bottom = ncread([ps_m_path, 'dp_bottom.nc'], 'dp_bottom'); 
dps = zeros(nlonf, nlatf, 12);
for i = 1:nlonf
    for j = 1:nlatf
        for nt = 1:12
            temp = find(player < ps_m(i, j, nt), 1, 'first');
            dps(i, j, nt) = ps_m(i, j, nt) - plevel(temp+1);
            % dps(i, j, nt) = dp_bottom(i,j)123;
            dp(i,j,temp,nt) = dps(i, j, nt);
            dp(i,j,temp+1:24,nt)=dp_raw(temp+1:24);
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
dps_ID = netcdf.defVar(ncid, 'dps', 'double', [dimidlon dimidlat dimidmonth]);
dp_ID = netcdf.defVar(ncid, 'dp', 'double', [dimidlon dimidlat dimidlevel dimidmonth]);

% netcdf.putAtt(ncid,pressure_ID,'standard_name','pressure');
netcdf.putAtt(ncid, dps_ID, 'long_name', 'surface thickness dps');
netcdf.putAtt(ncid, dps_ID, 'units', 'hPa');

netcdf.putAtt(ncid, dp_ID, 'long_name', 'revise dp in all layers, under the lowest layer are zeros');
netcdf.putAtt(ncid, dp_ID, 'units', 'hPa');


%We are done defining the NetCdf
netcdf.endDef(ncid);

%Then store the dimension variables in
netcdf.putVar(ncid, longitude_ID, lonf);
netcdf.putVar(ncid, latitude_ID, latf);

%Then store my main variable
netcdf.putVar(ncid, dps_ID, dps);
netcdf.putVar(ncid, dp_ID, dp);

%We're done, close the netcdf
netcdf.close(ncid)