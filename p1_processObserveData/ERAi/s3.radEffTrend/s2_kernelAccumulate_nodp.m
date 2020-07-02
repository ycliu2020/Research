%kernel_accumulate_nodp.m
% this program is aimed to calculate every level's accumulate quantitative value of kernels(nodp)
% ------------------------previous tips--------------------------------
% in the file dp.nc
% player denotes the middle pressure of each layer we perturbed
% plevel denotes the boundary pressure of each layer
% dp denotes the thickness (in hPa) of each layer
% ------------------------scheme--------------------------------
% to use the kernels:
% 1. read kernels of atmospheric temperature (Kt), water vapor (Kq),surface temperature (Kts) and albedo (Kalb)
%    read player from dp.nc

% 2. read anomaly of Ta (dTa), q (dq), Ts (dTs), and albedo (dalb), and Ta, surface pressure (sp) 
%    then interpolate to the resolution of the kernels

% 3. calculate the feedbacks
% dr_q =  dim_sum_n(Kq*dqn,0),0)
% dr_Ta=  dim_sum_n(Kt*dTa,0),0)(toa);dr_Ta=  dim_sum_n(Kt(1:,:,:)*dTa,0),0)+Kt(0,:,:)*dTs(surface)
% dr_Ts=  Kts*dTs
% dr_alb= Kalb*dalb

clc; clear; tic;
vector_var_str = {'toa','sfc'};
vector_var_str2 = {'all','clr'};
lonf = ncread('/data/pub/kernel-nodp-highR/toa/RRTMG_wv_lw_toa_cld_highR.nc', 'lon'); nlonf = length(lonf);
latf = ncread('/data/pub/kernel-nodp-highR/toa/RRTMG_wv_lw_toa_cld_highR.nc', 'lat'); nlatf = length(latf);
filepath1='/data/pub/kernel-nodp-highR/toa/';
filepath2='/data/pub/kernel-nodp-highR/surface/';
% modify path first
for ii = 1:2
    for jj = 1:2
        var = char(vector_var_str{ii});
        var2 = char(vector_var_str2{jj});

        if strcmp(var, 'toa')
            % filepath = 'E:\xieyan\kernels_YiH\toa\';
            dp = ncread([filepath1, 'dp.nc'], 'dp');
            numdp = length(dp);

            if strcmp(var2, 'all')
                wv_lwkernel = ncread([filepath1, 'RRTMG_wv_lw_toa_cld_nodp_highR.nc'], 'lwkernel');
                wv_swkernel = ncread([filepath1, 'RRTMG_wv_sw_toa_cld_nodp_highR.nc'], 'swkernel');
                t_lwkernel = ncread([filepath1, 'RRTMG_t_toa_cld_nodp_highR.nc'], 'lwkernel');
                ts_lwkernel = ncread([filepath1, 'RRTMG_ts_toa_cld_nodp_highR.nc'], 'lwkernel');
                alb_swkernel = ncread([filepath1, 'RRTMG_alb_toa_cld_nodp_highR.nc'], 'swkernel');
            else
                wv_lwkernel = ncread([filepath1, 'RRTMG_wv_lw_toa_clr_nodp_highR.nc'], 'lwkernel');
                wv_swkernel = ncread([filepath1, 'RRTMG_wv_sw_toa_clr_nodp_highR.nc'], 'swkernel');
                t_lwkernel = ncread([filepath1, 'RRTMG_t_toa_clr_nodp_highR.nc'], 'lwkernel');
                ts_lwkernel = ncread([filepath1, 'RRTMG_ts_toa_clr_nodp_highR.nc'], 'lwkernel');
                alb_swkernel = ncread([filepath1, 'RRTMG_alb_toa_clr_nodp_highR.nc'], 'swkernel');
            end
            
            wv_lwkernel(isnan(wv_lwkernel)) = 0;
            wv_swkernel(isnan(wv_swkernel)) = 0;
            t_lwkernel(isnan(t_lwkernel)) = 0;
            ts_lwkernel(isnan(ts_lwkernel)) = 0;
            alb_swkernel(isnan(alb_swkernel)) = 0;
            
            
        elseif strcmp(var, 'sfc')
            %note that this has an extra dimension which holds upward-downward-net
            %fluxes (we should probs use net)
            % filepath = 'E:\xieyan\kernels_YiH\surface\';

            dp = ncread([filepath2, 'dp.nc'], 'dp');
            numdp = length(dp);

            if strcmp(var2, 'all')
                wv_lwkernel = ncread([filepath2, 'RRTMG_wv_lw_surface_cld_nodp_highR.nc'], 'lwkernel');
                wv_swkernel = ncread([filepath2, 'RRTMG_wv_sw_surface_cld_nodp_highR.nc'], 'swkernel');
                t_lwkernel = ncread([filepath2, 'RRTMG_t_surface_cld_nodp_highR.nc'], 'lwkernel');
                ts_lwkernel = ncread([filepath2, 'RRTMG_ts_surface_cld_nodp_highR.nc'], 'lwkernel');
                alb_swkernel = ncread([filepath2, 'RRTMG_alb_surface_cld_nodp_highR.nc'], 'swkernel');
            else
                wv_lwkernel = ncread([filepath2, 'RRTMG_wv_lw_surface_clr_nodp_highR.nc'], 'lwkernel');
                wv_swkernel = ncread([filepath2, 'RRTMG_wv_sw_surface_clr_nodp_highR.nc'], 'swkernel');
                t_lwkernel = ncread([filepath2, 'RRTMG_t_surface_clr_nodp_highR.nc'], 'lwkernel');
                ts_lwkernel = ncread([filepath2, 'RRTMG_ts_surface_clr_nodp_highR.nc'], 'lwkernel');
                alb_swkernel = ncread([filepath2, 'RRTMG_alb_surface_clr_nodp_highR.nc'], 'swkernel');
            end
            
            wv_lwkernel(isnan(wv_lwkernel)) = 0;
            wv_swkernel(isnan(wv_swkernel)) = 0;
            t_lwkernel(isnan(t_lwkernel)) = 0;
            ts_lwkernel(isnan(ts_lwkernel)) = 0;
            alb_swkernel(isnan(alb_swkernel)) = 0;


        end

        %% create a .nc file from the regridded data (reference: p-martineau.com/saving-netcdf-in-matlab/)
        directory_name = '/home/lyc/repeat_Kernel_experiment/testdata/ERAi/';
        % directory_name = 'E:\Repeat_the_experiment\testdata\';
        %saveto = strcat(directory_name,'eraiRadAnom_1988',upper(var),var2,'sky.nc');
        saveto = strcat(directory_name, 'kernel_accumulate_', upper(var), '_', var2, 'sky_nodp.nc');
        ncid = netcdf.create(saveto, 'NC_WRITE');

        %Define the dimensions
        dimidlon = netcdf.defDim(ncid, 'longitude', nlonf);
        dimidlat = netcdf.defDim(ncid, 'latitude', nlatf);
        dimidlevel = netcdf.defDim(ncid, 'level', 24);
        dimidtlevel = netcdf.defDim(ncid, 'tlevel', size(t_lwkernel, 3));
        dimidmonth = netcdf.defDim(ncid, 'month', 12);

        %Define IDs for the dimension variables
        longitude_ID = netcdf.defVar(ncid, 'longitude', 'double', dimidlon);
        latitude_ID = netcdf.defVar(ncid, 'latitude', 'double', dimidlat);
        level_ID = netcdf.defVar(ncid, 'level', 'double', dimidlevel);
        tlevel_ID = netcdf.defVar(ncid, 'tlevel', 'double', dimidtlevel);
        month_ID = netcdf.defVar(ncid, 'month', 'double', dimidmonth);

        %Define units for the dimension variables
        netcdf.putAtt(ncid, longitude_ID, 'standard_name', 'longitude');
        netcdf.putAtt(ncid, latitude_ID, 'standard_name', 'latitude');

        %Define the main variable
        t_ID = netcdf.defVar(ncid, 't_lwkernel', 'double', [dimidlon dimidlat dimidtlevel dimidmonth]);
        wvlw_ID = netcdf.defVar(ncid, 'wv_lwkernel', 'double', [dimidlon dimidlat dimidlevel dimidmonth]);
        wvsw_ID = netcdf.defVar(ncid, 'wv_swkernel', 'double', [dimidlon dimidlat dimidlevel dimidmonth]);
        ts_ID = netcdf.defVar(ncid, 'ts_lwkernel', 'double', [dimidlon dimidlat dimidmonth]);
        alb_ID = netcdf.defVar(ncid, 'alb_swkernel', 'double', [dimidlon dimidlat dimidmonth]);

        % netcdf.putAtt(ncid,pressure_ID,'standard_name','pressure');

        netcdf.putAtt(ncid, t_ID, 'long_name', 'T kernel (longwave)');
        netcdf.putAtt(ncid, t_ID, 'units', 'W/m2/K');

        netcdf.putAtt(ncid, wvlw_ID, 'long_name', 'water vapor kernel (longwave)');
        netcdf.putAtt(ncid, wvlw_ID, 'units', 'W/m2/K');

        netcdf.putAtt(ncid, wvsw_ID, 'long_name', 'water vapor kernel (shortwave)');
        netcdf.putAtt(ncid, wvsw_ID, 'units', 'W/m2/K');

        netcdf.putAtt(ncid, ts_ID, 'long_name', ' Ts kernel (longwave)');
        netcdf.putAtt(ncid, ts_ID, 'units', 'W/m2/K');

        netcdf.putAtt(ncid, alb_ID, 'long_name', 'alb kernel (shortwave)');
        netcdf.putAtt(ncid, alb_ID, 'units', 'W/m2/0.01');

        %We are done defining the NetCdf
        netcdf.endDef(ncid);

        %Then store the dimension variables in
        netcdf.putVar(ncid, longitude_ID, lonf);
        netcdf.putVar(ncid, latitude_ID, latf);

        %Then store my main variable
        netcdf.putVar(ncid, t_ID, t_lwkernel);
        netcdf.putVar(ncid, wvlw_ID, wv_lwkernel);
        netcdf.putVar(ncid, wvsw_ID, wv_swkernel);
        netcdf.putVar(ncid, ts_ID, ts_lwkernel);
        netcdf.putVar(ncid, alb_ID, alb_swkernel);

        %We're done, close the netcdf
        netcdf.close(ncid)
    end

end
