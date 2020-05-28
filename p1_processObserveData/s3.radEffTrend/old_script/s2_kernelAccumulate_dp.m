%kernel_accumulate_lyc.m
% this program is aimed to calculate every level's accumulate quantitative value of kernels(nodp)
% ------------------------previous tips--------------------------------
% in the file dp.nc
% player denotes the middle pressure of each layer we perturbed
% plevel denotes the boundary pressure of each layer
% dp denotes the thickness (in hPa) of each layer
clc; clear; tic;
VECTOR_VAR_STR = {'toa','sfc'};
VECTOR_VAR_STR2 = {'all','clr'};
lonf = ncread('/data1/liuyincheng/y_kernels/kernels_YiH/toa/RRTMG_wv_lw_toa_cld_highR.nc', 'lon'); nlonf = length(lonf);
latf = ncread('/data1/liuyincheng/y_kernels/kernels_YiH/toa/RRTMG_wv_lw_toa_cld_highR.nc', 'lat'); nlatf = length(latf);

% modify path first
for ii = 1:2
    for jj = 1:2
        var = char(VECTOR_VAR_STR{ii});
        var2 = char(VECTOR_VAR_STR2{jj});

        if strcmp(var, 'toa')
            % inputpath = 'E:\xieyan\kernels_YiH\toa\';
            inputpath='/data1/liuyincheng/y_kernels/kernels_YiH/toa/';
            dp = ncread([inputpath, 'dp.nc'], 'dp');
            numdp = length(dp);

            if strcmp(var2, 'all')
                wv_lwkernel = ncread([inputpath, 'RRTMG_wv_lw_toa_cld_highR.nc'], 'lwkernel');
                wv_swkernel = ncread([inputpath, 'RRTMG_wv_sw_toa_cld_highR.nc'], 'swkernel');
                t_lwkernel = ncread([inputpath, 'RRTMG_t_toa_cld_highR.nc'], 'lwkernel');
                ts_lwkernel = ncread([inputpath, 'RRTMG_ts_toa_cld_highR.nc'], 'lwkernel');
                alb_swkernel = ncread([inputpath, 'RRTMG_alb_toa_cld_highR.nc'], 'swkernel');
            else
                wv_lwkernel = ncread([inputpath, 'RRTMG_wv_lw_toa_clr_highR.nc'], 'lwkernel');
                wv_swkernel = ncread([inputpath, 'RRTMG_wv_sw_toa_clr_highR.nc'], 'swkernel');
                t_lwkernel = ncread([inputpath, 'RRTMG_t_toa_clr_highR.nc'], 'lwkernel');
                ts_lwkernel = ncread([inputpath, 'RRTMG_ts_toa_clr_highR.nc'], 'lwkernel');
                alb_swkernel = ncread([inputpath, 'RRTMG_alb_toa_clr_highR.nc'], 'swkernel');
            end
            
            wv_lwkernel(isnan(wv_lwkernel)) = 0;
            wv_swkernel(isnan(wv_swkernel)) = 0;
            t_lwkernel(isnan(t_lwkernel)) = 0;
            ts_lwkernel(isnan(ts_lwkernel)) = 0;
            alb_swkernel(isnan(alb_swkernel)) = 0;
            
            
            dp1_path = '/data1/liuyincheng/y_kernels/Ps/';
            dp1 = ncread([dp1_path, 'ERAi_dp1_month.nc'], 'dp1'); %read surface pressure
            dps = ncread([dp1_path, 'ERAi_dp1_month.nc'], 'dps');
            % handle and unified unit: W/m2
            for i = 1:numdp
                wv_lwkernel(:, :, i, :) = wv_lwkernel(:, :, i, :) .* dps(:,:,i,:) / 100;
                wv_swkernel(:, :, i, :) = wv_swkernel(:, :, i, :) .* dps(:,:,i,:)  / 100;
                t_lwkernel(:, :, i , :) = t_lwkernel(:, :, i , :) .* dps(:,:,i,:)  / 100;
            end

        elseif strcmp(var, 'sfc')
            %note that this has an extra dimension which holds upward-downward-net
            %fluxes (we should probs use net)
            % inputpath = 'E:\xieyan\kernels_YiH\surface\';
            inputpath='/data1/liuyincheng/y_kernels/kernels_YiH/surface/';
            dp = ncread([inputpath, 'dp.nc'], 'dp');
            numdp = length(dp);

            if strcmp(var2, 'all')
                wv_lwkernel = ncread([inputpath, 'RRTMG_wv_lw_surface_cld_highR.nc'], 'lwkernel');
                wv_swkernel = ncread([inputpath, 'RRTMG_wv_sw_surface_cld_highR.nc'], 'swkernel');
                t_lwkernel = ncread([inputpath, 'RRTMG_t_surface_cld_highR.nc'], 'lwkernel');
                ts_lwkernel = ncread([inputpath, 'RRTMG_ts_surface_cld_highR.nc'], 'lwkernel');
                alb_swkernel = ncread([inputpath, 'RRTMG_alb_surface_cld_highR.nc'], 'swkernel');
            else
                wv_lwkernel = ncread([inputpath, 'RRTMG_wv_lw_surface_clr_highR.nc'], 'lwkernel');
                wv_swkernel = ncread([inputpath, 'RRTMG_wv_sw_surface_clr_highR.nc'], 'swkernel');
                t_lwkernel = ncread([inputpath, 'RRTMG_t_surface_clr_highR.nc'], 'lwkernel');
                ts_lwkernel = ncread([inputpath, 'RRTMG_ts_surface_clr_highR.nc'], 'lwkernel');
                alb_swkernel = ncread([inputpath, 'RRTMG_alb_surface_clr_highR.nc'], 'swkernel');
            end
            
            wv_lwkernel(isnan(wv_lwkernel)) = 0;
            wv_swkernel(isnan(wv_swkernel)) = 0;
            t_lwkernel(isnan(t_lwkernel)) = 0;
            ts_lwkernel(isnan(ts_lwkernel)) = 0;
            alb_swkernel(isnan(alb_swkernel)) = 0;

            dp1_path = '/data1/liuyincheng/Observe-process/ERAi/kernelsCal/';
            dp1 = ncread([dp1_path, 'ERAi_dp1_month.nc'], 'dp1'); %read surface pressure
            dps = ncread([dp1_path, 'ERAi_dp1_month.nc'], 'dps');
            % handle and unified unit: W/m2, note variable T is special
                t_lwkernel(:, :, 1, :) = squeeze(t_lwkernel(:, :, 1, :)) .* dp1/ 100;

            for i = 1:numdp
                wv_lwkernel(:, :, i, :) = wv_lwkernel(:, :, i, :) .* dps(:,:,i,:) / 100;
                wv_swkernel(:, :, i, :) = wv_swkernel(:, :, i, :) .* dps(:,:,i,:)  / 100;
                t_lwkernel(:, :, i + 1, :) = t_lwkernel(:, :, i + 1, :) .* dps(:,:,i,:)  / 100;
            end

        end

        %% create a .nc file from the regridded data (reference: p-martineau.com/saving-netcdf-in-matlab/)
        outpath = '/data1/liuyincheng/Observe-process/ERAi/kernelsCal/';
        % outpath = 'E:\Repeat_the_experiment\testdata\';
        %saveto = strcat(outpath,'eraiRadAnom_1988',upper(var),var2,'sky.nc');
        saveto = strcat(outpath, 'kernels_accumulate_', upper(var), '_', var2, 'sky_r.nc');
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
