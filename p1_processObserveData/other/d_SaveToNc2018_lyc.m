%% Radiation Data
clear

vector_var_str{1} = 'sfc';
vector_var_str{2} = 'toa';
vector_var_str2{1} = 'clr';
vector_var_str2{2} = 'all';

lon = 0:2.5:357.5; nlon = length(lon);
lat = 88.75:-2.5:-88.75; nlat = length(lat);
%%
for ii = 1:2
    for jj = 1:2
        %% Read in Radiation Anomalies
        var = char(vector_var_str{ii});
        var2 = char(vector_var_str2{jj});
        
        filename_vari = strcat('/home/lyc/repeat_Kernel_experiment/testdata/ERAi/eraiRadAnom_2018',upper(var),var2,'sky_r.nc');
        filename_ceres = strcat('/home/lyc/repeat_Kernel_experiment/testdata/CERES_EBAF/ceresAnom',upper(var),'.nc');
%         filename_erai = '/data1/xieyan/ERAi/Rad_ori/ERAi_radAnom.nc';
        filename_ceres1 = strcat('/home/lyc/repeat_Kernel_experiment/testdata/CERES_EBAF/ceresAnomclr',upper(var),'.nc');  % EBAF for clear-sky
        
        if strcmp(var2,'all')
            dR = ncread(filename_ceres,'dR');
            dR_lw = ncread(filename_ceres, 'dR_lw');
            dR_sw = ncread(filename_ceres, 'dR_sw');
            Clim_dR_lw = ncread(filename_ceres,'Clim_dR_lw');
            Clim_dR_sw = ncread(filename_ceres,'Clim_dR_sw');
            Clim_dR = ncread(filename_ceres,'Clim_dR');
        elseif strcmp(var2,'clr')
%             varName_lw = strcat('dR_lw_',var,'_clr');
%             varName_sw = strcat('dR_sw_',var,'_clr');
%             varClim_lw = strcat('ClimdR_lw_',var,'_clr');
%             varClim_sw = strcat('ClimdR_sw_',var,'_clr');
%             dR_lw = ncread(filename_erai,varName_lw);
%             dR_sw = ncread(filename_erai,varName_sw);
%             dR = dR_lw + dR_sw;
%             Clim_dR_lw = ncread(filename_erai,varClim_lw);
%             Clim_dR_sw = ncread(filename_erai,varClim_sw);
%             Clim_dR = Clim_dR_lw+Clim_dR_sw;
            dR = ncread(filename_ceres1,'dR');
            dR_lw = ncread(filename_ceres1, 'dR_lw');
            dR_sw = ncread(filename_ceres1, 'dR_sw');
            Clim_dR_lw = ncread(filename_ceres1,'Clim_dR_lw');
            Clim_dR_sw = ncread(filename_ceres1,'Clim_dR_sw');
            Clim_dR = ncread(filename_ceres1,'Clim_dR');
        end
  
        
        dR_ta = ncread(filename_vari,'tRadEff');
        dR_ts = ncread(filename_vari,'tsRadEff');
        dR_T = dR_ta + dR_ts;
        dR_lw_q = ncread(filename_vari,'wvlwRadEff');
        dR_sw_q = ncread(filename_vari,'wvswRadEff');
        dR_q = dR_lw_q + dR_sw_q;
        dR_alb = ncread(filename_vari,'albRadEff');
        dR_noncloud = ncread(filename_vari,'totalRadEff');
        dR_lw_noncloud = dR_lw_q + dR_T;
        dR_sw_noncloud = dR_sw_q + dR_alb;
        if strcmp(var2,'clr')
            dR_res = dR - dR_noncloud;
            dR_lw_res = dR_lw - dR_lw_noncloud;
            dR_sw_res = dR_sw - dR_sw_noncloud;
        elseif strcmp(var2,'all')
            dR_c = dR - dR_noncloud - dR_res;            
            dR_lw_c = dR_lw - dR_lw_noncloud - dR_lw_res;
            dR_sw_c = dR_sw - dR_sw_noncloud - dR_sw_res;
        end 
        time = ncread(filename_vari,'time'); 
        ntime = length(time);
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        %% save to nc
%         filename = strcat('/data1/xieyan/',var,'_',var2,'.nc');
        filename = strcat('/home/lyc/repeat_Kernel_experiment/testdata/EBAF_clear/',var,'_',var2,'.nc');
        ncid = netcdf.create(filename,'NC_WRITE');
        
        %Define the dimensions of the variable.
        time_dimID = netcdf.defDim(ncid,'time',ntime);
        lat_dimID = netcdf.defDim(ncid,'latitude',nlat);
        lon_dimID = netcdf.defDim(ncid,'longitude',nlon);
        month_dimID = netcdf.defDim(ncid,'month',12);
        %Define IDs for the dimension variables
        time_ID = netcdf.defVar(ncid,'time','double',time_dimID);
        lat_ID = netcdf.defVar(ncid,'latitude','double',lat_dimID);
        lon_ID = netcdf.defVar(ncid,'longitude','double',lon_dimID);
        month_ID = netcdf.defVar(ncid,'month','double',month_dimID);
        %Define a new variable in the file.
        dR_varID = netcdf.defVar(ncid,'dR','double',[lon_dimID,lat_dimID,time_dimID]);
        dR_lw_varID = netcdf.defVar(ncid,'dR_lw','double',[lon_dimID,lat_dimID,time_dimID]);
        dR_sw_varID = netcdf.defVar(ncid,'dR_sw','double',[lon_dimID,lat_dimID,time_dimID]);
        dR_T_varID = netcdf.defVar(ncid,'dR_T','double',[lon_dimID,lat_dimID,time_dimID]);
        dR_ta_varID = netcdf.defVar(ncid,'dR_ta','double',[lon_dimID,lat_dimID,time_dimID]);
        dR_ts_varID = netcdf.defVar(ncid,'dR_ts','double',[lon_dimID,lat_dimID,time_dimID]);
        dR_q_varID = netcdf.defVar(ncid,'dR_q','double',[lon_dimID,lat_dimID,time_dimID]);
        dR_lw_q_varID = netcdf.defVar(ncid,'dR_lw_q','double',[lon_dimID,lat_dimID,time_dimID]);
        dR_sw_q_varID = netcdf.defVar(ncid,'dR_sw_q','double',[lon_dimID,lat_dimID,time_dimID]);
        dR_alb_varID = netcdf.defVar(ncid,'dR_alb','double',[lon_dimID,lat_dimID,time_dimID]);
        dR_noncloud_varID = netcdf.defVar(ncid,'dR_noncloud','double',[lon_dimID,lat_dimID,time_dimID]);
        dR_lw_noncloud_varID = netcdf.defVar(ncid,'dR_lw_noncloud','double',[lon_dimID,lat_dimID,time_dimID]);
        dR_sw_noncloud_varID = netcdf.defVar(ncid,'dR_sw_noncloud','double',[lon_dimID,lat_dimID,time_dimID]);
        dR_res_varID = netcdf.defVar(ncid,'dR_res','double',[lon_dimID,lat_dimID,time_dimID]);
        dR_lw_res_varID = netcdf.defVar(ncid,'dR_lw_res','double',[lon_dimID,lat_dimID,time_dimID]);
        dR_sw_res_varID = netcdf.defVar(ncid,'dR_sw_res','double',[lon_dimID,lat_dimID,time_dimID]);
        Clim_dR_lw_varID = netcdf.defVar(ncid,'Clim_dR_lw','double',[lon_dimID,lat_dimID,month_dimID]);
        Clim_dR_sw_varID = netcdf.defVar(ncid,'Clim_dR_sw','double',[lon_dimID,lat_dimID,month_dimID]);
        Clim_dR_varID = netcdf.defVar(ncid,'Clim_dR','double',[lon_dimID,lat_dimID,month_dimID]);
        
        if strcmp(var2,'all')
            dR_c_varID = netcdf.defVar(ncid,'dR_c','double',[lon_dimID,lat_dimID,time_dimID]);
            dR_lw_c_varID = netcdf.defVar(ncid,'dR_lw_c','double',[lon_dimID,lat_dimID,time_dimID]);
            dR_sw_c_varID = netcdf.defVar(ncid,'dR_sw_c','double',[lon_dimID,lat_dimID,time_dimID]);
        end
        
        
        %Define units for the dimension variables
        netcdf.putAtt(ncid,lon_ID,'standard_name','longitude');
        netcdf.putAtt(ncid,lat_ID,'standard_name','latitude');
        
        % netcdf.putAtt(ncid,pressure_ID,'standard_name','pressure');
        netcdf.putAtt(ncid,time_ID,'long_name','time');
        netcdf.putAtt(ncid,time_ID,'units','days since 0000-00-00 00:00:0.0');
        netcdf.putAtt(ncid,time_ID,'calendar','gregorian');
        netcdf.putAtt(ncid,month_ID,'long_name','month');
        netcdf.putAtt(ncid,month_ID,'units','JFMAMJJASOND');
        
        netcdf.putAtt(ncid, dR_varID,'long_name','Total Net Radiation Anomaly');
        netcdf.putAtt(ncid, dR_varID,'units', 'W m-2');
        netcdf.putAtt(ncid, dR_lw_varID,'long_name','Total LW Radiation Anomaly');
        netcdf.putAtt(ncid, dR_lw_varID,'units', 'W m-2');
        netcdf.putAtt(ncid, dR_sw_varID,'long_name','Total SW Radiation Anomaly');
        netcdf.putAtt(ncid, dR_sw_varID,'units', 'W m-2');
        netcdf.putAtt(ncid, dR_T_varID,'long_name','Total Radiation Anomaly due to Effect of Temperature');
        netcdf.putAtt(ncid, dR_T_varID,'units', 'W m-2');
        netcdf.putAtt(ncid, dR_ta_varID,'long_name','Radiation Anomaly due to Effect of Temperature (in the atmosphere)');
        netcdf.putAtt(ncid, dR_ta_varID,'units', 'W m-2');
        netcdf.putAtt(ncid, dR_ts_varID,'long_name','Radiation Anomaly due to Effect of Temperature (at surface)');
        netcdf.putAtt(ncid, dR_ts_varID,'units', 'W m-2');
        netcdf.putAtt(ncid, dR_q_varID,'long_name','Net Radiation Anomaly due to Effect of Water Vapor');
        netcdf.putAtt(ncid, dR_q_varID,'units', 'W m-2');
        netcdf.putAtt(ncid, dR_lw_q_varID,'long_name','LW Radiation Anomaly due to Effect of Water Vapor');
        netcdf.putAtt(ncid, dR_lw_q_varID,'units', 'W m-2');
        netcdf.putAtt(ncid, dR_sw_q_varID,'long_name','SW Radiation Anomaly due to Effect of Water Vapor');
        netcdf.putAtt(ncid, dR_sw_q_varID,'units', 'W m-2');
        netcdf.putAtt(ncid, dR_alb_varID,'long_name','Radiation Anomaly due to Effect of Albedo');
        netcdf.putAtt(ncid, dR_alb_varID,'units', 'W m-2');
        netcdf.putAtt(ncid, dR_noncloud_varID,'long_name','Total Net Noncloud Radiation Anomaly');
        netcdf.putAtt(ncid, dR_noncloud_varID,'units', 'W m-2');
        netcdf.putAtt(ncid, dR_lw_noncloud_varID,'long_name','Total LW Noncloud Radiation Anomaly');
        netcdf.putAtt(ncid, dR_lw_noncloud_varID,'units', 'W m-2');
        netcdf.putAtt(ncid, dR_sw_noncloud_varID,'long_name','Total SW Noncloud Radiation Anomaly');
        netcdf.putAtt(ncid, dR_sw_noncloud_varID,'units', 'W m-2');
        netcdf.putAtt(ncid, dR_res_varID,'long_name','Net Residual Radiation Anomaly');
        netcdf.putAtt(ncid, dR_res_varID,'units', 'W m-2');
        netcdf.putAtt(ncid, dR_lw_res_varID,'long_name','LW Residual Radiation Anomaly');
        netcdf.putAtt(ncid, dR_lw_res_varID,'units', 'W m-2');
        netcdf.putAtt(ncid, dR_sw_res_varID,'long_name','SW Residual Radiation Anomaly');
        netcdf.putAtt(ncid, dR_sw_res_varID,'units', 'W m-2');
        netcdf.putAtt(ncid, Clim_dR_lw_varID,'long_name','LW Radiation Climatology');
        netcdf.putAtt(ncid, Clim_dR_lw_varID,'units', 'W m-2');
        netcdf.putAtt(ncid, Clim_dR_sw_varID,'long_name','SW Radiation Climatology');
        netcdf.putAtt(ncid, Clim_dR_sw_varID,'units', 'W m-2');
        netcdf.putAtt(ncid, Clim_dR_varID,'long_name','NET Radiation Climatology');
        netcdf.putAtt(ncid, Clim_dR_varID,'units', 'W m-2');
        
        if strcmp(var2,'all')
            netcdf.putAtt(ncid, dR_c_varID,'long_name','Net Radiation Anomaly due to Effect of Cloud');
            netcdf.putAtt(ncid, dR_c_varID,'units', 'W m-2');
            netcdf.putAtt(ncid, dR_lw_c_varID,'long_name','LW Radiation Anomaly due to Effect of Cloud');
            netcdf.putAtt(ncid, dR_lw_c_varID,'units', 'W m-2');
            netcdf.putAtt(ncid, dR_sw_c_varID,'long_name','SW Radiation Anomaly due to Effect of Cloud');
            netcdf.putAtt(ncid, dR_sw_c_varID,'units', 'W m-2');
        end
        
        %Leave define mode and enter data mode to write data.
        netcdf.endDef(ncid);
        
        
        %Write data to variable.
        netcdf.putVar(ncid,time_ID,time);
        netcdf.putVar(ncid,month_ID,1:12);
        netcdf.putVar(ncid,lat_ID,lat);
        netcdf.putVar(ncid,lon_ID,lon);
        
        %Write data to variable.
        netcdf.putVar(ncid,dR_varID,dR);
        netcdf.putVar(ncid,dR_lw_varID,dR_lw);
        netcdf.putVar(ncid,dR_sw_varID,dR_sw);
        netcdf.putVar(ncid,dR_T_varID,dR_T);
        netcdf.putVar(ncid,dR_ta_varID,dR_ta);
        netcdf.putVar(ncid,dR_ts_varID,dR_ts);
        netcdf.putVar(ncid,dR_q_varID,dR_q);
        netcdf.putVar(ncid,dR_lw_q_varID,dR_lw_q);
        netcdf.putVar(ncid,dR_sw_q_varID,dR_sw_q);
        netcdf.putVar(ncid,dR_alb_varID,dR_alb);
        netcdf.putVar(ncid,dR_noncloud_varID,dR_noncloud);
        netcdf.putVar(ncid,dR_lw_noncloud_varID,dR_lw_noncloud);
        netcdf.putVar(ncid,dR_sw_noncloud_varID,dR_sw_noncloud);
        netcdf.putVar(ncid,Clim_dR_lw_varID,Clim_dR_lw);
        netcdf.putVar(ncid,Clim_dR_sw_varID,Clim_dR_sw);
        netcdf.putVar(ncid,Clim_dR_varID,Clim_dR);
        netcdf.putVar(ncid,dR_res_varID,dR_res);
        netcdf.putVar(ncid,dR_lw_res_varID,dR_lw_res);
        netcdf.putVar(ncid,dR_sw_res_varID,dR_sw_res);
        

        if strcmp(var2,'all')
            netcdf.putVar(ncid,dR_c_varID,dR_c);
            netcdf.putVar(ncid,dR_lw_c_varID,dR_lw_c);
            netcdf.putVar(ncid,dR_sw_c_varID,dR_sw_c);
        end
        
        %Close the file.
        netcdf.close(ncid)
    end
end
        

        