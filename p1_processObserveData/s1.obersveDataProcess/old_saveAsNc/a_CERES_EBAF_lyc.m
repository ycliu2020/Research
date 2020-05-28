% This program is to calculate the radiation anomalies from CERES/EBAF4.0
clc;clear; tic;pre_load;
vector_var_str{1} = 'sfc';
vector_var_str{2} = 'toa';
%modify path first
filepath = 'H:\xieyan\CERES_EBAF\';
testpath= 'H:\Repeat_the_experiment\testdata\';
filename1 = 'CERES_EBAF-Surface_Ed4.0_Subset_200003-201802.nc';
filename2 = 'CERES_EBAF-TOA_Ed4.0_Subset_200003-201802.nc';

for level = 1:2 % 1=sfc, 2=toa
   var = char(vector_var_str{level});
   % define positive as downward
   if level == 1
       filename = strcat(filepath, filename1);
%  % All-sky
%        dR = ncread(filename,'sfc_net_tot_all_mon');                     
%        dR_lw = ncread(filename,'sfc_net_lw_all_mon');
%        dR_sw = ncread(filename,'sfc_net_sw_all_mon');
%        dR_cld = ncread(filename,'sfc_cre_net_tot_mon');                   % Cloud radiative effect
%        dR_cld_lw = ncread(filename,'sfc_cre_net_lw_mon');
%        dR_cld_sw = ncread(filename,'sfc_cre_net_sw_mon');
%  % Clear-sky
       dR = ncread(filename,'sfc_net_tot_clr_mon');
       dR_lw = ncread(filename,'sfc_net_lw_clr_mon');
       dR_sw = ncread(filename,'sfc_net_sw_clr_mon');
       lon = -0.5:1:360.5; nlon = length(lon); lon = double(lon);
       lat = -89.5:1:89.5; nlat = length(lat); lat = double(lat);
       time = double(ncread(filename,'time'));                            % days since 2000-03-01 00:00:00
   elseif level == 2
       filename = strcat(filepath, filename2);
%        dR = ncread(filename,'toa_net_all_mon');
%        dR_lw = -ncread(filename,'toa_lw_all_mon');                        % define positive as downward
%        dR_sw = -ncread(filename,'toa_sw_all_mon');
%        dR_cld = ncread(filename,'toa_cre_net_mon');                       % Cloud radiative effect
%        dR_cld_lw = ncread(filename,'toa_cre_lw_mon');
%        dR_cld_sw = ncread(filename,'toa_cre_sw_mon');
       dR = ncread(filename,'toa_net_clr_mon');
       dR_lw = -ncread(filename,'toa_lw_clr_mon');                        % define positive as downward
       dR_sw = -ncread(filename,'toa_sw_clr_mon');
       lon = -0.5:1:360.5; nlon = length(lon); lon = double(lon);  
       lat = -89.5:1:89.5; nlat = length(lat); lat = double(lat);
       time = double(ncread(filename,'time'));                                    % days since 2000-03-01 00:00:00
   end
   time = datenum(2000,03,1 + double(time),00,00,00);
   ntime = length(time);
   month = 1:12;
   
   %% Regrid (Because of kernels)
   lonf = 0:2.5:357.5; nlonf = length(lonf);
   latf = 88.75:-2.5:-88.75; nlatf = length(latf);
   % For interpolate(此步作用是将lonf变成lon的子集(lon范围扩大), latf原来就是lat的子集所以不用扩充)
   for ll = 360:-1:1
       dR(ll+1,:,:) = dR(ll,:,:); 
       dR_lw(ll+1,:,:) = dR_lw(ll,:,:); 
       dR_sw(ll+1,:,:) = dR_sw(ll,:,:);
%        dR_cld(ll+1,:,:) = dR_cld(ll,:,:);
%        dR_cld_lw(ll+1,:,:) = dR_cld_lw(ll,:,:);
%        dR_cld_sw(ll+1,:,:) = dR_cld_sw(ll,:,:);
   end
   dR(1,:,:) = dR(361,:,:);
   dR_lw(1,:,:) = dR_lw(361,:,:);
   dR_sw(1,:,:) = dR_sw(361,:,:);
%    dR_cld(1,:,:) = dR_cld(361,:,:);
%    dR_cld_lw(1,:,:) = dR_cld_lw(361,:,:);
%    dR_cld_sw(1,:,:) = dR_cld_sw(361,:,:);
   
   dR(362,:,:) = dR(2,:,:);
   dR_lw(362,:,:) = dR_lw(2,:,:);
   dR_sw(362,:,:) = dR_sw(2,:,:); 
%    dR_cld(362,:,:) = dR_cld(2,:,:);
%    dR_cld_lw(362,:,:) = dR_cld_lw(2,:,:);
%    dR_cld_sw(362,:,:) = dR_cld_sw(2,:,:); 
   
   % Regraded to 2.5*2.5
   time2 = 1:ntime;
   [Xlon, Ylat, Ttime] = meshgrid(lat,lon,time2);
   [Xlonf, Ylatf, Ttimef] = meshgrid(latf,lonf,time2);
   
   dR = interp3(Xlon,Ylat,Ttime,dR,Xlonf,Ylatf,Ttimef);
   dR_lw = interp3(Xlon,Ylat,Ttime,dR_lw,Xlonf,Ylatf,Ttimef);
   dR_sw = interp3(Xlon,Ylat,Ttime,dR_sw,Xlonf,Ylatf,Ttimef);
%    dR_cld = interp3(Xlon,Ylat,Ttime,dR_cld,Xlonf,Ylatf,Ttimef);
%    dR_cld_lw = interp3(Xlon,Ylat,Ttime,dR_cld_lw,Xlonf,Ylatf,Ttimef);
%    dR_cld_sw = interp3(Xlon,Ylat,Ttime,dR_cld_sw,Xlonf,Ylatf,Ttimef);
   
   %% Deseasonalized Anomaly
   [dR,Clim_dR] = monthlyAnomaly3D(nlonf,nlatf,time,dR);
   [dR_lw,Clim_dR_lw] = monthlyAnomaly3D(nlonf,nlatf,time,dR_lw);
   [dR_sw,Clim_dR_sw] = monthlyAnomaly3D(nlonf,nlatf,time,dR_sw);
%    [dR_cld,Clim_dR_cld] = monthlyAnomaly3D(nlonf,nlatf,time,dR_cld);
%    [dR_cld_lw,Clim_dR_cld_lw] = monthlyAnomaly3D(nlonf,nlatf,time,dR_cld_lw);
%    [dR_cld_sw,Clim_dR_cld_sw] = monthlyAnomaly3D(nlonf,nlatf,time,dR_cld_sw);
      
   %% create a .nc file (reference: p-martineau.com/saving-netcdf-in-matlab/)
   saveto = strcat(testpath,'ceresAnomclr',upper(var),'.nc');
   ncid = netcdf.create(saveto,'NC_WRITE');
   %Define the dimensions
   dimidlon = netcdf.defDim(ncid,'longitude',nlonf);
   dimidlat = netcdf.defDim(ncid,'latitude',nlatf);
   dimidtime = netcdf.defDim(ncid,'time',ntime);
   dimidmonth = netcdf.defDim(ncid,'month',12);
   %Define IDs for the dimension variables
   longitude_ID=netcdf.defVar(ncid,'longitude','double',dimidlon);
   latitude_ID=netcdf.defVar(ncid,'latitude','double',dimidlat);
   time_ID=netcdf.defVar(ncid,'time','double',dimidtime);
   month_ID=netcdf.defVar(ncid,'month','double',dimidmonth);
   %Define the main variable
   dR_ID = netcdf.defVar(ncid,'dR','double',[dimidlon dimidlat dimidtime]);
   dR_lw_ID = netcdf.defVar(ncid,'dR_lw','double',[dimidlon dimidlat dimidtime]);
   dR_sw_ID = netcdf.defVar(ncid,'dR_sw','double',[dimidlon dimidlat dimidtime]);
   Clim_dR_ID = netcdf.defVar(ncid,'Clim_dR','double',[dimidlon dimidlat dimidmonth]);
   Clim_dR_lw_ID = netcdf.defVar(ncid,'Clim_dR_lw','double',[dimidlon dimidlat dimidmonth]);
   Clim_dR_sw_ID = netcdf.defVar(ncid,'Clim_dR_sw','double',[dimidlon dimidlat dimidmonth]);
   
%    dR_cld_ID = netcdf.defVar(ncid,'dR_cld','double',[dimidlon dimidlat dimidtime]);
%    dR_cld_lw_ID = netcdf.defVar(ncid,'dR_cld_lw','double',[dimidlon dimidlat dimidtime]);
%    dR_cld_sw_ID = netcdf.defVar(ncid,'dR_cld_sw','double',[dimidlon dimidlat dimidtime]);
%    Clim_dR_cld_ID = netcdf.defVar(ncid,'Clim_dR_cld','double',[dimidlon dimidlat dimidmonth]);
%    Clim_dR_cld_lw_ID = netcdf.defVar(ncid,'Clim_dR_cld_lw','double',[dimidlon dimidlat dimidmonth]);
%    Clim_dR_cld_sw_ID = netcdf.defVar(ncid,'Clim_dR_cld_sw','double',[dimidlon dimidlat dimidmonth]);


   %Define units for the dimension variables
   netcdf.putAtt(ncid,longitude_ID,'standard_name','longitude');
   netcdf.putAtt(ncid,latitude_ID,'standard_name','latitude');
   netcdf.putAtt(ncid,time_ID,'long_name','time');
   netcdf.putAtt(ncid,time_ID,'units','days since 0000-00-00 00:00:0.0');
   netcdf.putAtt(ncid,time_ID,'calendar','gregorian');
   
   netcdf.putAtt(ncid,month_ID,'long_name','month');
   netcdf.putAtt(ncid,month_ID,'units','Month (JFMAMJJASOND)');
   netcdf.putAtt(ncid, dR_ID,'long_name','Total Radiation Anomaly, Monthly Means, All-Sky Conditions');
   netcdf.putAtt(ncid, dR_ID,'units', 'W m^-2');
   netcdf.putAtt(ncid, dR_lw_ID,'long_name','LW Radiation Anomaly, Monthly Means, All-Sky Conditions');
   netcdf.putAtt(ncid, dR_lw_ID,'units', 'W m^-2');
   netcdf.putAtt(ncid, dR_sw_ID,'long_name','SW Radiation Anomaly, Monthly Means, All-Sky Conditions');
   netcdf.putAtt(ncid, dR_sw_ID,'units', 'W m^-2');
   netcdf.putAtt(ncid, Clim_dR_ID,'long_name','Total Radiation Climatology, All-Sky Conditions');
   netcdf.putAtt(ncid, Clim_dR_ID,'units', 'W m^-2');
   netcdf.putAtt(ncid, Clim_dR_lw_ID,'long_name','LW Radiation Climatology, All-Sky Conditions');
   netcdf.putAtt(ncid, Clim_dR_lw_ID,'units', 'W m^-2');
   netcdf.putAtt(ncid, Clim_dR_sw_ID,'long_name','SW Radiation Climatology, All-Sky Conditions');
   netcdf.putAtt(ncid, Clim_dR_sw_ID,'units', 'W m^-2');
   
%    netcdf.putAtt(ncid, dR_cld_ID,'long_name','Total Cloud Radiative Effects Anomaly, Monthly Means');
%    netcdf.putAtt(ncid, dR_cld_ID,'units', 'W m^-2');
%    netcdf.putAtt(ncid, dR_cld_lw_ID,'long_name','LW Cloud Radiative Effects Anomaly, Monthly Means');
%    netcdf.putAtt(ncid, dR_cld_lw_ID,'units', 'W m^-2');
%    netcdf.putAtt(ncid, dR_cld_sw_ID,'long_name','SW Cloud Radiative Effects Anomaly, Monthly Means');
%    netcdf.putAtt(ncid, dR_cld_sw_ID,'units', 'W m^-2');
%    netcdf.putAtt(ncid, Clim_dR_cld_ID,'long_name','Total Cloud Radiative Effects Climatology');
%    netcdf.putAtt(ncid, Clim_dR_cld_ID,'units', 'W m^-2');
%    netcdf.putAtt(ncid, Clim_dR_cld_lw_ID,'long_name','LW Cloud Radiative Effects Climatology');
%    netcdf.putAtt(ncid, Clim_dR_cld_lw_ID,'units', 'W m^-2');
%    netcdf.putAtt(ncid, Clim_dR_cld_sw_ID,'long_name','SW Cloud Radiative Effects Climatology');
%    netcdf.putAtt(ncid, Clim_dR_cld_sw_ID,'units', 'W m^-2');
   
   %We are done defining the NetCdf
   netcdf.endDef(ncid);
   
   %Then store the dimension variables in
   netcdf.putVar(ncid,longitude_ID,lonf);
   netcdf.putVar(ncid,latitude_ID,latf);
   netcdf.putVar(ncid,time_ID,time);
   netcdf.putVar(ncid,month_ID,month);
   %Then store my main variable
   netcdf.putVar(ncid,dR_ID,dR);
   netcdf.putVar(ncid,dR_lw_ID,dR_lw);
   netcdf.putVar(ncid,dR_sw_ID,dR_sw);
   netcdf.putVar(ncid,Clim_dR_ID,Clim_dR);
   netcdf.putVar(ncid,Clim_dR_lw_ID,Clim_dR_lw);
   netcdf.putVar(ncid,Clim_dR_sw_ID,Clim_dR_sw);
   
%    netcdf.putVar(ncid,dR_cld_ID,dR_cld);
%    netcdf.putVar(ncid,dR_cld_lw_ID,dR_cld_lw);
%    netcdf.putVar(ncid,dR_cld_sw_ID,dR_cld_sw);
%    netcdf.putVar(ncid,Clim_dR_cld_ID,Clim_dR_cld);
%    netcdf.putVar(ncid,Clim_dR_cld_lw_ID,Clim_dR_cld_lw);
%    netcdf.putVar(ncid,Clim_dR_cld_sw_ID,Clim_dR_cld_sw);
   %We're done, close the netcdf
   netcdf.close(ncid)

end

t=toc;disp(t)
