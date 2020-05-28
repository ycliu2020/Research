%% Combine the ERAi climatological data from 2000-03 to 2018-02
clc;clear; tic;
%modify path first
filepath1 = '/data1/xieyan/ERAi/Surface_ori/';
filepath2 = '/data1/xieyan/ERAi/TQ_ori/';
a0 = zeros(360,181,228);                                                  % Forecast albedo
ts0 = zeros(360,181,228);                                                 % Skin temperature
q0 = zeros(360,181,37,228);                                               % Specific humidity
t0 = zeros(360,181,37,228);                                               % Temperature
time0 = zeros(228,1);
for ii = 0:18
    filename1 = strcat(filepath1,'sfc',num2str((ii+2000),'%04i'),'.nc'); 
    filename2 = strcat(filepath2,'tq',num2str((ii+2000),'%04i'),'.nc');
    temp1 = ncread(filename1,'fal');
    temp2 = ncread(filename1,'skt');
    temp3 = ncread(filename2,'q');
    temp4 = ncread(filename2,'t');
    temp5 = ncread(filename2,'time');
    a0(:,:,ii*12+1:(ii+1)*12) = temp1;
    ts0(:,:,ii*12+1:(ii+1)*12) = temp2;
    q0(:,:,:,ii*12+1:(ii+1)*12) = temp3;
    t0(:,:,:,ii*12+1:(ii+1)*12) = temp4;
    time0(ii*12+1:(ii+1)*12,1) = temp5;
end
a = a0(:,:,3:218);
ts = ts0(:,:,3:218);
q = q0(:,:,:,3:218);
t = t0(:,:,:,3:218);
time = time0(3:218,1); 
time = time./24;
time = datenum(1900,01,01)+time; ntime = length(time);

%% Regrid to 2.5*2.5
lon = 0:359; lat = 90:-1:-90;
plev = ncread('/data1/xieyan/ERAi/TQ_ori/tq2017.nc','level');
plev = double(plev);
% time2 = 1:348;                                                   %29*12 month
time2 = 1:216;                                                       %18*12 month
[Xlon, Ylat, Zplev, Ttime] = ndgrid(lon,lat,plev,time2);      
lonf = 0:2.5:357.5; nlonf = length(lonf);
latf = 90:-2.5:-90; nlatf = length(latf);
plevf1 = ncread('/data/pub/kernels_YiH/toa/dp.nc','player');  % pay attention to plevf's range must smaller than plev's
% kernel��������ϵ���ô�McGill University Prof.Yi Huang Email: yi.huang@mcgill.ca
plevf = plevf1(end:-1:1);
plevfnum=length(plevf);
[Xlonf, Ylatf, Zplevf, Ttimef] = ndgrid(lonf,latf,plevf,time2);      
q = interpn(Xlon,Ylat,Zplev,Ttime,q,Xlonf,Ylatf,Zplevf,Ttimef);
t = interpn(Xlon,Ylat,Zplev,Ttime,t,Xlonf,Ylatf,Zplevf,Ttimef);

[Xlon, Ylat, Ttime] = meshgrid(lat,lon,time2);
[Xlonf, Ylatf, Ttimef] = meshgrid(latf,lonf,time2);
ts = interp3(Xlon,Ylat,Ttime,ts,Xlonf,Ylatf,Ttimef);
a = interp3(Xlon,Ylat,Ttime,a,Xlonf,Ylatf,Ttimef);

%% Deseasonalize
[Anom_alb,Clim_alb] = monthlyAnomaly3D(144,73,time,a);
[Anom_ts,Clim_ts] = monthlyAnomaly3D(144,73,time,ts);
[Anom_q,Clim_q] = monthlyAnomaly4D(144,73,plevfnum,time,q);
[Anom_t,Clim_t] = monthlyAnomaly4D(144,73,plevfnum,time,t);

month = 1:12;

%% Save to nc
directory_name = '/home/lyc/repeat_Kernel_experiment/testdata/';
% saveto = strcat(directory_name,'ERAi_variAnom1988.nc');
saveto = strcat(directory_name,'ERAi_variAnom2018.nc');
ncid = netcdf.create(saveto,'NC_WRITE');
%Define the dimensions
dimidlon = netcdf.defDim(ncid,'longitude',nlonf);
dimidlat = netcdf.defDim(ncid,'latitude',nlatf);
dimidlevel = netcdf.defDim(ncid,'level',length(plevf));
dimidtime = netcdf.defDim(ncid,'time',ntime);
dimidmonth = netcdf.defDim(ncid,'month',12);

%Define IDs for the dimension variables 
longitude_ID=netcdf.defVar(ncid,'longitude','double',dimidlon);
latitude_ID=netcdf.defVar(ncid,'latitude','double',dimidlat);
level_ID = netcdf.defVar(ncid,'level','double',dimidlevel);
time_ID=netcdf.defVar(ncid,'time','double',dimidtime);
month_ID=netcdf.defVar(ncid,'month','double',dimidmonth);

%Define the main variable ()
Anom_alb_ID = netcdf.defVar(ncid,'Anom_alb','double',[dimidlon dimidlat dimidtime]);
Anom_ts_ID = netcdf.defVar(ncid,'Anom_ts','double',[dimidlon dimidlat dimidtime]);
Anom_q_ID = netcdf.defVar(ncid,'Anom_q','double',[dimidlon dimidlat dimidlevel dimidtime]);
Anom_t_ID = netcdf.defVar(ncid,'Anom_t','double',[dimidlon dimidlat dimidlevel dimidtime]);
Clim_alb_ID = netcdf.defVar(ncid,'Clim_alb','double',[dimidlon dimidlat dimidmonth]);
Clim_ts_ID = netcdf.defVar(ncid,'Clim_ts','double',[dimidlon dimidlat dimidmonth]);
Clim_q_ID = netcdf.defVar(ncid,'Clim_q','double',[dimidlon dimidlat dimidlevel dimidmonth]);
Clim_t_ID = netcdf.defVar(ncid,'Clim_t','double',[dimidlon dimidlat dimidlevel dimidmonth]);
q_ID = netcdf.defVar(ncid,'q','double',[dimidlon dimidlat dimidlevel dimidtime]);

%Define units for the dimension variables
netcdf.putAtt(ncid,longitude_ID,'standard_name','longitude');
netcdf.putAtt(ncid,latitude_ID,'standard_name','latitude');
netcdf.putAtt(ncid,level_ID,'standard_name','Pressure levels');
netcdf.putAtt(ncid,time_ID,'long_name','time');
netcdf.putAtt(ncid,time_ID,'units','days since 0000-00-00 00:00:0.0');
netcdf.putAtt(ncid,time_ID,'calendar','gregorian');
netcdf.putAtt(ncid,month_ID,'long_name','month');
netcdf.putAtt(ncid,month_ID,'units','Month (JFMAMJJASOND)');

netcdf.putAtt(ncid,Anom_alb_ID,'long_name','Forecast Albedo Anomaly');
netcdf.putAtt(ncid,Anom_alb_ID,'units', '1');
netcdf.putAtt(ncid,Anom_ts_ID,'long_name','Skin Temperature Anomaly');
netcdf.putAtt(ncid,Anom_ts_ID,'units', 'K');
netcdf.putAtt(ncid,Anom_q_ID,'long_name','Specific Humidity Anomaly');
netcdf.putAtt(ncid,Anom_q_ID,'units', 'kg kg**-1');
netcdf.putAtt(ncid,Anom_t_ID,'long_name','Air Temperature Anomaly');
netcdf.putAtt(ncid,Anom_t_ID,'units', 'K');
netcdf.putAtt(ncid,q_ID,'long_name','Specific Humidity');
netcdf.putAtt(ncid,q_ID,'units', 'kg kg**-1');

netcdf.putAtt(ncid,Clim_alb_ID,'long_name','Forecast Albedo Climatology');
netcdf.putAtt(ncid,Clim_alb_ID,'units', '1');
netcdf.putAtt(ncid,Clim_ts_ID,'long_name','Skin Temperature Climatology');
netcdf.putAtt(ncid,Clim_ts_ID,'units', 'K');
netcdf.putAtt(ncid,Clim_q_ID,'long_name','Specific Humidity Climatology');
netcdf.putAtt(ncid,Clim_q_ID,'units', 'kg kg**-1');
netcdf.putAtt(ncid,Clim_t_ID,'long_name','Air Temperature Climatology');
netcdf.putAtt(ncid,Clim_t_ID,'units', 'K');
%We are done defining the NetCdf
netcdf.endDef(ncid);

%Then store the dimension variables in
netcdf.putVar(ncid,longitude_ID,lonf);
netcdf.putVar(ncid,latitude_ID,latf);
netcdf.putVar(ncid,level_ID,plevf);
netcdf.putVar(ncid,time_ID,time);
netcdf.putVar(ncid,month_ID,month);

%Then store my main variable
netcdf.putVar(ncid,Anom_alb_ID,Anom_alb);
netcdf.putVar(ncid,Anom_ts_ID,Anom_ts);
netcdf.putVar(ncid,Anom_q_ID,Anom_q);
netcdf.putVar(ncid,Anom_t_ID,Anom_t);
netcdf.putVar(ncid,q_ID,q);
netcdf.putVar(ncid,Clim_alb_ID,Clim_alb);
netcdf.putVar(ncid,Clim_ts_ID,Clim_ts);
netcdf.putVar(ncid,Clim_q_ID,Clim_q);
netcdf.putVar(ncid,Clim_t_ID,Clim_t);
%We're done, close the netcdf
netcdf.close(ncid)