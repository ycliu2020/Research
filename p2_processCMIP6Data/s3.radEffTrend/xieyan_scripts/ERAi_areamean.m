% clc
clear
%This program is to calculate ERAi surface temperature area means
%% Read in data
filename1 = '/data1/xieyan/ERAi/Ts/ERAi_tsAnom.nc'; 
ts0 = ncread(filename1,'Anom_ts');                                        % 360*181*468 surface temperature anomaly, Global land
lon = ncread(filename1,'longitude'); nlon = length(lon);
lat = ncread(filename1,'latitude'); nlat = length(lat);
time = ncread(filename1,'time'); ntime = length(time);                    % 1979/03 - 2018/02
%% Time series of area mean
% Define the following regions for area mean calculation
% Global land/ Northern Hemisphere(NH) land/ 20-50N land/ China mainland
% 0 - Global land
ts0mean = zeros(ntime,1);
for tt = 1:ntime
   temp1 = squeeze(ts0(:,:,tt));
   ts0mean(tt,1) = area_mean(temp1',ndgrid(lat,lon)); 
end

% 1 - Northern Hemisphere(NH) land
indNH = lat'>0;                                                          % Transfer to 1*nlat or else exceed matrix dimensions
lat1 = lat(lat>0); lon1 = lon;
ts1 = ts0(:,indNH,:); 
ts1mean = zeros(ntime,1); 
for tt = 1:ntime
   temp1 = squeeze(ts1(:,:,tt));
   ts1mean(tt,1) = area_mean(temp1',ndgrid(lat1,lon1)); 
end

% 2 - 20-50N land
ind2 = ((lat>=20)&(lat<=50))';                                                          % Transfer to 1*nlat or else exceed matrix dimensions
lat2 = lat(lat>=20 & lat<=50); lon2 = lon;
ts2 = ts0(:,ind2,:); 
ts2mean = zeros(ntime,1); 
for tt = 1:ntime
   temp1 = squeeze(ts2(:,:,tt));
   ts2mean(tt,1) = area_mean(temp1',ndgrid(lat2,lon2)); 
end

% 3 - China mainland
% longitude70-140 latitude55-15
indlat3 = ((lat>=15)&(lat<=55))'; indlon3 = ((lon>=70)&(lon<=140))';
lat3 = lat(lat>=15 & lat<=55); nlat3 =length(lat3);
lon3 = lon(lon>=70 & lon<=140); nlon3 = length(lon3);
ts3 = ts0(indlon3,indlat3,:); 
map = '/data1/xieyan/maps/worldcountry.shp';
readmap = shaperead(map);
[Xlon3,Ylat3]=meshgrid(lon3,lat3);                                  % [Xlon Ylat] = meshgrid(lon,lat);
isin = inpolygon(Xlon3,Ylat3,readmap(49,1).X,readmap(49,1).Y);
for tt = 1:ntime
    temp1 = squeeze(ts3(:,:,tt));
    temp1 = temp1';
    temp1(~isin)=NaN;
    ts3(:,:,tt) = temp1';
end
ts3mean = zeros(ntime,1); 
for tt = 1:ntime
   temp1 = squeeze(ts3(:,:,tt));
   ts3mean(tt,1) = area_mean(temp1',ndgrid(lat3,lon3)); 
end
%% Save to NC file
saveto = strcat('/data1/xieyan/ERAi_Tsmean.nc');
ncid = netcdf.create(saveto,'NC_WRITE');
%Define the dimensions
dimidlon = netcdf.defDim(ncid,'longitude',nlon3);
dimidlat = netcdf.defDim(ncid,'latitude',nlat3);
dimidtime = netcdf.defDim(ncid,'time',ntime);
%Define IDs for the dimension variables 
longitude_ID=netcdf.defVar(ncid,'longitude','double',dimidlon);
latitude_ID=netcdf.defVar(ncid,'latitude','double',dimidlat);
time_ID=netcdf.defVar(ncid,'time','double',dimidtime);
%Define the main variable ()
tslandm_ID = netcdf.defVar(ncid,'tslandm','double',dimidtime);
tsnhm_ID = netcdf.defVar(ncid,'tsnhm','double',dimidtime);
tsmidm_ID = netcdf.defVar(ncid,'tsmidm','double',dimidtime);
tschm_ID = netcdf.defVar(ncid,'tschm','double',dimidtime);
tsch_ID = netcdf.defVar(ncid,'tsch','double',[dimidlon dimidlat dimidtime]);

%Define units for the dimension variables
netcdf.putAtt(ncid,longitude_ID,'standard_name','longitude');
netcdf.putAtt(ncid,latitude_ID,'standard_name','latitude');
netcdf.putAtt(ncid,time_ID,'long_name','time');
netcdf.putAtt(ncid,time_ID,'units','days since 0000-00-00 00:00:0.0');
netcdf.putAtt(ncid,time_ID,'calendar','gregorian');

netcdf.putAtt(ncid,tslandm_ID,'long_name','Surface Temperature Anomaly, Global land area-mean');
netcdf.putAtt(ncid,tslandm_ID,'units', 'K');
netcdf.putAtt(ncid,tsnhm_ID,'long_name','Surface Temperature Anomaly, Northern Hemisphere land area-mean');
netcdf.putAtt(ncid,tsnhm_ID,'units', 'K');
netcdf.putAtt(ncid,tsmidm_ID,'long_name','Surface Temperature Anomaly, 20 - 50N land area-mean');
netcdf.putAtt(ncid,tsmidm_ID,'units', 'K');
netcdf.putAtt(ncid,tschm_ID,'long_name','Surface Temperature Anomaly, China land area-mean');
netcdf.putAtt(ncid,tschm_ID,'units', 'K');
netcdf.putAtt(ncid,tsch_ID,'long_name','Surface Temperature Anomaly, China land');
netcdf.putAtt(ncid,tsch_ID,'units', 'K');
%We are done defining the NetCdf
netcdf.endDef(ncid);

%Then store the dimension variables in
netcdf.putVar(ncid,longitude_ID,lon3);
netcdf.putVar(ncid,latitude_ID,lat3);
netcdf.putVar(ncid,time_ID,time);

%Then store my main variable
netcdf.putVar(ncid,tslandm_ID,ts0mean);
netcdf.putVar(ncid,tsnhm_ID,ts1mean);
netcdf.putVar(ncid,tsmidm_ID,ts2mean);
netcdf.putVar(ncid,tschm_ID,ts3mean);
netcdf.putVar(ncid,tsch_ID,ts3);
%We're done, close the netcdf
netcdf.close(ncid)