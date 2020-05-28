% clc
clear
%This program is to calculate the HadCRUT4 surface temperature area means
%% Read in data
filename1 = '/data1/xieyan/HadCRUT4/HadCRUT.4.6.0.0.median.nc'; 
filename2 = '/data1/xieyan/ERAi/Ts/t1979.nc'; 
% filename1 = 'E:\Undergraduate\大四下\本科毕业论文\data\HadCRUT.4.6.0.0.median.nc'; 
% filename2 = 'E:\Undergraduate\大四下\本科毕业论文\data\ERAi\Ts\t1979.nc';
%Had
ts0 = ncread(filename1,'temperature_anomaly');
ts0 = ts0(:,:,1551:2018);                                                 % 72*36*468 surface temperature anomaly, Global land
lon = ncread(filename1,'longitude'); nlon = length(lon);                  % -177.5:5:177.5
lat = ncread(filename1,'latitude'); nlat = length(lat);                   % -87.5:5:87.5
lat = double(lat); lon = double(lon);
time = ncread(filename1,'time');  time = time(1551:2018);                 % 1979/03 - 2018/02
ntime = length(time); 
time = time + datenum(1850,1,1);
%ERA
sst = ncread(filename2,'sst');
latf = -90:1:90; nlatf = length(latf); latf = double(latf);
lonf = -179:1:180; nlonf = length(lonf); lonf = double(lonf);
sst = squeeze(sst(:,:,1));                                                % 360*181
sst1 = NaN(size(sst));
sst1(1:180,:) = sst(181:360,:); sst1(181:360,:) = sst(1:180,:);
sst1 = flip(sst1,2);
land = isnan(sst1);
land = double(land);

% Regraded 
[Xlon,Ylat] = meshgrid(lat,lon);
[Xlonf,Ylatf] = meshgrid(latf,lonf);
land0 = interp2(Xlonf,Ylatf,land,Xlon,Ylat);
for ii = 1:nlon
   for jj = 1:nlat
      if land0(ii,jj)< 0.1
         ts0(ii,jj,:) = NaN;
      end
   end
end
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
map = '/gs/project/ktg-565-aa/xieyan/map/worldcountry.shp';
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
saveto = strcat('/data1/xieyan/HadCRUT4/HadCRUT4_Tsmean.nc');
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