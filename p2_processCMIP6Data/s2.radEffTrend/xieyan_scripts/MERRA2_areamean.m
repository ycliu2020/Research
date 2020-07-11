% clc
clear
%This program is to calculate ERAi surface temperature area means
%% Read in data
filepath1 = '/data1/xieyan/MERRA2/Ts/'; 
ts = zeros(576,361,456);                                                  % 38 years: from 1980/03 - 2018/02
time = zeros(456,1); ntime = length(time);
for ii = 0:37
   for jj = 1:12
      if (ii*12+jj)<143
          vv = 100;
      elseif (ii*12+jj)<251
          vv = 200;
      elseif (ii*12+jj)<371
          vv = 300;
      else
          vv = 400;
      end
      if jj<11
          yy = ii; mm = jj+2;
      else
          yy = ii+1; mm = jj - 10;
      end
      filename1 = strcat(filepath1,'MERRA2_',num2str(vv),'.instM_2d_asm_Nx.',num2str((1980+yy),'%04i'),num2str(mm,'%02i'),'.nc4.nc');
      ts(:,:,ii*12+jj) = ncread(filename1,'TS');
      time(ii*12+jj,1) = datenum((1980+yy),mm,1,0,0,0);
   end
end    
%% Interpolate MERRA2 data
lon = ncread(strcat(filepath1,'MERRA2_200.instM_2d_asm_Nx.200001.nc4.nc'),'lon'); nlon = length(lon);
lat = ncread(strcat(filepath1,'MERRA2_200.instM_2d_asm_Nx.200001.nc4.nc'),'lat'); nlat = length(lat);
time2 = 1:ntime;
ind_lon = lon<0; lon(ind_lon) = lon(ind_lon)+360;
t = lon(289:576); lon(289:576) = lon(1:288);lon(1:288)=t; lon(1) = 0;
t = ts(289:576,:,:); ts(289:576,:,:) = ts(1:288,:,:); ts(1:288,:,:) = t;
lonf = 0.0:1:359.0; nlonf = length(lonf); lonf = double(lonf);
latf = 90.0:-1:-90.0; nlatf = length(latf); latf = double(latf);

[Xlon, Ylat, Ttime] = meshgrid(lat,lon,time2);
[Xlonf, Ylatf, Ttimef] = meshgrid(latf,lonf,time2);
ts1 = interp3(Xlon,Ylat,Ttime,ts,Xlonf,Ylatf,Ttimef);

%% Deseasonalize
[tsAnom,tsClim] = monthlyAnomaly3D(360,181,time,ts1);
%% Land
filename2 = strcat('/data1/xieyan/ERAi/Ts/t1979.nc'); 
sst = ncread(filename2,'sst');
land = isnan(squeeze(sst(:,:,6)));                                                        % 360*181 logical: 1 means land/ 0 means ocean
land = int8(land);
for jj = 1:nlatf
   for ii = 1:nlonf
      if land(ii,jj) == 0
          tsAnom(ii,jj,:) = NaN;
          tsClim(ii,jj,:) = NaN;
      end
   end
end
ts0 = tsAnom;                                       % 360*181*456 surface temperature anomaly, Global land
%% Time series of area mean
% Define the following regions for area mean calculation
% Global land/ Northern Hemisphere(NH) land/ 20-50N land/ China mainland
% 0 - Global land
ts0mean = zeros(ntime,1);
for tt = 1:ntime
   temp1 = squeeze(ts0(:,:,tt));
   ts0mean(tt,1) = area_mean(temp1',ndgrid(latf,lonf)); 
end

% 1 - Northern Hemisphere(NH) land
indNH = latf'>0;                                                          % Transfer to 1*nlat or else exceed matrix dimensions
lat1 = latf(latf>0); lon1 = lonf;
ts11 = ts0(:,indNH,:); 
ts1mean = zeros(ntime,1); 
for tt = 1:ntime
   temp1 = squeeze(ts11(:,:,tt));
   ts1mean(tt,1) = area_mean(temp1',ndgrid(lat1,lon1)); 
end

% 2 - 20-50N land
ind2 = ((latf>=20)&(latf<=50))';                                                          % Transfer to 1*nlat or else exceed matrix dimensions
lat2 = latf(latf>=20 & latf<=50); lon2 = lonf;
ts2 = ts0(:,ind2,:); 
ts2mean = zeros(ntime,1); 
for tt = 1:ntime
   temp1 = squeeze(ts2(:,:,tt));
   ts2mean(tt,1) = area_mean(temp1',ndgrid(lat2,lon2)); 
end

% 3 - China mainland
% longitude70-140 latitude55-15
indlat3 = ((latf>=15)&(latf<=55))'; indlon3 = ((lonf>=70)&(lonf<=140))';
lat3 = latf(latf>=15 & latf<=55); nlat3 =length(lat3);
lon3 = lonf(lonf>=70 & lonf<=140); nlon3 = length(lon3);
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
saveto = strcat('/data1/xieyan/MERRA2_Tsmean.nc');
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