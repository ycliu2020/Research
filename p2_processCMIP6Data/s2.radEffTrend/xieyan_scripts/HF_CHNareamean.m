% clc
clear
% This program is to calculate ERAi heat flux and surface thermal radiation flux's area-means
% in three regions: Qinghai-Tibet Plateau/ China exclude(QTP)/ China
load /data1/xieyan/mask720.mat
% load E:\Undergraduate\大四下\本科毕业论文\data\mask720
lon0 = lon720; lon0 = double(lon0);
lat0 = lat720; lat0 = double(lat0);
% mask_qtp0 = int8(mask2);
mask_qtp0 = mask2;
% Deal with lon-lat
lon = zeros(size(lon0)); 
lon(1:360,1) = lon0(361:720,1); 
lon(361:720,1) = lon0(1:360,1)+360.0;
lat = flip(lat0,1);
nlon =length(lon);
nlat = length(lat); 
mask_qtp = zeros(size(mask_qtp0));
mask_qtp(1:360,:) = mask_qtp0(361:720,:);
mask_qtp(361:720,:) = mask_qtp0(1:360,:);
mask_qtp = flip(mask_qtp,2);

filepath1 = '/data1/xieyan/ERAi/Heat_Flux/';
% filepath1 = 'E:\Undergraduate\大四下\本科毕业论文\data\';
filename1 = strcat(filepath1,'ERAi_hfAnom.nc');

lon1 = ncread(filename1,'longitude'); 
lat1 = ncread(filename1,'latitude');
time = ncread(filename1,'time'); ntime = length(time);
lahf = ncread(filename1,'lahf'); 
sehf = ncread(filename1,'sehf'); 
slwr = ncread(filename1,'slwr'); 
slwrd = ncread(filename1,'slwrd'); 
slwru = ncread(filename1,'slwru');
%% Define China region and regrid
latch = lat(lat>0 & lat<55.5); lonch = lon(lon>70 & lon<141);
nlatch =length(latch);nlonch = length(lonch);
indlatch = ((lat>0)&(lat<55.5))'; indlonch = ((lon>70)&(lon<141))';
isinqtp = mask_qtp(indlonch,indlatch);
time2 = 1:ntime;
[Xlon1,Ylat1,Ttime1] = meshgrid(lat1,lon1,time2);
[Xlonch,Ylatch,Ttimech] = meshgrid(latch,lonch,time2); 
lahf = interp3(Xlon1,Ylat1,Ttime1,lahf,Xlonch,Ylatch,Ttimech);
sehf = interp3(Xlon1,Ylat1,Ttime1,sehf,Xlonch,Ylatch,Ttimech);
slwr = interp3(Xlon1,Ylat1,Ttime1,slwr,Xlonch,Ylatch,Ttimech);
slwrd = interp3(Xlon1,Ylat1,Ttime1,slwrd,Xlonch,Ylatch,Ttimech);
slwru = interp3(Xlon1,Ylat1,Ttime1,slwru,Xlonch,Ylatch,Ttimech);
% Combine variables together
dFlux = zeros(5,nlonch,nlatch,ntime);
dFlux(1,:,:,:) = lahf;
dFlux(2,:,:,:) = sehf;
dFlux(3,:,:,:) = slwr;
dFlux(4,:,:,:) = slwrd;
dFlux(5,:,:,:) = slwru;
dFlux_mean = zeros(5,3,ntime);                                            % Variable, Region
% China Land Map
map1 = '/data1/xieyan/maps/worldcountry.shp';
map2 = '/data1/xieyan/maps/Chinaregion.shp';
map3 = '/data1/xieyan/maps/bou2_4l.shp';
readmap1 = shaperead(map1);
readmap2 = shaperead(map2);
readmap3 = shaperead(map3);
[Xlonch,Ylatch]=meshgrid(lonch,latch);  
isin1 = inpolygon(Xlonch,Ylatch,readmap1(49,1).X,readmap1(49,1).Y);         % mainland
isin2 = inpolygon(Xlonch,Ylatch,readmap2(6,1).X,readmap2(6,1).Y);           % Taiwan Province
isin3 = inpolygon(Xlonch,Ylatch,readmap2(5,1).X,readmap2(5,1).Y);           % Xizang Province
isin12 = isin1|isin2; 
isinch = isin12|isin3;
isinqtp = isinqtp';
isinch_nonqtp = isinch & (~isinqtp);
isinch = isinch'; 
isinqtp = isinqtp';
isinch_nonqtp = isinch_nonqtp';
for xx = 1:5
    for tt = 1:ntime
        temp = squeeze(dFlux(xx,:,:,tt));
        temp1 = temp; temp2 = temp; temp3 = temp;
        temp1(~isinch) = NaN; 
        temp2(~isinqtp) = NaN;
        temp3(~isinch_nonqtp) = NaN;
        dFlux_mean(xx,1,tt) = area_mean(temp1',ndgrid(latch,lonch));
        dFlux_mean(xx,2,tt) = area_mean(temp2',ndgrid(latch,lonch));
        dFlux_mean(xx,3,tt) = area_mean(temp3',ndgrid(latch,lonch));
   end
end
% Verify we got the right region
temp = squeeze(dFlux(1,:,:,1));
temp1 = temp; temp2 = temp; temp3 = temp;
temp1(~isinch) = NaN; 
temp2(~isinqtp) = NaN;
temp3(~isinch_nonqtp) = NaN;
%% Detrend
trend_dFluxmean = NaN(5,3,12,2);
p_dFluxmean = NaN(5,3,12);
dFluxmean_detrend = NaN(size(dFlux_mean));
for xx = 1:5
    for region = 1:3
       [dFluxmean_detrend(xx,region,:),trend_dFluxmean(xx,region,:,:),~,p_dFluxmean(xx,region,:)] = detrend_yan(dFlux_mean(xx,region,:),time);
    end
end
%% Save to NC file
five= 1:5; three = 1:3; two = 1:2; month = 1:12;
saveto = strcat('/data1/xieyan/ERAiHF_Chinamean.nc'); 
ncid = netcdf.create(saveto,'NC_WRITE');
%Define the dimensions
dimidx = netcdf.defDim(ncid,'x',5);
dimidregion = netcdf.defDim(ncid,'region',3);
dimidcoef = netcdf.defDim(ncid,'coef',2);
dimidlon = netcdf.defDim(ncid,'longitude',nlonch);
dimidlat = netcdf.defDim(ncid,'latitude',nlatch);
dimidtime = netcdf.defDim(ncid,'time',ntime);
dimidmonth = netcdf.defDim(ncid,'month',12);
%Define IDs for the dimension variables 
x_ID=netcdf.defVar(ncid,'x','int',dimidx);
region_ID = netcdf.defVar(ncid,'region','int',dimidregion);
longitude_ID=netcdf.defVar(ncid,'longitude','double',dimidlon);
latitude_ID=netcdf.defVar(ncid,'latitude','double',dimidlat);
time_ID=netcdf.defVar(ncid,'time','double',dimidtime);
month_ID = netcdf.defVar(ncid,'month','double',dimidmonth);
coef_ID = netcdf.defVar(ncid,'coef','double',dimidcoef);
%Define the main variable ()
temp1_ID =  netcdf.defVar(ncid,'temp1','double',[dimidlon dimidlat]);
temp2_ID =  netcdf.defVar(ncid,'temp2','double',[dimidlon dimidlat]);
temp3_ID =  netcdf.defVar(ncid,'temp3','double',[dimidlon dimidlat]);
dFlux_mean_ID = netcdf.defVar(ncid,'dFlux_mean','double',[dimidx dimidregion dimidtime]);
dFluxmean_detrend_ID = netcdf.defVar(ncid,'dFluxmean_detrend','double',[dimidx dimidregion dimidtime]);
trend_dFluxmean_ID = netcdf.defVar(ncid,'trend_dFluxmean','double',[dimidx dimidregion dimidmonth dimidcoef]);
p_dFluxmean_ID = netcdf.defVar(ncid,'p_dFluxmean','double',[dimidx dimidregion dimidmonth]);
%Define units for the dimension variables
netcdf.putAtt(ncid,x_ID,'standard_name','x');
netcdf.putAtt(ncid,x_ID,'units','Latent Heat,Sensible Heat,SFC LW net,SFC LW downward,SFC LW upward');
netcdf.putAtt(ncid,region_ID,'standard_name','region');
netcdf.putAtt(ncid,region_ID,'units','China; Qinghai-Tibet Plateau; China exclude Qinghai-Tibet Plateau');
netcdf.putAtt(ncid,longitude_ID,'standard_name','longitude');
netcdf.putAtt(ncid,latitude_ID,'standard_name','latitude');
netcdf.putAtt(ncid,time_ID,'long_name','time');
netcdf.putAtt(ncid,time_ID,'units','days since 0000-00-00 00:00:0.0');
netcdf.putAtt(ncid,time_ID,'calendar','gregorian');
netcdf.putAtt(ncid,month_ID,'long_name','month');
netcdf.putAtt(ncid,month_ID,'units','Month (JFMAMJJASOND)');
netcdf.putAtt(ncid,coef_ID,'units','a; b');

netcdf.putAtt(ncid,dFlux_mean_ID,'long_name','ERA-interim Surface Flux, Deseasonalized, Area-mean');
netcdf.putAtt(ncid,dFlux_mean_ID,'units', 'W * m^2');
netcdf.putAtt(ncid,temp1_ID,'long_name','Surface Flux, Deseasonalized, China land(as a test)');
netcdf.putAtt(ncid,temp1_ID,'units', 'W * m^2');
netcdf.putAtt(ncid,temp2_ID,'long_name','Surface Flux, Deseasonalized, Qinghai-Tibet Plateau(as a test)');
netcdf.putAtt(ncid,temp2_ID,'units', 'W * m^2');
netcdf.putAtt(ncid,temp3_ID,'long_name','Surface Flux, Deseasonalized, China exclude Qinghai-Tibet Plateau(as a test)');
netcdf.putAtt(ncid,temp3_ID,'units', 'W * m^2');
netcdf.putAtt(ncid,dFluxmean_detrend_ID,'long_name','ERA-interim Surface Flux, Deseasonalized, Area-mean, then detrended');
netcdf.putAtt(ncid,dFluxmean_detrend_ID,'units', 'W * m^2');
netcdf.putAtt(ncid,trend_dFluxmean_ID,'long_name','Trend of Deseasonalized Area-mean Surface Flux');
netcdf.putAtt(ncid,trend_dFluxmean_ID,'units', 'W * m^2');
netcdf.putAtt(ncid,p_dFluxmean_ID,'long_name','P Value of trending');
netcdf.putAtt(ncid,p_dFluxmean_ID,'units', 'W * m^2');
%We are done defining the NetCdf
netcdf.endDef(ncid);
%Then store the dimension variables in
netcdf.putVar(ncid,x_ID,five);
netcdf.putVar(ncid,region_ID,three);
netcdf.putVar(ncid,longitude_ID,lonch);
netcdf.putVar(ncid,latitude_ID,latch);
netcdf.putVar(ncid,time_ID,time);
netcdf.putVar(ncid,month_ID,month);
netcdf.putVar(ncid,coef_ID,two);
%Then store my main variable
netcdf.putVar(ncid,dFlux_mean_ID,dFlux_mean);
netcdf.putVar(ncid,temp1_ID,temp1);
netcdf.putVar(ncid,temp2_ID,temp2);
netcdf.putVar(ncid,temp3_ID,temp3);
netcdf.putVar(ncid,dFluxmean_detrend_ID,dFluxmean_detrend);
netcdf.putVar(ncid,trend_dFluxmean_ID,trend_dFluxmean);
netcdf.putVar(ncid,p_dFluxmean_ID,p_dFluxmean);
%We're done, close the netcdf
netcdf.close(ncid)