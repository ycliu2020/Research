% clc
clear
% This program is to calculate dRx0/dRx/dRxs0/dRxs and their area-means
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

% filepath1 = '/data1/xieyan/';
filepath1 = '/data1/xieyan/EBAF_clear/';
% filepath1 = 'E:\Undergraduate\大四下\本科毕业论文\data\123';
filename1 = strcat(filepath1,'datasets_radAnom2018.nc');
filename2 = strcat(filepath1,'datasets_radTrend2018.nc');
% filename1 = strcat(filepath1,'datasets_radAnom2018supp.nc');
% filename2 = strcat(filepath1,'datasets_radTrend2018supp.nc');

lon1 = ncread(filename1,'longitude'); 
lat1 = ncread(filename1,'latitude');
time = ncread(filename1,'time'); ntime = length(time);
dRx0 = ncread(filename1,'dRx0'); dR1mean = zeros(5,2,3,2,3,ntime);
dRx = ncread(filename2,'dRx'); dR2mean = zeros(5,2,3,2,3,ntime);
% dRx0 = ncread(filename1,'dRxs0'); dR1mean = zeros(5,2,3,2,3,ntime);
% dRx = ncread(filename2,'dRxs'); dR2mean = zeros(5,2,3,2,3,ntime);
%% Define China region and regrid
latch = lat(lat>0 & lat<55.5); lonch = lon(lon>70 & lon<141);
nlatch =length(latch);nlonch = length(lonch);
indlatch = ((lat>0)&(lat<55.5))'; indlonch = ((lon>70)&(lon<141))';
isinqtp = mask_qtp(indlonch,indlatch);
dR1 = zeros(5,2,3,2,nlonch,nlatch,ntime); dR2 = zeros(5,2,3,2,nlonch,nlatch,ntime);
time2 = 1:ntime;
[Xlon1,Ylat1,Ttime1] = meshgrid(lat1,lon1,time2);
[Xlonch,Ylatch,Ttimech] = meshgrid(latch,lonch,time2); 
for xx = 1:5
   for sky = 1:2
      for spec = 1:3 
         for budg = 1:2
             temp1 = squeeze(dRx0(xx,sky,spec,budg,:,:,:));
             temp1 = interp3(Xlon1,Ylat1,Ttime1,temp1,Xlonch,Ylatch,Ttimech);
             dR1(xx,sky,spec,budg,:,:,:) = temp1;
             temp2 = squeeze(dRx(xx,sky,spec,budg,:,:,:));
             temp2 = interp3(Xlon1,Ylat1,Ttime1,temp2,Xlonch,Ylatch,Ttimech);
             dR2(xx,sky,spec,budg,:,:,:) = temp2;
%              temp3 = squeeze(dRxs0(xx,sky,spec,budg,:,:,:));
%              temp3 = interp3(Xlon1,Ylat1,Ttime1,temp3,Xlonch,Ylatch,Ttimech);
%              dR3(xx,sky,spec,budg,:,:,:) = temp3;
%              temp4 = squeeze(dRxs(xx,sky,spec,budg,:,:,:));
%              temp4 = interp3(Xlon1,Ylat1,Ttime1,temp4,Xlonch,Ylatch,Ttimech);
%              dR4(xx,sky,spec,budg,:,:,:) = temp4;
         end
      end
   end
end

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
   for sky = 1:2
      for spec = 1:3 
         for budg = 1:2
             for tt = 1:ntime
                 temp = squeeze(dR1(xx,sky,spec,budg,:,:,tt));
                 temp1 = temp; temp2 = temp; temp3 = temp;
                 temp1(~isinch) = NaN; 
                 temp2(~isinqtp) = NaN;
                 temp3(~isinch_nonqtp) = NaN;
                 dR1mean(xx,sky,spec,budg,1,tt) = area_mean(temp1',ndgrid(latch,lonch));
                 dR1mean(xx,sky,spec,budg,2,tt) = area_mean(temp2',ndgrid(latch,lonch));
                 dR1mean(xx,sky,spec,budg,3,tt) = area_mean(temp3',ndgrid(latch,lonch));
                 
                 temp = squeeze(dR2(xx,sky,spec,budg,:,:,tt));
                 temp1 = temp; temp2 = temp; temp3 = temp;
                 temp1(~isinch) = NaN; 
                 temp2(~isinqtp) = NaN;
                 temp3(~isinch_nonqtp) = NaN;
                 dR2mean(xx,sky,spec,budg,1,tt) = area_mean(temp1',ndgrid(latch,lonch));
                 dR2mean(xx,sky,spec,budg,2,tt) = area_mean(temp2',ndgrid(latch,lonch));
                 dR2mean(xx,sky,spec,budg,3,tt) = area_mean(temp3',ndgrid(latch,lonch));
             end
         end
      end
   end
end
% Verify we got the right region
temp = squeeze(dR1(1,1,1,1,:,:,1));
temp1 = temp; temp2 = temp; temp3 = temp;
temp1(~isinch) = NaN; 
temp2(~isinqtp) = NaN;
temp3(~isinch_nonqtp) = NaN;

trend_dR1mean = NaN(5,2,3,2,3,12,2);
p_dR1mean = NaN(5,2,3,2,3,12);
dR1mean_detrend = NaN(size(dR1mean));
for xx = 1:5
   for sky = 1:2
      for spec = 1:3 
         for budg = 1:2
             for region = 1:3
                [dR1mean_detrend(xx,sky,spec,budg,region,:),trend_dR1mean(xx,sky,spec,budg,region,:,:),~,p_dR1mean(xx,sky,spec,budg,region,:)] = detrend_yan(dR1mean(xx,sky,spec,budg,region,:),time);
             end
         end
      end
   end
end


%% Save to NC file
five= 1:5; three = 1:3; two = 1:2; month = 1:12;
% saveto = strcat('/data1/xieyan/Rad_Chinamean.nc'); 
% saveto = strcat('/data1/xieyan/Radsupp_Chinamean.nc');
% saveto = strcat('/data1/xieyan/EBAF_clear/Radsupp_Chinamean.nc');
saveto = strcat('/data1/xieyan/EBAF_clear/Rad_Chinamean.nc');
ncid = netcdf.create(saveto,'NC_WRITE');
%Define the dimensions
dimidx = netcdf.defDim(ncid,'x',5);
dimidsky = netcdf.defDim(ncid,'sky',2);
dimidspec = netcdf.defDim(ncid,'spectrum',3);
dimidbudget = netcdf.defDim(ncid,'budget',2);
dimidregion = netcdf.defDim(ncid,'region',3);
dimidcoef = netcdf.defDim(ncid,'coef',2);
dimidlon = netcdf.defDim(ncid,'longitude',nlonch);
dimidlat = netcdf.defDim(ncid,'latitude',nlatch);
dimidtime = netcdf.defDim(ncid,'time',ntime);
dimidmonth = netcdf.defDim(ncid,'month',12);
%Define IDs for the dimension variables 
x_ID=netcdf.defVar(ncid,'x','int',dimidx);
sky_ID =netcdf.defVar(ncid,'sky','int',dimidsky);
spec_ID=netcdf.defVar(ncid,'spectrum','int',dimidspec);
budget_ID =netcdf.defVar(ncid,'budget','int',dimidbudget);
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
dR1mean_ID = netcdf.defVar(ncid,'dR1mean','double',[dimidx dimidsky dimidspec dimidbudget dimidregion dimidtime]);
dR2mean_ID = netcdf.defVar(ncid,'dR2mean','double',[dimidx dimidsky dimidspec dimidbudget dimidregion dimidtime]);
dR1mean_detrend_ID = netcdf.defVar(ncid,'dR1mean_detrend','double',[dimidx dimidsky dimidspec dimidbudget dimidregion dimidtime]);
dR1meantrend_ID = netcdf.defVar(ncid,'dR1meantrend','double',[dimidx dimidsky dimidspec dimidbudget dimidregion dimidmonth dimidcoef]);
pvalue_ID = netcdf.defVar(ncid,'pvalue','double',[dimidx dimidsky dimidspec dimidbudget dimidregion dimidmonth]);
%Define units for the dimension variables
netcdf.putAtt(ncid,x_ID,'standard_name','x');
netcdf.putAtt(ncid,x_ID,'units','Total,Temperature,Water vapor,Albedo,Cloud');
% netcdf.putAtt(ncid,x_ID,'units','non-cloud, diagnosed all, residual, ta, ts');
netcdf.putAtt(ncid,sky_ID,'standard_name','sky');
netcdf.putAtt(ncid,sky_ID,'units','all-sky; clear-sky');
netcdf.putAtt(ncid,spec_ID,'standard_name','spectrum');
netcdf.putAtt(ncid,spec_ID,'units','NET;LW;SW');
netcdf.putAtt(ncid,budget_ID,'standard_name','budget');
netcdf.putAtt(ncid,budget_ID,'units','TOA; SFC');
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

netcdf.putAtt(ncid,dR1mean_ID,'long_name','Radiative Effect, Deseasonalized, Area-mean');
netcdf.putAtt(ncid,dR1mean_ID,'units', 'W * m^2');
netcdf.putAtt(ncid,dR2mean_ID,'long_name','Radiative Effect, Deseasonalized and Monthly Detrended, Area-mean');
netcdf.putAtt(ncid,dR2mean_ID,'units', 'W * m^2');
netcdf.putAtt(ncid,temp1_ID,'long_name','Radiative Effect, Deseasonalized, China land(as a test)');
netcdf.putAtt(ncid,temp1_ID,'units', 'W * m^2');
netcdf.putAtt(ncid,temp2_ID,'long_name','Radiative Effect, Deseasonalized, Qinghai-Tibet Plateau(as a test)');
netcdf.putAtt(ncid,temp2_ID,'units', 'W * m^2');
netcdf.putAtt(ncid,temp3_ID,'long_name','Radiative Effect, Deseasonalized, China exclude Qinghai-Tibet Plateau(as a test)');
netcdf.putAtt(ncid,temp3_ID,'units', 'W * m^2');
netcdf.putAtt(ncid,dR1mean_detrend_ID,'long_name','Radiative Effect, Deseasonalized, Area-mean, then detrended');
netcdf.putAtt(ncid,dR1mean_detrend_ID,'units', 'W * m^2');
netcdf.putAtt(ncid,dR1meantrend_ID,'long_name','Trend of Deseasonalized Area-mean Radiative Effect');
netcdf.putAtt(ncid,dR1meantrend_ID,'units', 'W * m^2');
netcdf.putAtt(ncid,pvalue_ID,'long_name','P Value of trending');
netcdf.putAtt(ncid,pvalue_ID,'units', 'W * m^2');
%We are done defining the NetCdf
netcdf.endDef(ncid);
%Then store the dimension variables in
netcdf.putVar(ncid,x_ID,five);
netcdf.putVar(ncid,sky_ID,two);
netcdf.putVar(ncid,spec_ID,three);
netcdf.putVar(ncid,budget_ID,two);
netcdf.putVar(ncid,region_ID,three);
netcdf.putVar(ncid,longitude_ID,lonch);
netcdf.putVar(ncid,latitude_ID,latch);
netcdf.putVar(ncid,time_ID,time);
netcdf.putVar(ncid,month_ID,month);
netcdf.putVar(ncid,coef_ID,two);
%Then store my main variable
netcdf.putVar(ncid,dR1mean_ID,dR1mean);
netcdf.putVar(ncid,dR2mean_ID,dR2mean);
netcdf.putVar(ncid,temp1_ID,temp1);
netcdf.putVar(ncid,temp2_ID,temp2);
netcdf.putVar(ncid,temp3_ID,temp3);
netcdf.putVar(ncid,dR1mean_detrend_ID,dR1mean_detrend);
netcdf.putVar(ncid,dR1meantrend_ID,trend_dR1mean);
netcdf.putVar(ncid,pvalue_ID,p_dR1mean);
%We're done, close the netcdf
netcdf.close(ncid)