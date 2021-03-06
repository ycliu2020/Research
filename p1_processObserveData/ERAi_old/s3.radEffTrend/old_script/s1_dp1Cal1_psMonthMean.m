%% calculate surface pressure from 2008-01 to 2012-12(5 years)
clc;clear; 
%modify path first
inputpath = '/data1/liuyincheng/Observe-rawdata/ERAi/Ps/';
outpath= '/data1/liuyincheng/Observe-process/ERAi/kernelsCal/';
%creat 19yrs to cal
ps_0 = zeros(360,181,228);
time0 = zeros(228,1);
for ii = 0:18
    filename = strcat(inputpath,'ps',num2str((ii+2000),'%04i'),'.nc');
    temp1 = ncread(filename,'sp');        % surface pressure
    temp2 = ncread(filename,'time');       % hours since 1900-01-01 00:00:00.0
        % sum up the cumulative radiation of time step 0-12 and 12-24
        ps_0(:,:,ii*12+1:(ii+1)*12) = temp1;
        time0(ii*12+1:(ii+1)*12,1) = temp2;
end

%以下为掐头去???2000???3月至2018???2月期???
ps = ps_0(:,:,85:144); 
time = time0(85:144,1);
time = datenum(1900,01,01,double(time),00,00);
ntime = length(time);
ps = ps./100; %转换为hp单位

%% Regrid to 2.5*2.5
lon = 0:359;
lat = 90:-1:-90;
lonf = 0:2.5:357.5; nlonf = length(lonf);
latf = 90:-2.5:-90; nlatf = length(latf);
time2 = 1:ntime; 
[Xlon, Ylat, Ttime] = meshgrid(lat,lon,time2);
[Xlonf, Ylatf, Ttimef] = meshgrid(latf,lonf,time2);
ps = interp3(Xlon,Ylat,Ttime,ps,Xlonf,Ylatf,Ttimef);
%calculate mean of months
ps_m=zeros(nlonf,nlatf,12);
[yy,mm,dd] = datevec(time);
for im = 1:12
    index_mm = mm == im;
    ps_m(:,:,im) = nanmean(ps(:,:,index_mm),3);%this is the average of each month
end
month = 1:12;

%% Save to nc
saveto = strcat(outpath,'ERAi_ps_month.nc');
ncid = netcdf.create(saveto,'NC_WRITE');
%Define the dimensions
dimidlon = netcdf.defDim(ncid,'longitude',nlonf);
dimidlat = netcdf.defDim(ncid,'latitude',nlatf);
dimidmonth = netcdf.defDim(ncid,'month',12);
%Define IDs for the dimension variables 
longitude_ID=netcdf.defVar(ncid,'longitude','double',dimidlon);
latitude_ID=netcdf.defVar(ncid,'latitude','double',dimidlat);
month_ID=netcdf.defVar(ncid,'month','double',dimidmonth);
%Define units for the dimension variables
netcdf.putAtt(ncid,longitude_ID,'standard_name','longitude');
netcdf.putAtt(ncid,latitude_ID,'standard_name','latitude');
netcdf.putAtt(ncid,month_ID,'long_name','month');
netcdf.putAtt(ncid,month_ID,'units','Month (JFMAMJJASOND)');

%Define the main variable ()
Ps_ID = netcdf.defVar(ncid,'Ps','double',[dimidlon dimidlat dimidmonth]);
%Define units for the main variables
netcdf.putAtt(ncid,Ps_ID,'long_name','Surface Pressure(monthly mean of 18 years)');
netcdf.putAtt(ncid,Ps_ID,'units', 'hPa');

%We are done defining the NetCdf
netcdf.endDef(ncid);

%Then store the dimension variables in
netcdf.putVar(ncid,longitude_ID,lonf);
netcdf.putVar(ncid,latitude_ID,latf);
netcdf.putVar(ncid,month_ID,month);

%Then store my main variable
netcdf.putVar(ncid,Ps_ID,ps_m);

%We're done, close the netcdf
netcdf.close(ncid)