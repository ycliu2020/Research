%% Combine the ERAi rad data from 2000-03 to 2018-02(18 years)
clc;clear; tic;pre_load;
%modify path first
filepath = 'H:\xieyan\ERAi\Rad_ori\';
testpath= 'H:\Repeat_the_experiment\testdata\ERAi\';
%creat 19yrs to cal
dR_lw_toa_clr0 = zeros(360,181,228);
dR_sw_toa_clr0 = zeros(360,181,228);
dR_lw_sfc_clr0 = zeros(360,181,228);
dR_sw_sfc_clr0 = zeros(360,181,228);
time0 = zeros(228,1);
for ii = 0:18
    filename = strcat('E:\xieyan\ERAi\Rad_ori\rad2004.nc');
    filename = strcat(filepath,'rad',num2str((ii+2000),'%04i'),'.nc');
    temp1 = ncread(filename,'ttrc');        % Top net thermal radiation, clear sky
    temp2 = ncread(filename,'tsrc');        % Top net solar radiation, clear sky
    temp3 = ncread(filename,'strc');        % Surface net thermal radiation, clear sky
    temp4 = ncread(filename,'ssrc');        % Surface net solar radiation, clear sky
    temp5 = ncread(filename,'time');       % hours since 1900-01-01 00:00:00.0
        % sum up the cumulative radiation of time step 0-12 and 12-24
        dR_lw_toa_clr0(:,:,ii*12+1:(ii+1)*12) = temp1(:,:,1:2:23);
        dR_sw_toa_clr0(:,:,ii*12+1:(ii+1)*12) = temp2(:,:,1:2:23);
        dR_lw_sfc_clr0(:,:,ii*12+1:(ii+1)*12) = temp3(:,:,1:2:23);
        dR_sw_sfc_clr0(:,:,ii*12+1:(ii+1)*12) = temp4(:,:,1:2:23);
        time0(ii*12+1:(ii+1)*12,1) = temp5(1:2:23);
end
%以下为掐头去尾2000年3月至2018年2月期间
dR_lw_toa_clr = dR_lw_toa_clr0(:,:,3:218); 
dR_sw_toa_clr = dR_sw_toa_clr0(:,:,3:218); 
dR_lw_sfc_clr = dR_lw_sfc_clr0(:,:,3:218); 
dR_sw_sfc_clr = dR_sw_sfc_clr0(:,:,3:218); 
time = time0(3:218,1);
time = datenum(1900,01,01,double(time),00,00);
ntime = length(time);

dR_lw_toa_clr = dR_lw_toa_clr./(3600*24); 
dR_sw_toa_clr = dR_sw_toa_clr./(3600*24); 
dR_lw_sfc_clr = dR_lw_sfc_clr./(3600*24); 
dR_sw_sfc_clr = dR_sw_sfc_clr./(3600*24); 
%% Regrid to 2.5*2.5
lon = 0:359;
lat = 90:-1:-90;
lonf = 0:2.5:357.5; nlonf = length(lonf);
latf = 88.75:-2.5:-88.75; nlatf = length(latf);
time2 = 1:ntime; 

[Xlon, Ylat, Ttime] = meshgrid(lat,lon,time2);
[Xlonf, Ylatf, Ttimef] = meshgrid(latf,lonf,time2);

dR_lw_toa_clr = interp3(Xlon,Ylat,Ttime,dR_lw_toa_clr,Xlonf,Ylatf,Ttimef);
dR_sw_toa_clr = interp3(Xlon,Ylat,Ttime,dR_sw_toa_clr,Xlonf,Ylatf,Ttimef);
dR_lw_sfc_clr = interp3(Xlon,Ylat,Ttime,dR_lw_sfc_clr,Xlonf,Ylatf,Ttimef);
dR_sw_sfc_clr = interp3(Xlon,Ylat,Ttime,dR_sw_sfc_clr,Xlonf,Ylatf,Ttimef);

% Deseasonalized(消除季节变动的数据即变量减去平均值)
[dR_lw_toa_clr,ClimdR_lw_toa_clr] = monthlyAnomaly3D(144,72,time,dR_lw_toa_clr);
[dR_sw_toa_clr,ClimdR_sw_toa_clr] = monthlyAnomaly3D(144,72,time,dR_sw_toa_clr);
[dR_lw_sfc_clr,ClimdR_lw_sfc_clr] = monthlyAnomaly3D(144,72,time,dR_lw_sfc_clr);
[dR_sw_sfc_clr,ClimdR_sw_sfc_clr] = monthlyAnomaly3D(144,72,time,dR_sw_sfc_clr);
month = 1:12;

%% Save to nc
saveto = strcat(testpath,'ERAi_radAnom.nc');
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
%Define units for the dimension variables
netcdf.putAtt(ncid,longitude_ID,'standard_name','longitude');
netcdf.putAtt(ncid,latitude_ID,'standard_name','latitude');
netcdf.putAtt(ncid,time_ID,'long_name','time');
netcdf.putAtt(ncid,time_ID,'units','days since 0000-00-00 00:00:0.0');
netcdf.putAtt(ncid,time_ID,'calendar','gregorian');
netcdf.putAtt(ncid,month_ID,'long_name','month');
netcdf.putAtt(ncid,month_ID,'units','Month (JFMAMJJASOND)');

%Define the main variable ()
dR_lw_toa_clr_ID = netcdf.defVar(ncid,'dR_lw_toa_clr','double',[dimidlon dimidlat dimidtime]);
dR_sw_toa_clr_ID = netcdf.defVar(ncid,'dR_sw_toa_clr','double',[dimidlon dimidlat dimidtime]);
dR_lw_sfc_clr_ID = netcdf.defVar(ncid,'dR_lw_sfc_clr','double',[dimidlon dimidlat dimidtime]);
dR_sw_sfc_clr_ID = netcdf.defVar(ncid,'dR_sw_sfc_clr','double',[dimidlon dimidlat dimidtime]);

ClimdR_lw_toa_clr_ID = netcdf.defVar(ncid,'ClimdR_lw_toa_clr','double',[dimidlon dimidlat dimidmonth]);
ClimdR_sw_toa_clr_ID = netcdf.defVar(ncid,'ClimdR_sw_toa_clr','double',[dimidlon dimidlat dimidmonth]);
ClimdR_lw_sfc_clr_ID = netcdf.defVar(ncid,'ClimdR_lw_sfc_clr','double',[dimidlon dimidlat dimidmonth]);
ClimdR_sw_sfc_clr_ID = netcdf.defVar(ncid,'ClimdR_sw_sfc_clr','double',[dimidlon dimidlat dimidmonth]);
%Define units for the main variables
netcdf.putAtt(ncid, dR_lw_toa_clr_ID,'long_name','Radiation Anomaly (lw) at TOA, clear sky');
netcdf.putAtt(ncid, dR_lw_toa_clr_ID,'units', 'W m**-2');
netcdf.putAtt(ncid, dR_sw_toa_clr_ID,'long_name','Radiation Anomaly (sw) at TOA, clear sky');
netcdf.putAtt(ncid, dR_sw_toa_clr_ID,'units', 'W m**-2');
netcdf.putAtt(ncid,dR_lw_sfc_clr_ID,'long_name','Radiation Anomaly (lw) at SFC, clear sky');
netcdf.putAtt(ncid,dR_lw_sfc_clr_ID,'units', 'W m**-2');
netcdf.putAtt(ncid,dR_sw_sfc_clr_ID,'long_name','Radiation Anomaly (sw) at SFC, clear sky');
netcdf.putAtt(ncid,dR_sw_sfc_clr_ID,'units', 'W m**-2');

netcdf.putAtt(ncid, ClimdR_lw_toa_clr_ID,'long_name','Radiation Climatology (lw) at TOA, clear sky');
netcdf.putAtt(ncid, ClimdR_lw_toa_clr_ID,'units', 'W m**-2');
netcdf.putAtt(ncid, ClimdR_sw_toa_clr_ID,'long_name','Radiation Climatology (sw) at TOA, clear sky');
netcdf.putAtt(ncid, ClimdR_sw_toa_clr_ID,'units', 'W m**-2');
netcdf.putAtt(ncid,ClimdR_lw_sfc_clr_ID,'long_name','Radiation Climatology (lw) at SFC, clear sky');
netcdf.putAtt(ncid,ClimdR_lw_sfc_clr_ID,'units', 'W m**-2');
netcdf.putAtt(ncid,ClimdR_sw_sfc_clr_ID,'long_name','Radiation Climatology (sw) at SFC, clear sky');
netcdf.putAtt(ncid,ClimdR_sw_sfc_clr_ID,'units', 'W m**-2');
%We are done defining the NetCdf
netcdf.endDef(ncid);

%Then store the dimension variables in
netcdf.putVar(ncid,longitude_ID,lonf);
netcdf.putVar(ncid,latitude_ID,latf);
netcdf.putVar(ncid,time_ID,time);
netcdf.putVar(ncid,month_ID,month);

%Then store my main variable
netcdf.putVar(ncid,dR_lw_toa_clr_ID,dR_lw_toa_clr);
netcdf.putVar(ncid,dR_sw_toa_clr_ID,dR_sw_toa_clr);
netcdf.putVar(ncid,dR_lw_sfc_clr_ID,dR_lw_sfc_clr);
netcdf.putVar(ncid,dR_sw_sfc_clr_ID,dR_sw_sfc_clr);

netcdf.putVar(ncid,ClimdR_lw_toa_clr_ID,ClimdR_lw_toa_clr);
netcdf.putVar(ncid,ClimdR_sw_toa_clr_ID,ClimdR_sw_toa_clr);
netcdf.putVar(ncid,ClimdR_lw_sfc_clr_ID,ClimdR_lw_sfc_clr);
netcdf.putVar(ncid,ClimdR_sw_sfc_clr_ID,ClimdR_sw_sfc_clr);
%We're done, close the netcdf
netcdf.close(ncid)
 
t=toc;disp(t)


