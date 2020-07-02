%% Combine the ERAi rad and Ts data from 2000-03 to 2018-02 (18 years)
% attention:  matlab's calendar doesnt consist with the nc file, so result nc show's wrong time in ncview, but doest affect result
% only matter that data from 2000-03 to 2018-02 (18 years), not from 2000-01 to 2017-12
% 
% this function only cal surface varity anomaly (including cloud and clc):
% strd(Surface thermal radiation downwards), 
% ssr(surface_net_downward_shortwave_flux),
% ts
% cal time seris:
% every month
% 

clc;clear; tic;pre_load;
%modify path first
filepath1 = '/data1/xieyan/ERAi/Rad_ori/';
filepath2 = '/data1/xieyan/ERAi/Heat_Flux/Heat_Flux_Accumlate/';
filepath3 = '/data1/xieyan/ERAi/Ts/';
outpath= '/home/lyc/repeat_Kernel_experiment/testdata/';
%creat 19yrs to cal
dR_downlw_sfc_cld0 = zeros(360,181,228);
dts_sfc0 = zeros(360,181,228);
dR_netsw_sfc_cld0 = zeros(360,181,228);
time0 = zeros(228,1);
for ii = 0:18
    % filename1 = strcat('E:\xieyan\ERAi\Rad_ori\rad2004.nc');
    filename1 = strcat(filepath1,'rad',num2str((ii+2000),'%04i'),'.nc');
    filename2 = strcat(filepath2,'hf',num2str((ii+2000),'%04i'),'.nc');
    filename3 = strcat(filepath3,'t',num2str((ii+2000),'%04i'),'.nc');
    temp1 = ncread(filename2,'strd');        % 
    temp2 = ncread(filename3,'skt');        % ts
    temp3 = ncread(filename1,'ssr');        % 
    temp4 = ncread(filename1,'time');       % hours since 1900-01-01 00:00:00.0
        % sum up the cumulative radiation of time step 0-12 and 12-24
        dR_downlw_sfc_cld0(:,:,ii*12+1:(ii+1)*12) = temp1(:,:,1:2:23);
        dts_sfc0(:,:,ii*12+1:(ii+1)*12) = temp2(:,:,:);
        dR_netsw_sfc_cld0(:,:,ii*12+1:(ii+1)*12) = temp3(:,:,1:2:23);
        time0(ii*12+1:(ii+1)*12,1) = temp4(1:2:23);
end
%???????2000?3??2018?2???
dR_mon_downlw_sfc_cld = dR_downlw_sfc_cld0(:,:,3:218); 
dts_mon_sfc = dts_sfc0(:,:,3:218); 
dR_mon_netsw_sfc_cld = dR_netsw_sfc_cld0(:,:,3:218); 
time = time0(3:218,1);
time = datenum(1900,01,01,double(time),00,00);
ntime = length(time);

dR_mon_downlw_sfc_cld = dR_mon_downlw_sfc_cld./(3600*24); 
dR_mon_netsw_sfc_cld = dR_mon_netsw_sfc_cld./(3600*24); 
%% Regrid to 2.5*2.5
lon = 0:359;
lat = 90:-1:-90;
lonf = 0:2.5:357.5; nlonf = length(lonf);
latf = 88.75:-2.5:-88.75; nlatf = length(latf);
time2 = 1:ntime; 

[Xlon, Ylat, Ttime] = meshgrid(lat,lon,time2);
[Xlonf, Ylatf, Ttimef] = meshgrid(latf,lonf,time2);

dR_mon_downlw_sfc_cld = interp3(Xlon,Ylat,Ttime,dR_mon_downlw_sfc_cld,Xlonf,Ylatf,Ttimef);
dts_mon_sfc = interp3(Xlon,Ylat,Ttime,dts_mon_sfc,Xlonf,Ylatf,Ttimef);
dR_mon_netsw_sfc_cld = interp3(Xlon,Ylat,Ttime,dR_mon_netsw_sfc_cld,Xlonf,Ylatf,Ttimef);

% Deseasonalized(????????)
[dR_mon_downlw_sfc_cld,ClimdR_mon_downlw_sfc_cld] = monthlyAnomaly3D(144,72,time,dR_mon_downlw_sfc_cld);
[dts_mon_sfc,Climdts_mon_sfc] = monthlyAnomaly3D(144,72,time,dts_mon_sfc);
[dR_mon_netsw_sfc_cld,ClimdR_mon_netsw_sfc_cld] = monthlyAnomaly3D(144,72,time,dR_mon_netsw_sfc_cld);
month = 1:12;



%% Save to nc
saveto = strcat(outpath,'ERAi_radcld_ts_Anom_mon.nc');
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
dR_downlw_sfc_cld_ID = netcdf.defVar(ncid,'dR_mon_downlw_sfc_cld','double',[dimidlon dimidlat dimidtime]);
dTs_sfc_ID = netcdf.defVar(ncid,'dts_mon_sfc','double',[dimidlon dimidlat dimidtime]);
dR_netsw_sfc_cld_ID = netcdf.defVar(ncid,'dR_mon_netsw_sfc_cld','double',[dimidlon dimidlat dimidtime]);

ClimdR_downlw_sfc_cld_ID = netcdf.defVar(ncid,'ClimdR_mon_downlw_sfc_cld','double',[dimidlon dimidlat dimidmonth]);
ClimdTs_mon_sfc_ID = netcdf.defVar(ncid,'Climdts_mon_sfc','double',[dimidlon dimidlat dimidmonth]);
ClimdR_mon_netsw_sfc_cld_ID = netcdf.defVar(ncid,'ClimdR_mon_netsw_sfc_cld','double',[dimidlon dimidlat dimidmonth]);
%Define units for the main variables
netcdf.putAtt(ncid, dR_downlw_sfc_cld_ID,'long_name','Radiation Anomaly (downwards lw) at SFC');
netcdf.putAtt(ncid, dR_downlw_sfc_cld_ID,'units', 'W m**-2');
netcdf.putAtt(ncid,dTs_sfc_ID,'long_name','Skin temperature');
netcdf.putAtt(ncid,dTs_sfc_ID,'units', 'K');
netcdf.putAtt(ncid,dR_netsw_sfc_cld_ID,'long_name','Radiation Anomaly (net sw) at SFC');
netcdf.putAtt(ncid,dR_netsw_sfc_cld_ID,'units', 'W m**-2');

netcdf.putAtt(ncid, ClimdR_downlw_sfc_cld_ID,'long_name','Radiation Climatology (downwards lw) at SFC');
netcdf.putAtt(ncid, ClimdR_downlw_sfc_cld_ID,'units', 'W m**-2');
netcdf.putAtt(ncid,ClimdTs_mon_sfc_ID,'long_name','Ts Climatology at SFC');
netcdf.putAtt(ncid,ClimdTs_mon_sfc_ID,'units', 'K');
netcdf.putAtt(ncid,ClimdR_mon_netsw_sfc_cld_ID,'long_name','Radiation Climatology (net sw) at SFC');
netcdf.putAtt(ncid,ClimdR_mon_netsw_sfc_cld_ID,'units', 'W m**-2');
%We are done defining the NetCdf
netcdf.endDef(ncid);

%Then store the dimension variables in
netcdf.putVar(ncid,longitude_ID,lonf);
netcdf.putVar(ncid,latitude_ID,latf);
netcdf.putVar(ncid,time_ID,time);
netcdf.putVar(ncid,month_ID,month);

%Then store my main variable
netcdf.putVar(ncid,dR_downlw_sfc_cld_ID,dR_mon_downlw_sfc_cld);
netcdf.putVar(ncid,dTs_sfc_ID,dts_mon_sfc);
netcdf.putVar(ncid,dR_netsw_sfc_cld_ID,dR_mon_netsw_sfc_cld);

netcdf.putVar(ncid,ClimdR_downlw_sfc_cld_ID,ClimdR_mon_downlw_sfc_cld);
netcdf.putVar(ncid,ClimdTs_mon_sfc_ID,Climdts_mon_sfc);
netcdf.putVar(ncid,ClimdR_mon_netsw_sfc_cld_ID,ClimdR_mon_netsw_sfc_cld);
%We're done, close the netcdf
netcdf.close(ncid)
 
t=toc;disp(t)


