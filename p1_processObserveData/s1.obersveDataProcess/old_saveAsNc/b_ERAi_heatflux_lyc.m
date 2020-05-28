%% Combine the ERAi Heat Flux data from 2000-03 to 2018-02
% Include Latent Heat and Sensible Heat
clear
filepath = '/data1/xieyan/ERAi/Heat_Flux/Heat_Flux_Accumlate/';
lahf0 = zeros(360,181,228);
sehf0 = zeros(360,181,228);
slwr0 = zeros(360,181,228);
slwrd0 = zeros(360,181,228);
time0 = zeros(228,1);
for ii = 0:18
    filename = strcat(filepath,'hf',num2str((ii+2000),'%04i'),'.nc');
%     temp1 = ncread(filename,'slhf');                                      % surface_upward_latent_heat_flux
%     temp2 = ncread(filename,'sshf');                                      % surface_upward_sensible_heat_flux
%     temp3 = ncread(filename,'str');                                       % surface_net_upward_longwave_flux
%     temp4 = ncread(filename,'strd');                                      % Surface thermal radiation downwards
    temp5 = ncread(filename,'time');                                      % hours since 1900-01-01 00:00:00.0
        % sum up the cumulative radiation of time step 0-12 and 12-24
%         lahf0(:,:,ii*12+1:(ii+1)*12) = temp1(:,:,1:2:23);
%         sehf0(:,:,ii*12+1:(ii+1)*12) = temp2(:,:,1:2:23);
%         slwr0(:,:,ii*12+1:(ii+1)*12) = temp3(:,:,1:2:23);
%         slwrd0(:,:,ii*12+1:(ii+1)*12) = temp4(:,:,1:2:23);
        time0(ii*12+1:(ii+1)*12,1) = temp5(1:2:23);
end
% slwru0 = slwr0 + slwrd0;                                                  % surface upward longwave flux
% lahf = -lahf0(:,:,3:218);                                                 % define downward as positive
% sehf = -sehf0(:,:,3:218);                                                 % define downward as positive
% slwrd = slwrd0(:,:,3:218); 
% slwr = - slwr0(:,:,3:218);                                                % define downward as positive
% slwru = -slwru0(:,:,3:218);                                               % define downward as positive
time = time0(3:218,1);
time = datenum(1900,01,01,double(time),00,00);
ntime = length(time);

lahf = lahf./(3600*24); 
sehf = sehf./(3600*24); 
slwr = slwr./(3600*24); 
slwrd  = slwrd./(3600*24); 
slwru = slwru./(3600*24); 
%% Regrid to 2.5*2.5
lon = 0:1:359; lon = double(lon); nlon = length(lon);
lat = 90:-1:-90; lat = double(lat); nlat = length(lat);
lonf = 0:2.5:357.5; nlonf = length(lonf);
latf = 88.75:-2.5:-88.75; nlatf = length(latf);
time2 = 1:ntime; 

[Xlon,Ylat,Ttime] = meshgrid(lat,lon,time2);
[Xlonf,Ylatf,Ttimef] = meshgrid(latf,lonf,time2);

lahf = interp3(Xlon,Ylat,Ttime,lahf,Xlonf,Ylatf,Ttimef);
sehf = interp3(Xlon,Ylat,Ttime,sehf,Xlonf,Ylatf,Ttimef);
slwr = interp3(Xlon,Ylat,Ttime,slwr,Xlonf,Ylatf,Ttimef);
slwrd = interp3(Xlon,Ylat,Ttime,slwrd,Xlonf,Ylatf,Ttimef);
slwru = interp3(Xlon,Ylat,Ttime,slwru,Xlonf,Ylatf,Ttimef);
% Deseasonalized
[lahf,Clim_lahf] = monthlyAnomaly3D(144,72,time,lahf);
[sehf,Clim_sehf] = monthlyAnomaly3D(144,72,time,sehf);
[slwr,Clim_slwr] = monthlyAnomaly3D(144,72,time,slwr);
[slwrd,Clim_slwrd] = monthlyAnomaly3D(144,72,time,slwrd);
[slwru,Clim_slwru] = monthlyAnomaly3D(144,72,time,slwru);
month = 1:12;


%% Save to nc
filepath1='/home/lyc/repeat_Kernel_experiment/testdata/ERAi/';
saveto = strcat(filepath1,'ERAi_hfAnom.nc');
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
%Define the main variable ()
lahf_ID = netcdf.defVar(ncid,'lahf','double',[dimidlon dimidlat dimidtime]);
sehf_ID = netcdf.defVar(ncid,'sehf','double',[dimidlon dimidlat dimidtime]);
slwr_ID = netcdf.defVar(ncid,'slwr','double',[dimidlon dimidlat dimidtime]);
slwrd_ID = netcdf.defVar(ncid,'slwrd','double',[dimidlon dimidlat dimidtime]);
slwru_ID = netcdf.defVar(ncid,'slwru','double',[dimidlon dimidlat dimidtime]);
Clim_lahf_ID = netcdf.defVar(ncid,'Clim_lahf','double',[dimidlon dimidlat dimidmonth]);
Clim_sehf_ID = netcdf.defVar(ncid,'Clim_sehf','double',[dimidlon dimidlat dimidmonth]);
Clim_slwr_ID = netcdf.defVar(ncid,'Clim_slwr','double',[dimidlon dimidlat dimidmonth]);
Clim_slwrd_ID = netcdf.defVar(ncid,'Clim_slwrd','double',[dimidlon dimidlat dimidmonth]);
Clim_slwru_ID = netcdf.defVar(ncid,'Clim_slwru','double',[dimidlon dimidlat dimidmonth]);

%Define units for the dimension variables
netcdf.putAtt(ncid,longitude_ID,'standard_name','longitude');
netcdf.putAtt(ncid,latitude_ID,'standard_name','latitude');
netcdf.putAtt(ncid,time_ID,'long_name','time');
netcdf.putAtt(ncid,time_ID,'units','days since 0000-00-00 00:00:0.0');
netcdf.putAtt(ncid,time_ID,'calendar','gregorian');
netcdf.putAtt(ncid,month_ID,'long_name','month');
netcdf.putAtt(ncid,month_ID,'units','Month (JFMAMJJASOND)');

netcdf.putAtt(ncid,lahf_ID,'long_name','Surface latent heat flux, Anomaly');
netcdf.putAtt(ncid,lahf_ID,'units', 'W m**-2');
netcdf.putAtt(ncid,sehf_ID,'long_name','Surface sensible heat flux, Anomaly');
netcdf.putAtt(ncid,sehf_ID,'units', 'W m**-2');
netcdf.putAtt(ncid,slwr_ID,'long_name','Surface LW Radiation Anomaly, net');
netcdf.putAtt(ncid,slwr_ID,'units', 'W m**-2');
netcdf.putAtt(ncid,slwrd_ID,'long_name','Surface LW Radiation Anomaly, downward');
netcdf.putAtt(ncid,slwrd_ID,'units', 'W m**-2');
netcdf.putAtt(ncid,slwru_ID,'long_name','Surface LW Radiation Anomaly, upward');
netcdf.putAtt(ncid,slwru_ID,'units', 'W m**-2');

netcdf.putAtt(ncid,Clim_lahf_ID,'long_name','Surface latent heat flux, Climatology');
netcdf.putAtt(ncid,Clim_lahf_ID,'units', 'W m**-2');
netcdf.putAtt(ncid,Clim_sehf_ID,'long_name','Surface sensible heat flux, Climatology');
netcdf.putAtt(ncid,Clim_sehf_ID,'units', 'W m**-2');
netcdf.putAtt(ncid,Clim_slwr_ID,'long_name','Surface LW Radiation Climatology, net');
netcdf.putAtt(ncid,Clim_slwr_ID,'units', 'W m**-2');
netcdf.putAtt(ncid,Clim_slwrd_ID,'long_name','Surface LW Radiation Climatology, downward');
netcdf.putAtt(ncid,Clim_slwrd_ID,'units', 'W m**-2');
netcdf.putAtt(ncid,Clim_slwru_ID,'long_name','Surface LW Radiation Climatology, upward');
netcdf.putAtt(ncid,Clim_slwru_ID,'units', 'W m**-2');
%We are done defining the NetCdf
netcdf.endDef(ncid);

%Then store the dimension variables in
netcdf.putVar(ncid,longitude_ID,lonf);
netcdf.putVar(ncid,latitude_ID,latf);
netcdf.putVar(ncid,time_ID,time);
netcdf.putVar(ncid,month_ID,month);

%Then store my main variable
netcdf.putVar(ncid,lahf_ID,lahf);
netcdf.putVar(ncid,sehf_ID,sehf);
netcdf.putVar(ncid,slwr_ID,slwr);
netcdf.putVar(ncid,slwrd_ID,slwrd);
netcdf.putVar(ncid,slwru_ID,slwru);
netcdf.putVar(ncid,Clim_lahf_ID,Clim_lahf);
netcdf.putVar(ncid,Clim_sehf_ID,Clim_sehf);
netcdf.putVar(ncid,Clim_slwr_ID,Clim_slwr);
netcdf.putVar(ncid,Clim_slwrd_ID,Clim_slwrd);
netcdf.putVar(ncid,Clim_slwru_ID,Clim_slwru);
%We're done, close the netcdf
netcdf.close(ncid)


