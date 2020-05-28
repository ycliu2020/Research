%% Combine the ERAi climatological data from 2000-03 to 2018-02
clc;clear; tic;
%modify path first
filepath1 = '/home/lyc/data/era_download/water_vapor/';
tcwv0 =zeros(360,181,228);      %Total column water vapour

time0 = zeros(228,1);
for ii = 0:18
    filename1 = strcat(filepath1,'wv',num2str((ii+2000),'%04i'),'.nc'); 
    temp1 = ncread(filename1,'tcwv');
    temp5 = ncread(filename1,'time');
    tcwv0 (:,:,ii*12+1:(ii+1)*12) = temp1;
    time0(ii*12+1:(ii+1)*12,1) = temp5;
end
tcwv = tcwv0(:,:,3:218);
time = time0(3:218,1); 
time = time./24;
time = datenum(1900,01,01)+time; ntime = length(time);

%% Regrid to 2.5*2.5
lon = 0:359; lat = 90:-1:-90;
% time2 = 1:348;                                                   %29*12 month
time2 = 1:216;                                                       %18*12 month
    
lonf = 0:2.5:357.5; nlonf = length(lonf);
latf = 90:-2.5:-90; nlatf = length(latf);

[Xlon, Ylat, Ttime] = meshgrid(lat,lon,time2);
[Xlonf, Ylatf, Ttimef] = meshgrid(latf,lonf,time2);
tcwv = interp3(Xlon,Ylat,Ttime,tcwv,Xlonf,Ylatf,Ttimef);

%% Deseasonalize
[Anom_tcwv,Clim_tcwv] = monthlyAnomaly3D(144,73,time,tcwv);

outfilename=[filepath1,'wv_variAnom2018.mat'];

save(outfilename,'lon','lat','tcwv','Anom_tcwv','Clim_tcwv','time')
t=toc;
disp(t)







