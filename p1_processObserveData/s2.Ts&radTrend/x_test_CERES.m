% plot CERES's mean rad to compare the long downward and net short rad why is different. 

clc; clear; tic; pre_load;
load('D:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\tools\map\a_world_map\mask\mask_cp144.mat') % load word land mask
load('D:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\tools\map\a_world_map\mask\mask_ce72.mat') % load word land mask
load('D:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\tools\map\a_world_map\correct_worldmap.mat') % 订正后的世界地图文件word_mapx(:),word_mapy(:)
load('D:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\tools\map\a_china_map\mask14472.mat')
% CERES rad
filepath = 'E:\Repeat_the_experiment\testdata\CERES_EBAF\';
filename1 = 'CERES_EBAF-Surface_Ed4.1_Subset_200003-201802.nc';
filename = strcat(filepath, filename1);

% All-sky
dRnetCE = ncread(filename,'sfc_net_tot_all_mon');dRnetCE_lw_cld = ncread(filename,'sfc_net_lw_all_mon');dRdwnCE_lw_cld = ncread(filename,'sfc_lw_down_all_mon');dRnetCE_sw_cld = ncread(filename,'sfc_net_sw_all_mon');
lon = -0.5:1:359.5; nlon = length(lon); lon = double(lon);
lat = -89.5:1:89.5; nlat = length(lat); lat = double(lat);
time = double(ncread(filename,'time'));                            % days since 2000-03-01 00:00:00

time = datenum(2000,03,1 + double(time),00,00,00);
ntime = length(time);
month = 1:12;

%% Regrid (Because of kernels)
lone = 0:2.5:357.5; nlone = length(lone);
late = 88.75:-2.5:-88.75; nlate = length(late);
% For interpolate(此步作用是将lonf变成lon的子集(lon范围扩大), latf原来就是lat的子集所以不用扩充)
for ll = 360:-1:1
    dRnetCE(ll+1,:,:) = dRnetCE(ll,:,:);
    dRnetCE_lw_cld(ll+1,:,:) = dRnetCE_lw_cld(ll,:,:);
    dRdwnCE_lw_cld(ll+1,:,:) = dRdwnCE_lw_cld(ll,:,:);
    dRnetCE_sw_cld(ll+1,:,:) = dRnetCE_sw_cld(ll,:,:);
end
dRnetCE(1,:,:) = dRnetCE(361,:,:);dRnetCE_lw_cld(1,:,:) = dRnetCE_lw_cld(361,:,:);dRdwnCE_lw_cld(1,:,:) = dRdwnCE_lw_cld(361,:,:);dRnetCE_sw_cld(1,:,:) = dRnetCE_sw_cld(361,:,:);

% Regraded to 2.5*2.5(144*72*216)
time2 = 1:ntime;
[Xlon, Ylat, Ttime] = meshgrid(lat,lon,time2);
[Xlonf, Ylatf, Ttimef] = meshgrid(late,lone,time2);
dRnetCE = interp3(Xlon,Ylat,Ttime,dRnetCE,Xlonf,Ylatf,Ttimef);
dRnetCE_lw_cld = interp3(Xlon,Ylat,Ttime,dRnetCE_lw_cld,Xlonf,Ylatf,Ttimef);
dRdwnCE_lw_cld = interp3(Xlon,Ylat,Ttime,dRdwnCE_lw_cld,Xlonf,Ylatf,Ttimef);
dRnetCE_sw_cld = interp3(Xlon,Ylat,Ttime,dRnetCE_sw_cld,Xlonf,Ylatf,Ttimef);
[dRnetCE,Clim_dRnet] = monthlyAnomaly3D(nlone,nlate,time,dRnetCE);
[dRnetCE_lw_cld,Clim_dRnet_lw] = monthlyAnomaly3D(nlone,nlate,time,dRnetCE_lw_cld);
[dRdwnCE_lw_cld,Clim_dRdwn_lw] = monthlyAnomaly3D(nlone,nlate,time,dRdwnCE_lw_cld);
[dRnetCE_sw_cld,Clim_dRnet_sw] = monthlyAnomaly3D(nlone,nlate,time,dRnetCE_sw_cld);


mean_RdwnCE_lw_cld=nanmean(dRdwnCE_lw_cld,3)*365*10;
mean_RnetCE_sw_cld=nanmean(dRnetCE_sw_cld,3)*365*10;


% ERA rad
filepath1 = 'E:\xieyan\ERAi\Rad_ori\';
filepath2 = 'E:\xieyan\ERAi\Heat_Flux\Heat_Flux_Accumlate\';
%creat 19yrs to cal
dRdwnERA_lw_cld0 = zeros(360,181,228);
dRnetERA_sw_cld0 = zeros(360,181,228);
time0 = zeros(228,1);
for ii = 0:18
    % filename1 = strcat('E:\xieyan\ERAi\Rad_ori\rad2004.nc');
    filename1 = strcat(filepath1,'rad',num2str((ii+2000),'%04i'),'.nc');
    filename2 = strcat(filepath2,'hf',num2str((ii+2000),'%04i'),'.nc');
    temp1 = ncread(filename2,'strd');        % 
    temp3 = ncread(filename1,'ssr');        % 
    temp4 = ncread(filename1,'time');       % hours since 1900-01-01 00:00:00.0
        % sum up the cumulative radiation of time step 0-12 and 12-24
        dRdwnERA_lw_cld0(:,:,ii*12+1:(ii+1)*12) = temp1(:,:,1:2:23);
        dRnetERA_sw_cld0(:,:,ii*12+1:(ii+1)*12) = temp3(:,:,1:2:23);
        time0(ii*12+1:(ii+1)*12,1) = temp4(1:2:23);
end
%???????2000?3??2018?2???
dRdwnERA_lw_cld = dRdwnERA_lw_cld0(:,:,3:218); 
dRnetERA_sw_cld = dRnetERA_sw_cld0(:,:,3:218); 
time_era = time0(3:218,1);
time_era = datenum(1900,01,01,double(time_era),00,00);
ntime_era = length(time_era);
dRdwnERA_lw_cld = dRdwnERA_lw_cld./(3600*24); 
dRnetERA_sw_cld = dRnetERA_sw_cld./(3600*24); 
%% Regrid to 2.5*2.5
time2 = 1:ntime_era; 
lon = 0:359;
lat = 90:-1:-90;
[Xlon, Ylat, Ttime] = meshgrid(lat,lon,time2);
[Xlonf, Ylatf, Ttimef] = meshgrid(late,lone,time2);
dRdwnERA_lw_cld = interp3(Xlon,Ylat,Ttime,dRdwnERA_lw_cld,Xlonf,Ylatf,Ttimef);
dRnetERA_sw_cld = interp3(Xlon,Ylat,Ttime,dRnetERA_sw_cld,Xlonf,Ylatf,Ttimef);

[dRdwnERA_lw_cld,ClimdRdwnERA_lw_cld] = monthlyAnomaly3D(144,72,time,dRdwnERA_lw_cld);
[dRnetERA_sw_cld,ClimdRnetERA_sw_cld] = monthlyAnomaly3D(144,72,time,dRnetERA_sw_cld);

mean_RdwnERA_lw_cld = nanmean(dRdwnERA_lw_cld,3)*365*10; 
mean_RnetERA_sw_cld = nanmean(dRnetERA_sw_cld,3)*365*10; 


mean_all(:,:,1,1)=mean_RdwnERA_lw_cld;
mean_all(:,:,1,2)=mean_RdwnCE_lw_cld;
mean_all(:,:,2,1)=mean_RnetERA_sw_cld;
mean_all(:,:,2,2)=mean_RnetCE_sw_cld;
title_name{1,1}={'RdwnERA lw'};title_name{1,2}={'RdwnCE lw'};
title_name{2,1}={'RnetERA sw'};title_name{2,2}={'RnetCE sw'};
t=toc;disp(t)
