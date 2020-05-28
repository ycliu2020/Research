%% Combine the surface CERES cldrad and HadCRUT4 Ts data from 2000-03 to 2018-02 (18 years)
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

% -------section 1: calculate the radiation anomalies from CERES/EBAF4.0-----------
vector_var_str{1} = 'sfc';
% vector_var_str{2} = 'toa';

%modify path first
filepath = 'E:\Repeat_the_experiment\testdata\CERES_EBAF\';
outpath= 'E:\Repeat_the_experiment\testdata\CERES_EBAF\';
filename1 = 'CERES_EBAF-Surface_Ed4.1_Subset_200003-201802.nc';


var = char(vector_var_str{1});
filename = strcat(filepath, filename1);
% All-sky
dRnetCE = ncread(filename,'sfc_net_tot_all_mon');
dRnetCE_lw_cld = ncread(filename,'sfc_net_lw_all_mon');
dRdwnCE_lw_cld = ncread(filename,'sfc_lw_down_all_mon');
dRnetCE_sw_cld = ncread(filename,'sfc_net_sw_all_mon');
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
dRnetCE(1,:,:) = dRnetCE(361,:,:);
dRnetCE_lw_cld(1,:,:) = dRnetCE_lw_cld(361,:,:);
dRdwnCE_lw_cld(1,:,:) = dRdwnCE_lw_cld(361,:,:);
dRnetCE_sw_cld(1,:,:) = dRnetCE_sw_cld(361,:,:);

% dRnetCE(362,:,:) = dRnetCE(2,:,:);
% dRnetCE_lw_cld(362,:,:) = dRnetCE_lw_cld(2,:,:);
% dRdwnCE_lw_cld(362,:,:) = dRdwnCE_lw_cld(2,:,:);
% dRnetCE_sw_cld(362,:,:) = dRnetCE_sw_cld(2,:,:);

% Regraded to 2.5*2.5
time2 = 1:ntime;
[Xlon, Ylat, Ttime] = meshgrid(lat,lon,time2);
[Xlonf, Ylatf, Ttimef] = meshgrid(late,lone,time2);

dRnetCE = interp3(Xlon,Ylat,Ttime,dRnetCE,Xlonf,Ylatf,Ttimef);
dRnetCE_lw_cld = interp3(Xlon,Ylat,Ttime,dRnetCE_lw_cld,Xlonf,Ylatf,Ttimef);
dRdwnCE_lw_cld = interp3(Xlon,Ylat,Ttime,dRdwnCE_lw_cld,Xlonf,Ylatf,Ttimef);
dRnetCE_sw_cld = interp3(Xlon,Ylat,Ttime,dRnetCE_sw_cld,Xlonf,Ylatf,Ttimef);

%% Deseasonalized Anomaly
[dRnetCE,Clim_dRnet] = monthlyAnomaly3D(nlone,nlate,time,dRnetCE);
[dRnetCE_lw_cld,Clim_dRnet_lw] = monthlyAnomaly3D(nlone,nlate,time,dRnetCE_lw_cld);
[dRdwnCE_lw_cld,Clim_dRdwn_lw] = monthlyAnomaly3D(nlone,nlate,time,dRdwnCE_lw_cld);
[dRnetCE_sw_cld,Clim_dRnet_sw] = monthlyAnomaly3D(nlone,nlate,time,dRnetCE_sw_cld);

% -------section 2: calculate the ts anomalies from HadCRUT4-----------
filename2 = 'E:\xieyan\HadCRUT4\HadCRUT.4.6.0.0.median.nc';
ts0 = ncread(filename2,'temperature_anomaly');
dts0_had = ts0(:,:,1803:2018);    % 1803:2000.03, 2018:2018.02, 72*36*216 surface temperature anomaly, Global land
lonh = ncread(filename2,'longitude'); nlonh = length(lonh);                  % -177.5:5:177.5
lath = ncread(filename2,'latitude'); nlath = length(lath);
timeh = ncread(filename2,'time');  timeh = timeh(1803:2018);                 % 1979/03 - 2018/02
ntimeh = length(timeh);
timeh = timeh + datenum(1850,1,1);
timeh=double(timeh);
%缺测太多将在时间序列中只要缺测一次就剔除
test=dts0_had;
test1=sum(isnan(dts0_had),3);
for ii=1:216
    temp1=test(:,:,ii);
   temp1(test1~=0)=0;
   test(:,:,ii)=temp1;
end
dts0_had=test;

% -------section 3: calculate rhs(LW↓+SW) and ts anomalies' trend(mon,season,year)----------
% cal the rhs(LW↓+SW)
CErhs = dRdwnCE_lw_cld + dRnetCE_sw_cld;
% cal the ts trend
mon_trend_hadTs = zeros(nlonh, nlath, 12, 2); mon_trend_CErhs = zeros(nlone, nlate, 12, 2); % trend, 12 months
% im==1 is MAM, im==2 is JJA, im==3 is SON, im==4 is DJF
sea_trend_hadTs = zeros(nlonh, nlath, 4);sea_trend_CErhs = zeros(nlone, nlate, 4);% trend, 4 seasons

p_hadTs = zeros(nlonh, nlath, 12); p_CErhs = zeros(nlone, nlate, 12); % f test
cons_hadTs = zeros(nlonh, nlath, 12); cons_CErhs = zeros(nlone, nlate, 12); % mean
% months mean trend(unite:per day)
for i = 1:nlone
    for j = 1:nlate
        [~, mon_trend_CErhs(i, j, :, :), cons_CErhs(i, j, :), p_CErhs(i, j, :)] = detrend_yan(CErhs(i, j, :), time);
    end
end

for i = 1:nlonh
    for j = 1:nlath
        [~, mon_trend_hadTs(i, j, :, :), cons_hadTs(i, j, :), p_hadTs(i, j, :)] = detrend_yan(dts0_had(i, j, :), timeh);
    end
end

% seasons mean trend (unite:per day)(1-4代表3-5,.....)
for ii = 1:4
    if ii ~= 4
        sea_trend_hadTs(:,:,ii) = squeeze(mean(mon_trend_hadTs(:,:,ii*3:ii*3+2,1),3));
        sea_trend_CErhs(:,:,ii) = squeeze(mean(mon_trend_CErhs(:,:,ii*3:ii*3+2,1),3));
    elseif ii == 4
        ii1=[1,2,12];
        sea_trend_hadTs(:,:,ii) = squeeze(mean(mon_trend_hadTs(:,:,ii1,1),3));
        sea_trend_CErhs(:,:,ii) = squeeze(mean(mon_trend_CErhs(:,:,ii1,1),3));
    end
end
% years mean trend (unite:per day)
yr_trend_hadTs = squeeze(mean(mon_trend_hadTs(:,:,:,1),3));
yr_trend_CErhs = squeeze(mean(mon_trend_CErhs(:,:,:,1),3));

% ----save-----
outfilename1=[outpath,'CERESrad_HadTs4.1.mat'];
save(outfilename1,'lone','late','dRnetCE','dRnetCE_lw_cld','dRdwnCE_lw_cld','dRnetCE_sw_cld','Clim_dRnet','Clim_dRnet_lw','Clim_dRdwn_lw','Clim_dRnet_sw','time',...
    'lonh','lath','dts0_had','timeh')

outfilename2=[outpath,'trend_CERESrad_HadTs4.1.mat'];
save(outfilename2,'lone','late','lonh','lath','time','timeh',...
    'mon_trend_hadTs','mon_trend_CErhs','sea_trend_hadTs','sea_trend_CErhs','yr_trend_hadTs','yr_trend_CErhs',...
    'p_hadTs','p_CErhs','cons_hadTs','cons_CErhs')

t=toc;disp(t)




