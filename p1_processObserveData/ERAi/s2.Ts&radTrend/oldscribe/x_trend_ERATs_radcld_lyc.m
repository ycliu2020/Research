% x_trend_Ts_radcld_lyc.m
% cal the trend, if figure show the result is good, then cal the the coef
% note monthly data from 2000-03 to 2018-02
% trend1:(12,2) k and b, from 1-12 month
% cons: (12,1) means of month in 18 years
% pvalue: (12,1) f检�? 检验结�?
% time seris:
% 1.every month
% 2.MAM,JJA,SON,DJF
% 3.year

clc; clear; tic; pre_load;
inputfile = '/home/lyc/repeat_Kernel_experiment/testdata/ERAi/ERAi_radcld_ts_Anom_mon.nc';
time = ncread(inputfile, 'time'); ntime = length(time);
lon = ncread(inputfile, 'longitude'); nlon = length(lon);
lat = ncread(inputfile, 'latitude'); nlat = length(lat);

dts = ncread(inputfile, 'dts_mon_sfc'); % skin temperature K
lwdow_anom = ncread(inputfile, 'dR_mon_downlw_sfc_cld'); % Radiation Anomaly (downwards lw) at SFC W/m*2
swnet_anom = ncread(inputfile, 'dR_mon_netsw_sfc_cld'); % Radiation Anomaly (net sw) at SFC
drhs = lwdow_anom + swnet_anom; % R.H.S = lwdow_anom + swnet_anom

trendm_dts = zeros(nlon, nlat, 12, 2); trendm_drhs = zeros(nlon, nlat, 12, 2); % trend, 12 months
% im==1 is MAM, im==2 is JJA, im==3 is SON, im==4 is DJF
trends_dts = zeros(nlon, nlat, 4);trends_drhs = zeros(nlon, nlat, 4);% trend, 4 seasons
% trendyr_dts = zeros(nlon, nlat); trendyr_drhs = zeros(nlon, nlat); % trend, yearly

p_ts = zeros(nlon, nlat, 12); p_rhs = zeros(nlon, nlat, 12); % f test
cons_ts = zeros(nlon, nlat, 12); cons_rhs = zeros(nlon, nlat, 12); % mean
% months mean trend(unite:per day)
for i = 1:nlon
    for j = 1:nlat
        [~, trendm_dts(i, j, :, :), cons_ts(i, j, :), p_ts(i, j, :)] = detrend_yan(dts(i, j, :), time);
        [~, trendm_drhs(i, j, :, :), cons_rhs(i, j, :), p_rhs(i, j, :)] = detrend_yan(drhs(i, j, :), time);
    end
end
% seasons mean trend (unite:per day)
for ii = 1:4
 if ii ~= 4
  trends_dts(:,:,ii) = squeeze(mean(trendm_dts(:,:,ii*3:ii*3+2,1),3));
  trends_drhs(:,:,ii) = squeeze(mean(trendm_drhs(:,:,ii*3:ii*3+2,1),3));
 elseif ii == 4
    ii1=[1,2,12];
    trends_dts(:,:,ii) = squeeze(mean(trendm_dts(:,:,ii1,1),3));
    trends_drhs(:,:,ii) = squeeze(mean(trendm_drhs(:,:,ii1,1),3));
 end
end
% years mean trend (unite:per day)
trendyr_dts = squeeze(mean(trendm_dts(:,:,:,1),3));
trendyr_drhs = squeeze(mean(trendm_drhs(:,:,:,1),3));


%----------------------save varity----------------------------------
outfilename = '/home/lyc/repeat_Kernel_experiment/testdata/trend_tsrad';
save(outfilename, 'lon', 'lat', 'time', 'trendm_dts', 'trendm_drhs', 'cons_ts', 'cons_rhs', 'p_ts', 'p_rhs','trends_dts','trends_drhs',...
'trendyr_dts','trendyr_drhs')

t = toc; disp(t)