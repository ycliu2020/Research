% x_trend_Ts_radcld_lyc.m
% cal the trend, if figure show the result is good, then cal the the coef
% note monthly data from 2000-03 to 2018-02
% trend1:(12,2) k and b, from 1-12 month
% cons: (12,1) means of month in 18 years
% pvalue: (12,1) f检验 检验结果
% time seris:
% 1.every month
% 2.MAM,JJA,SON,DJF
% 3.year

clc; clear; tic; pre_load;
inputfile = '/home/lyc/repeat_Kernel_experiment/testdata/ERAi/ERAi_radcld_ts_Anom_mon.nc';
time = ncread(inputfile, 'time'); ntime = length(time);
lon = ncread(inputfile, 'longitude'); nlon = length(lon);
lat = ncread(inputfile, 'latitude'); nlat = length(lat);

ts_anom = ncread(inputfile, 'dts_mon_sfc'); % skin temperature K
lwdow_anom = ncread(inputfile, 'dR_mon_downlw_sfc_cld'); % Radiation Anomaly (downwards lw) at SFC W/m*2
swnet_anom = ncread(inputfile, 'dR_mon_netsw_sfc_cld'); % Radiation Anomaly (net sw) at SFC
sumrhs = lwdow_anom + swnet_anom; % R.H.S = lwdow_anom + swnet_anom

mon_trend_Ts = zeros(nlon, nlat, 12, 2); mon_trend_rhs = zeros(nlon, nlat, 12, 2); % trend, 12 months
% im==1 is MAM, im==2 is JJA, im==3 is SON, im==4 is DJF
sea_trend_Ts = zeros(nlon, nlat, 4);sea_trend_rhs = zeros(nlon, nlat, 4);% trend, 4 seasons
% yr_trend_Ts = zeros(nlon, nlat); yr_trend_rhs = zeros(nlon, nlat); % trend, yearly

p_Ts = zeros(nlon, nlat, 12); p_rhs = zeros(nlon, nlat, 12); % f test
cons_Ts = zeros(nlon, nlat, 12); cons_rhs = zeros(nlon, nlat, 12); % mean
% months mean trend(unite:per day)
for i = 1:nlon

    for j = 1:nlat
        [~, mon_trend_Ts(i, j, :, :), cons_Ts(i, j, :), p_Ts(i, j, :)] = detrend_yan(ts_anom(i, j, :), time);
        [~, mon_trend_rhs(i, j, :, :), cons_rhs(i, j, :), p_rhs(i, j, :)] = detrend_yan(sumrhs(i, j, :), time);
    end

end
% seasons mean trend (unite:per day)
for ii = 1:4
 if ii ~= 4
  sea_trend_Ts(:,:,ii) = squeeze(mean(mon_trend_Ts(:,:,ii*3:ii*3+2,1),3));
  sea_trend_rhs(:,:,ii) = squeeze(mean(mon_trend_rhs(:,:,ii*3:ii*3+2,1),3));
 elseif ii == 4
    ii1=[1,2,12];
    sea_trend_Ts(:,:,ii) = squeeze(mean(mon_trend_Ts(:,:,ii1,1),3));
    sea_trend_rhs(:,:,ii) = squeeze(mean(mon_trend_rhs(:,:,ii1,1),3));
 end
end
% years mean trend (unite:per day)
yr_trend_Ts = squeeze(mean(mon_trend_Ts(:,:,:,1),3));
yr_trend_rhs = squeeze(mean(mon_trend_rhs(:,:,:,1),3));


%----------------------save varity----------------------------------
outfilename = '/home/lyc/repeat_Kernel_experiment/testdata/trend_tsrad';
save(outfilename, 'lon', 'lat', 'time', 'mon_trend_Ts', 'mon_trend_rhs', 'cons_Ts', 'cons_rhs', 'p_Ts', 'p_rhs','sea_trend_Ts','sea_trend_rhs',...
'yr_trend_Ts','yr_trend_rhs')

t = toc; disp(t)