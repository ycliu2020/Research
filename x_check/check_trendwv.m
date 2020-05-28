%% Combine the ERAi climatological data from 2000-03 to 2018-02
% attention: ts date used by xieyan's method, one popuse is to check the ts trend 
% cal the trend of q and ts to check the rad result right or not
inputpath1 = '/home/lyc/research/02.Ts_change_research/testdata/ERAi/ERAi_variAnom2018.nc';
inputpath2 ='/home/lyc/data/era_download/water_vapor/wv_variAnom2018.mat';% 'lon','lat','tcwv','Anom_tcwv','Clim_tcwv','time'
load(inputpath2)
time = ncread(inputpath1, 'time'); ntime = length(time);
lon = ncread(inputpath1, 'longitude'); nlon = length(lon);
lat = ncread(inputpath1, 'latitude'); nlat = length(lat);

ts_anom = ncread(inputpath1, 'Anom_ts'); % Skin temperature
wv_anom = Anom_tcwv; % water vapor content 

mon_trend_Ts = zeros(nlon, nlat, 12, 2); mon_trend_wv = zeros(nlon, nlat, 12, 2); % trend, 12 months
% im==1 is MAM, im==2 is JJA, im==3 is SON, im==4 is DJF
sea_trend_Ts = zeros(nlon, nlat, 4); sea_trend_wv = zeros(nlon, nlat, 4); % trend, 4 seasons
% yr_trend_Ts = zeros(nlon, nlat); yr_trend_wv = zeros(nlon, nlat); % trend, yearly

p_Ts = zeros(nlon, nlat, 12); p_wv = zeros(nlon, nlat, 12); % f test
cons_Ts = zeros(nlon, nlat, 12); cons_wv = zeros(nlon, nlat, 12); % mean

for i = 1:nlon

    for j = 1:nlat
        [~, mon_trend_Ts(i, j, :, :), cons_Ts(i, j, :), p_Ts(i, j, :)] = detrend_yan(ts_anom(i, j, :), time);
        [~, mon_trend_wv(i, j, :, :), cons_wv(i, j, :), p_wv(i, j, :)] = detrend_yan(wv_anom(i, j, :), time);
    end

end

% seasons mean trend (unite:per day)
for ii = 1:4

    if ii ~= 4
        sea_trend_Ts(:, :, ii) = squeeze(mean(mon_trend_Ts(:, :, ii * 3:ii * 3 + 2, 1), 3));
        sea_trend_wv(:, :, ii) = squeeze(mean(mon_trend_wv(:, :, ii * 3:ii * 3 + 2, 1), 3));
    elseif ii == 4
        ii1 = [1, 2, 12];
        sea_trend_Ts(:, :, ii) = squeeze(mean(mon_trend_Ts(:, :, ii1, 1), 3));
        sea_trend_wv(:, :, ii) = squeeze(mean(mon_trend_wv(:, :, ii1, 1), 3));
    end

end

% years mean trend (unite:per day)
yr_trend_Ts = squeeze(mean(mon_trend_Ts(:, :, :, 1), 3));
yr_trend_wv = squeeze(mean(mon_trend_wv(:, :, :, 1), 3));

%----------------------save varity----------------------------------
outfilename = '/home/lyc/research/02.Ts_change_research/testdata/trend_tswv';
save(outfilename, 'lon', 'lat', 'time', 'mon_trend_Ts', 'mon_trend_wv', 'cons_Ts', 'cons_wv', 'p_Ts', 'p_wv','sea_trend_Ts','sea_trend_wv',...
'yr_trend_Ts','yr_trend_wv')

t = toc; disp(t)