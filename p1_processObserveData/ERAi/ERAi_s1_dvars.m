%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-07-07 12:51:42
% LastEditTime : 2020-07-08 10:03:21
% LastEditors  : LYC
% Description  : 
% FilePath     : /code/p1_processObserveData/ERAi/ERAi_s1_dvars.m
%  
%%---------------------------------------------------------
%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-07-02 15:23:02
% LastEditTime : 2020-07-08 09:27:29
% LastEditors  : LYC
% Description  : process ERAi data to anomaly (includ meto vars and rad)
%                time line: 1. 2000-03 to 2018-02(18*12) 2. 200207-201706(15*12)
%                note that all vertical fluxes is positive downwards.
% FilePath     : /code/p1_processObserveData/ERAi/ERAi_s1_dvars.m
%
%%---------------------------------------------------------
clc; clear; tic;
% constant
sigma = 5.67e-8; % Stefan-Boltzmann constant: Wm-2K-4
kernels_path = '/data1/liuyincheng/y_kernels/kernels_YiH/toa/dp.nc';
lon_k = 0:2.5:357.5; nlonk = length(lon_k); % kernel lat lon
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);
var_state = {'d', 'clim_', 'trendm_d', 'trends_d', 'trendyr_d'};
% modify path first
obsPath = '/data1/liuyincheng/Observe-rawdata/ERA/ERAi';
metoVarsPath = fullfile(obsPath, 'meto_vars/');
sfc_radVarsPath = fullfile(obsPath, 'rad_vars/SFC/');
toa_radVarsPath = fullfile(obsPath, 'rad_vars/TOA/');
[readme, level, tLin, vars] = obsParameters('ERAi');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read raw data from 2000-2018 for futher process to 200003-201802,,,,
% 1.1 连续时间的变量
alb_Path = fullfile(metoVarsPath, 'alb_1979-2018.nc');
ps_Path = fullfile(metoVarsPath, 'ps_1979-2018.nc');
ts_Path = fullfile(metoVarsPath, 'ts_1979-2018.nc');
% read time
% locate first year
formatOut = 'yyyy-mm';
temp_dir = dir(alb_Path);
startT = cdftime2loc(temp_dir, formatOut, '2000-01');
% read time and transfor mat time
startLoc = startT; count = 19 * 12; stride = 1;
time0 = cmipTimeRead(alb_Path, startLoc, count, stride); %time.date, Units, Calendar, length
ntime0 = length(time0.date);
time1 = 1:ntime0;
% read original lat and lon
lat_ori = double(ncread(alb_Path, 'latitude'));
lon_ori = double(ncread(alb_Path, 'longitude'));
% read vars
startLoc = [1, 1, startT]; count = [inf, inf, 19 * 12]; stride = [1, 1, 1];
alb0 = squeeze(ncread(alb_Path, 'fal', startLoc, count, stride));
ps0 = squeeze(ncread(ps_Path, 'sp', startLoc, count, stride));
ts0 = squeeze(ncread(ts_Path, 'skt', startLoc, count, stride));
% Regrid to 2.5*2.5
% 3D meto vars
ps0 = autoRegrid3(lon_ori, lat_ori, time1, ps0, lon_k, lat_k, time1);
alb0 = autoRegrid3(lon_ori, lat_ori, time1, alb0, lon_k, lat_k, time1);
ts0 = autoRegrid3(lon_ori, lat_ori, time1, ts0, lon_k, lat_k, time1);

% 1.2 按时间划分文件的变量
for yearNum = 0:18% (read raw range: 2000-2018)
    % meto data
    TQ_Path = strcat(metoVarsPath, 'tq', num2str((yearNum + 2000), '%04i'), '.nc');
    ta0(:, :, :, yearNum * 12 + 1:(yearNum + 1) * 12) = ncread(TQ_Path, 't');
    hus0(:, :, :, yearNum * 12 + 1:(yearNum + 1) * 12) = ncread(TQ_Path, 'q');

    % rad data({'tsr', 'ttr', 'tsrcs', 'ttrcs', 'slhf', 'sshf', 'ssr', 'str', 'strd', 'ssrd', 'ssrc', 'strc', 'strdc', 'ssrdc'};)
    sfcRad_Path = strcat(sfc_radVarsPath, 'rad_sfc', num2str((yearNum + 2000), '%04i'), '.nc');
    toaRad_Path = strcat(toa_radVarsPath, 'rad_toa', num2str((yearNum + 2000), '%04i'), '.nc');

    for radVarNum = 1:length(vars.sfcRad)
        rad_sfc_ori(:, :, yearNum * 12 + 1:(yearNum + 1) * 12, radVarNum) = ncread(sfcRad_Path, vars.sfcRad{radVarNum});
    end

    for radVarNum = 1:length(vars.toaRad)
        rad_toa_ori(:, :, yearNum * 12 + 1:(yearNum + 1) * 12, radVarNum) = ncread(toaRad_Path, vars.toaRad{radVarNum});
    end

end

% read level
plev_ori = ncread(TQ_Path, 'level');
plev_ori = double(plev_ori);
plev_ori=plev_ori(end:-1:1);
plev_k = ncread('/data1/liuyincheng/y_kernels/kernels_YiH/toa/dp.nc', 'player'); % player 最大1000符合条件pay attention to plev_k's range must smaller than plev_ori's
nplevk = length(plev_k);

% transform plev direction
hus0=hus0(:,:,end:-1:1,:);
ta0=ta0(:,:,end:-1:1,:);

% Regrid to 2.5*2.5
% 4D meto vars
hus0 = autoRegrid4(lon_ori, lat_ori, plev_ori, time1, hus0, lon_k, lat_k, plev_k, time1);
ta0 = autoRegrid4(lon_ori, lat_ori, plev_ori, time1, ta0, lon_k, lat_k, plev_k, time1);

% rad vars
rad_sfc0 = zeros(nlonk, nlatk, ntime0, length(vars.sfcRad));
rad_toa0 = zeros(nlonk, nlatk, ntime0, length(vars.toaRad));

for radVarNum = 1:length(vars.toaRad)
    rad_toa0(:, :, :, radVarNum) = autoRegrid3(lon_ori, lat_ori, time1, rad_toa_ori(:, :, :, radVarNum), lon_k, lat_k, time1);
end

clear rad_toa_ori

for radVarNum = 1:length(vars.sfcRad)
    rad_sfc0(:, :, :, radVarNum) = autoRegrid3(lon_ori, lat_ori, time1, rad_sfc_ori(:, :, :, radVarNum), lon_k, lat_k, time1);
end

clear rad_sfc_ori

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% start to process
for p_1 = 1:2% 1 mean 200003-201802
    regridPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi', level.standVarPath{1});
    anomPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi', level.standVarPath{2});
    anomTrendPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi', level.standVarPath{3});
    kernelCalPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi', level.standVarPath{4});
    radEffectPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi', level.standVarPath{5});
    radEffectTrendPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi', level.standVarPath{6});
    auto_mkdir(regridPath); auto_mkdir(anomPath); auto_mkdir(anomTrendPath);
    auto_mkdir(kernelCalPath); auto_mkdir(radEffectPath); auto_mkdir(radEffectTrendPath);
    % find start and end month
    timeSpec = tLin.start{p_1};
    formatOut = 'yyyy-mm';
    timeStr = string(datestr(time0.date, formatOut));
    startMon = find(timeStr == timeSpec);
    endMon = startMon + tLin.inter{p_1} - 1;

    % read specific time series
    ps = ps0(:, :, startMon:endMon);
    alb = alb0(:, :, startMon:endMon);
    ts = ts0(:, :, startMon:endMon);
    hus = hus0(:, :, :, startMon:endMon);
    ta = ta0(:, :, :, startMon:endMon);
    time = time0.date(startMon:endMon, 1);
    ntime = length(time);
    rad_sfc = rad_sfc0(:, :, startMon:endMon, :);
    rad_toa = rad_toa0(:, :, startMon:endMon, :);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Part1: handle ps
    % cal average ps
    ps = ps ./ 100; % transform hp unit
    ps_m = zeros(nlonk, nlatk, 12);
    [yy, mm, dd] = datevec(time);

    for im = 1:12
        index_mm = mm == im;
        ps_m(:, :, im) = nanmean(ps(:, :, index_mm), 3); %this is the average of each month
    end

    month = 1:12;
    % cal dp and dps(bottom level height)
    plevel = ncread(kernels_path, 'plevel');
    player = ncread(kernels_path, 'player');
    dp_raw = ncread(kernels_path, 'dp');
    dp = zeros(nlonk, nlatk, 24, 12); %revise dp in all layers, under the lowest layer are zeros
    dp_level2 = zeros(nlonk, nlatk, 24, 12); % only contain the near sfc
    dps = zeros(nlonk, nlatk, 12);

    for i = 1:nlonk

        for j = 1:nlatk

            for nt = 1:12
                temp = find(player < ps_m(i, j, nt), 1, 'first');
                dps(i, j, nt) = ps_m(i, j, nt) - plevel(temp + 1);
                % dps(i, j, nt) = dp_bottom(i,j)123;
                dp(i, j, temp, nt) = dps(i, j, nt);
                dp(i, j, temp + 1:24, nt) = dp_raw(temp + 1:24);
                % near surface level2
                dp_level2(i, j, temp, nt) = dps(i, j, nt);
                dp_level2(i, j, temp + 1, nt) = dp_raw(temp + 1);
            end

        end

    end

    % Save
    readme_ps_m.unite = 'hpa';
    readme_ps_m.state = 'contain 12 months';
    readme_ps_m.longName = 'Surface Pressure(monthly mean of mutiple years)';
    save([kernelCalPath, 'global_vars.mat'], 'lon_k', 'lat_k', 'time', 'plev_k')
    save([kernelCalPath, 'kernel_ps_m.mat'], 'ps_m', 'readme_ps_m')
    save([kernelCalPath, 'kernel_dp.mat'], 'dps', 'dp', 'dp_level2');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Part2: deseasonlize vars
    % meto vars
    startMonth = tLin.startMonth{p_1};
    [dalb, clim_alb] = monthlyAnomaly3D(nlonk, nlatk, time, alb, startMonth);
    [dts, clim_ts] = monthlyAnomaly3D(nlonk, nlatk, time, ts, startMonth);
    [dhus, clim_hus] = monthlyAnomaly4D(nlonk, nlatk, nplevk, time, hus, startMonth);
    [dta, clim_ta] = monthlyAnomaly4D(nlonk, nlatk, nplevk, time, ta, startMonth);
    % Save
    save([anomPath, 'global_vars.mat'], 'lon_k', 'lat_k', 'time', 'plev_k')
    save([anomPath, 'meto_dvars.mat'], 'dalb', 'dts', 'dhus', 'dta' ...
        , 'clim_alb', 'clim_ts', 'clim_hus', 'clim_ta')

    save([regridPath, 'global_vars.mat'], 'lon_k', 'lat_k', 'time', 'plev_k')
    save([regridPath, 'meto_vars.mat'], 'alb', 'ts', 'hus', 'ta')
    % rad vars
    % transform unite: W/m2
    rad_sfc = rad_sfc ./ (3600 * 24);
    rad_toa = rad_toa ./ (3600 * 24);
    % deseasonlize
    dvars.sfcRad = strcat(var_state{1}, vars.sfcRad); clim_vars.sfcRad = strcat(var_state{2}, vars.sfcRad);
    dvars.toaRad = strcat(var_state{1}, vars.toaRad); clim_vars.toaRad = strcat(var_state{2}, vars.toaRad);

    for radVarNum = 1:length(vars.sfcRad)
        eval([vars.sfcRad{radVarNum},'=squeeze(rad_sfc(:,:,:,radVarNum));']); %
        save([regridPath, vars.sfcRad{radVarNum}], vars.sfcRad{radVarNum});

        eval(['[', dvars.sfcRad{radVarNum}, ', ', clim_vars.sfcRad{radVarNum}, ']=monthlyAnomaly3D(nlonk, nlatk, time, squeeze(rad_sfc(:,:,:,radVarNum)), startMonth);']); %
        save([anomPath, dvars.sfcRad{radVarNum}], dvars.sfcRad{radVarNum}, clim_vars.sfcRad{radVarNum});
        
    end

    for radVarNum = 1:length(vars.toaRad)
        eval([vars.toaRad{radVarNum},'=squeeze(rad_toa(:,:,:,radVarNum));']); %
        save([regridPath, vars.toaRad{radVarNum}], vars.toaRad{radVarNum});

        eval(['[', dvars.toaRad{radVarNum}, ', ', clim_vars.toaRad{radVarNum}, ']=monthlyAnomaly3D(nlonk, nlatk, time, squeeze(rad_toa(:,:,:,radVarNum)), startMonth);']); %
        save([anomPath, dvars.toaRad{radVarNum}], dvars.toaRad{radVarNum}, clim_vars.toaRad{radVarNum});
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Part3: cal specific rad componet (rhs, rhsPlus) dim positive downwards.dR_swnet=dssr;dR_lwdow=dstrd
    %% cal the deltaE(4*epis*T^3^delatT)
    epis = 1; % emissivity
    pnum = startMonth - 2;
    dDeltaE = zeros(nlonk, nlatk, ntime);

    for ii = 1:length(dts)
        dDeltaE(:, :, ii) = 4 * epis * sigma * squeeze(clim_ts(:, :, mod(ii + pnum, 12) + 1)).^3 .* dts(:, :, ii);
    end

    dDeltaE = autoRegrid3(lon_k, lat_k, time, dDeltaE, lon_f, lat_f, time);
    save([radEffectPath, 'dDeltaE.mat'], 'dDeltaE', 'lon_f', 'lat_f', 'time')
    [trendm_dDeltaE, trends_dDeltaE, trendyr_dDeltaE, ~, ~] = autoCalTrend(dDeltaE, nlonf, nlatf, time, startMonth);
    save([radEffectTrendPath, 'trend_dDeltaE.mat'], 'lon_f', 'lat_f', 'time', ...
        'trendm_dDeltaE', 'trends_dDeltaE', 'trendyr_dDeltaE')

    %% sfc rad Flux
    dhFlux = dslhf + dsshf;
    drhs = dstrd + dssr;
    drhsPlus = drhs + dhFlux;
    save([anomPath, 'drhs.mat'], 'drhs', 'drhsPlus', 'dhFlux');
    %regrid 144x72(unite grids)
    dhFlux = autoRegrid3(lon_k, lat_k, time, dhFlux, lon_f, lat_f, time);
    dts = autoRegrid3(lon_k, lat_k, time, dts, lon_f, lat_f, time);
    drhs = autoRegrid3(lon_k, lat_k, time, drhs, lon_f, lat_f, time);
    drhsPlus = autoRegrid3(lon_k, lat_k, time, drhsPlus, lon_f, lat_f, time);
    % cal the trend
    [trendm_dts, trends_dts, trendyr_dts, p_dts, cons_dts] = autoCalTrend(dts, nlonf, nlatf, time, startMonth);
    [trendm_drhs, trends_drhs, trendyr_drhs, p_drhs, cons_drhs] = autoCalTrend(drhs, nlonf, nlatf, time, startMonth);
    [trendm_drhsPlus, trends_drhsPlus, trendyr_drhsPlus, p_drhsPlus, cons_drhsPlus] = autoCalTrend(drhsPlus, nlonf, nlatf, time, startMonth);
    [trendm_dhFlux, trends_dhFlux, trendyr_dhFlux, p_dhFlux, cons_dhFlux] = autoCalTrend(dhFlux, nlonf, nlatf, time, startMonth);
    % save the vars
    save([anomTrendPath, 'trend_dts.mat'], 'trendm_dts', 'cons_dts', 'p_dts', 'trends_dts', 'trendyr_dts');
    save([anomTrendPath, 'trend_drhs.mat'], 'trendm_drhs', 'cons_drhs', 'p_drhs', 'trends_drhs', 'trendyr_drhs');
    save([anomTrendPath, 'trend_drhsPlus.mat'], 'trendm_drhsPlus', 'cons_drhsPlus', 'p_drhsPlus', 'trends_drhsPlus', 'trendyr_drhsPlus');
    save([anomTrendPath, 'trend_dhFlux.mat'], 'trendm_dhFlux', 'cons_dhFlux', 'p_dhFlux', 'trends_dhFlux', 'trendyr_dhFlux');
    save([anomTrendPath, 'global_vars.mat'], 'lon_f', 'lat_f', 'time', 'plev_k')
    %% toa rad Flux
    dnetTOA = dtsr + dttr;
    save([anomPath, 'dnetTOA.mat'], 'dnetTOA');
    % cal the trend
    dnetTOA = autoRegrid3(lon_k, lat_k, time, dnetTOA, lon_f, lat_f, time);
    [trendm_dnetTOA, trends_dnetTOA, trendyr_dnetTOA, p_dnetTOA, cons_dnetTOA] = autoCalTrend(dnetTOA, nlonf, nlatf, time, startMonth);
    save([anomTrendPath, 'trend_dnetTOA.mat'], 'trendm_dnetTOA', 'cons_dnetTOA', 'p_dnetTOA', 'trends_dnetTOA', 'trendyr_dnetTOA');

end

t = toc; disp(t)
