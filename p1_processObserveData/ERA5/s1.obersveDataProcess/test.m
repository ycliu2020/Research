%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-07-03 16:52:42
% LastEditTime : 2020-07-05 20:52:01
% LastEditors  : LYC
% Description  : 
% FilePath     : /code/p1_processObserveData/ERA5/s1.obersveDataProcess/test.m
%  
%%---------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% start to process
for p_1 = 1:2% 1 mean 200003-201802
    regridPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{1});
    anomPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{2});
    anomTrendPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{3});
    kernelCalPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{4});
    radEffectPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{5});
    radEffectTrendPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{6});
    auto_mkdir(regridPath);auto_mkdir(anomPath);auto_mkdir(anomTrendPath);
    auto_mkdir(kernelCalPath);auto_mkdir(radEffectPath);auto_mkdir(radEffectTrendPath);
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
    q = q0(:, :, :, startMon:endMon);
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
    startmonth = tLin.startmonth{p_1};
    [dalb, clim_alb] = monthlyAnomaly3D(nlonk, nlatk, time, alb, startmonth);
    [dts, clim_ts] = monthlyAnomaly3D(nlonk, nlatk, time, ts, startmonth);
    [dq, clim_q] = monthlyAnomaly4D(nlonk, nlatk, nplevk, time, q, startmonth);
    [dta, clim_ta] = monthlyAnomaly4D(nlonk, nlatk, nplevk, time, ta, startmonth);
    % Save
    save([anomPath, 'global_vars.mat'], 'lon_k', 'lat_k', 'time', 'plev_k')
    save([anomPath, 'meto_dvars.mat'], 'dalb', 'dts', 'dq', 'dta' ...
        , 'clim_alb', 'clim_ts', 'clim_q', 'clim_ta')

    save([regridPath, 'global_vars.mat'], 'lon_k', 'lat_k', 'time', 'plev_k')
    save([regridPath, 'meto_vars.mat'], 'alb', 'ts', 'q', 'ta')
    % rad vars
    % transform unite: W/m2
    rad_sfc = rad_sfc ./ (3600 * 24);
    rad_toa = rad_toa ./ (3600 * 24);
    % deseasonlize
    dvars.sfcRad = strcat(var_state{1}, vars.sfcRad); clim_vars.sfcRad = strcat(var_state{2}, vars.sfcRad);
    dvars.toaRad = strcat(var_state{1}, vars.toaRad); clim_vars.toaRad = strcat(var_state{2}, vars.toaRad);
    var_state = {'d', 'clim_', 'trendm_d', 'trends_d', 'trendyr_d'};

    for radVarNum = 1:length(vars.sfcRad)
        eval(['[', dvars.sfcRad{radVarNum}, ', ', clim_vars.sfcRad{radVarNum}, ']=monthlyAnomaly3D(nlonk, nlatk, time, squeeze(rad_sfc(:,:,:,radVarNum)), startmonth);']); %
        save([anomPath, dvars.sfcRad{radVarNum}], dvars.sfcRad{radVarNum}, clim_vars.sfcRad{radVarNum});
    end

    for radVarNum = 1:length(vars.toaRad)
        eval(['[', dvars.toaRad{radVarNum}, ', ', clim_vars.toaRad{radVarNum}, ']=monthlyAnomaly3D(nlonk, nlatk, time, squeeze(rad_toa(:,:,:,radVarNum)), startmonth);']); %
        save([anomPath, dvars.toaRad{radVarNum}], dvars.toaRad{radVarNum}, clim_vars.toaRad{radVarNum});
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Part3: cal specific rad componet (rhs, rhsPlus) dim positive downwards.dR_swnet=dssr;dR_lwdow=dstrd
    %% cal the deltaE(4*epis*T^3^delatT)
    epis = 1; % emissivity
    pnum = startmonth - 2;
    dDeltaE = zeros(nlonk, nlatk, ntime);
    for ii = 1:length(dts)
        dDeltaE(:, :, ii) = 4 * epis * sigma * squeeze(clim_ts(:, :, mod(ii + pnum, 12) + 1)).^3 .* dts(:, :, ii);
    end
    dDeltaE = autoRegrid3(lon_k, lat_k, time, dDeltaE, lon_f, lat_f, time);
    save([radEffectPath, 'dDeltaE.mat'], 'dDeltaE', 'lon_f', 'lat_f', 'time')
    [trendm_dDeltaE, trends_dDeltaE, trendyr_dDeltaE, ~, ~] = autoCalTrend(dDeltaE, nlonf, nlatf, time, startmonth);
    save([radEffectTrendPath, 'trend_dDeltaE.mat'], 'lon_f', 'lat_f', 'time',...
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
    [trendm_dts, trends_dts, trendyr_dts, p_dts, cons_dts] = autoCalTrend(dts, nlonf, nlatf, time, startmonth);
    [trendm_drhs, trends_drhs, trendyr_drhs, p_drhs, cons_drhs] = autoCalTrend(drhs, nlonf, nlatf, time, startmonth);
    [trendm_drhsPlus, trends_drhsPlus, trendyr_drhsPlus, p_drhsPlus, cons_drhsPlus] = autoCalTrend(drhsPlus, nlonf, nlatf, time, startmonth);
    [trendm_dhFlux, trends_dhFlux, trendyr_dhFlux, p_dhFlux, cons_dhFlux] = autoCalTrend(dhFlux, nlonf, nlatf, time, startmonth);
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
    [trendm_dnetTOA, trends_dnetTOA, trendyr_dnetTOA, p_dnetTOA, cons_dnetTOA] = autoCalTrend(dnetTOA, nlonf, nlatf, time, startmonth);
    save([anomTrendPath, 'trend_dnetTOA.mat'], 'trendm_dnetTOA', 'cons_dnetTOA', 'p_dnetTOA', 'trends_dnetTOA', 'trendyr_dnetTOA');

end

t = toc; disp(t)
