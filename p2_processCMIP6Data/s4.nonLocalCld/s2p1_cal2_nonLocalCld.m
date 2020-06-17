%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-17 13:16:32
% LastEditTime : 2020-06-17 15:02:09
% LastEditors  : LYC
% Description  : 计算非局地云的辐射效应和温度贡献(用云辐射贡献的ts变化求)
% FilePath     : /code/p2_processCMIP6Data/s4.nonLocalCld/s2p1_cal2_nonLocalCld.m
% Attention    : both amip and ssp on surface caled
%%---------------------------------------------------------

clear; clc; tic;
p_3 = 88.75; % Latitude range
lon1 = [2.5 357.5]; lat1 = [-p_3 + 1 p_3 - 1]; % world area
toaSfc = {'toa', 'sfc'};

for p_1 = 1:4%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    nowpath = pwd;
    [readme, Experiment, level, tLin, mPlev, vars] = modelParameters(p_1);
    % inputPath
    inputPath = ['/data1/liuyincheng/cmip6-process/', level.time1{p_1}]; %/data1/liuyincheng/cmip6-process/amip_2000-2014/
    dvarsPath = [inputPath, level.model2{1}, '/', level.process3{6}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
    load([dvarsPath, 'global_vars.mat'])% lat lon time plevf readme
    nlat = length(lat); nlon = length(lon); ntime = length(time);

    for level1 = 1:length(level.model2)
        % data path
        varsPath = [inputPath, level.model2{level1}, '/', level.process3{1}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/rawdata
        dvarsPath = [inputPath, level.model2{level1}, '/', level.process3{2}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
        dvarsTrendPath = [inputPath, level.model2{level1}, '/', level.process3{3}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly_trend
        kernelPath = [inputPath, level.model2{level1}, '/', level.process3{5}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/kernelsCal
        dradEffectPath = [inputPath, level.model2{level1}, '/', level.process3{6}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
        dnonLocalCldPath = [inputPath, level.model2{level1}, '/', level.process3{8}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/non_localCld/
        % load dtsg and DeltaTsg
        load([dvarsPath, 'dtsg_nomask.mat'])
        load([dvarsTrendPath, 'DeltaTsg_nomask.mat'])% DeltaTsg, DeltaTsgMeaning

        varUsed = zeros(nlon, nlat, ntime, 3); % dR
        varKerUsed = varUsed; % dR/kernels
        %% Step1: calculate Heating(J/m2) and DeltaTsg_cld
        month_second = 30 * 24 * 60 * 60; % sum seconds in 1 month
        % cal total heating
        load([varsPath, 'global_vars.mat'])
        load([varsPath, 'rlut.mat'])% toa_outgoing_longwave_flux
        load([varsPath, 'rsut.mat'])% toa_outgoing_shortwave_flux
        load([varsPath, 'rsdt.mat'])% toa_incoming_shortwave_flux
        % unite define a vector component which is positive when directed downward
        netTOA = rsdt - rlut - rsut; % net radiation in toa (lonxlatxmonth): W/m2
        %regrid to 144x72(unite grids)
        lat = 88.75:-2.5:-88.75; nlat = length(lat);
        lon = lonf; nlon = length(lon);
        nlonf = length(lonf); nlatf = length(latf);
        [Xlon, Ylat, Ttime] = meshgrid(lat, lon, time);
        [Xlonf, Ylatf, Ttimef] = meshgrid(latf, lonf, time);
        netTOA = interp3(Xlonf, Ylatf, Ttimef, netTOA, Xlon, Ylat, Ttime);

        Heat_tot = sum(netTOA, 3) * month_second;

        % cal cloud heating
        load([dradEffectPath, 'dR_cloud_toa.mat'])% dR_cloud_toa
        % control state = The average state in the first five years
        ctrl_Rcloud = mean(dR_cloud_toa(:, :, 1:5 * 12), 3);
        ctrl_Rcloud = repmat(ctrl_Rcloud, [1 1 length(time)]);
        delta_Rcloud = dR_cloud_toa - ctrl_Rcloud;
        Heat_cld = sum(delta_Rcloud, 3) * month_second;

        % global Zonal weighted average (time, vars)
        jiaquan = cosd(lat);
        wei = ones(144, 72); %格点纬度加权

        for latiNum = 1:72
            wei(:, latiNum) = wei(:, latiNum) * jiaquan(latiNum); %格点相对大小
        end

        Heat_totAverage = nansum(nansum(Heat_tot .* wei)) / nansum(nansum(wei));
        Heat_cldAverage = nansum(nansum(Heat_cld .* wei)) / nansum(nansum(wei));
        % cal DeltaTsg_cld
        DeltaTsg_cld = Heat_cldAverage / Heat_totAverage * DeltaTsg;
        % save
        save([dradEffectPath, 'DeltaTsg_cld.mat'], 'DeltaTsg_cld')

        %% Step2: cal dRnonLocalCld2 and trendyr_dRnonLocalCld2(unite: per day)
        load([dradEffectPath, 'dradEfect_sfc_cld.mat'])%albEffect, husEffect, taEffect, mainEffect, totalEffect, tsEffect, talwEffect, taswEffect
        varUsed(:, :, :, 1) = husEffect;
        varUsed(:, :, :, 2) = taEffect;
        varUsed(:, :, :, 3) = albEffect;
        varParitial = zeros(nlon, nlat, 3);
        varKerParitial = varParitial;
        % Part1: cal trendyr_dRnonLocalCld_toa/sfc(unite: per day)
        for varNum = 1:3
            tempTemp = squeeze(varUsed(:, :, :, varNum));

            for latNum = 1:nlat

                for lonNum = 1:nlon
                    tempUsed = squeeze(squeeze(tempTemp(lonNum, latNum, :)));
                    k = polyfit(dtsg, tempUsed, 1); % 一元一次拟合
                    varParitial(lonNum, latNum, varNum) = k(1);
                end

            end

        end

        dRnonLocalCld2_hus = squeeze(varParitial(:, :, 1)) * DeltaTsg_cld;
        dRnonLocalCld2_ta = squeeze(varParitial(:, :, 2)) * DeltaTsg_cld;
        dRnonLocalCld2_alb = squeeze(varParitial(:, :, 3)) * DeltaTsg_cld;
        dRnonLocalCld2 = dRnonLocalCld2_hus + dRnonLocalCld2_ta + dRnonLocalCld2_alb;

        trendyr_dRnonLocalCld2 = dRnonLocalCld2 ./ (time(end) - time(1));
        save([dnonLocalCldPath, 'dRnonLocalCld2_sfc.mat'], 'dRnonLocalCld2_hus', 'dRnonLocalCld2_ta', 'dRnonLocalCld2_alb', 'dRnonLocalCld2')
        save([dnonLocalCldPath, 'trendyr_dRnonLocalCld2_sfc.mat'], 'trendyr_dRnonLocalCld2')

        % Part2: cal trendyr_dTsnonLocalCld2_toa/sfc(unite: per day)
        % load ts kernel
        load([kernelPath, '/kernels_sfc_cld'], 'ts_lwkernel');
        load([kernelPath, '/global_vars.mat']); % latf,lonf,time
        % regrid ts_lwkernel to 144x72(unite grids)
        kernelTime = 1:12;
        [Xlon, Ylat, Ttime] = meshgrid(lat, lon, kernelTime);
        [Xlonf, Ylatf, Ttimef] = meshgrid(latf, lonf, kernelTime);
        ts_lwkernel = interp3(Xlonf, Ylatf, Ttimef, ts_lwkernel, Xlon, Ylat, Ttime);
        % extend to the whole time series
        startmonth = 1;

        if startmonth ~= 1
            tempK = zeros(nlon, nlat, 12);
            tempK(1:12 - startmonth + 1) = ts_lwkernel(:, :, startmonth:12);
            tempK(12 - startmonth + 2:12) = ts_lwkernel(:, :, 1:startmonth - 1);
            ts_lwkernel = tempK;
        end

        nyear = ntime / 12;
        ts_lwkernelSeries = repmat(ts_lwkernel, [1 1 nyear]);

        for varNum = 1:3
            varKerUsed(:, :, :, varNum) = squeeze(varUsed(:, :, :, varNum)) ./ -ts_lwkernelSeries;
        end

        for varNum = 1:3
            tempTemp = squeeze(varKerUsed(:, :, :, varNum));

            for latNum = 1:nlat

                for lonNum = 1:nlon
                    tempUsed = squeeze(squeeze(tempTemp(lonNum, latNum, :)));
                    k = polyfit(dtsg, tempUsed, 1); % 一元一次拟合
                    varKerParitial(lonNum, latNum, varNum) = k(1);
                end

            end

        end

        dTsnonLocalCld2_hus = squeeze(varKerParitial(:, :, 1)) * DeltaTsg_cld;
        dTsnonLocalCld2_ta = squeeze(varKerParitial(:, :, 2)) * DeltaTsg_cld;
        dTsnonLocalCld2_alb = squeeze(varKerParitial(:, :, 3)) * DeltaTsg_cld;

        dTsnonLocalCld2 = dTsnonLocalCld2_hus + dTsnonLocalCld2_ta + dTsnonLocalCld2_alb;
        trendyr_dTsnonLocalCld2 = dTsnonLocalCld2 ./ (time(end) - time(1));

        save([dnonLocalCldPath, 'dTsnonLocalCld2_sfc.mat'], 'dTsnonLocalCld2_hus', 'dTsnonLocalCld2_ta', 'dTsnonLocalCld2_alb', 'dTsnonLocalCld2')
        save([dnonLocalCldPath, 'trendyr_dTsnonLocalCld2_sfc.mat'], 'trendyr_dTsnonLocalCld2')

        disp([level.model2{level1}, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{p_1}, ' era is done!'])
    disp(' ')

end

t = toc; disp(t)
