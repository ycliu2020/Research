%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-07-08 13:09:46
% LastEditTime : 2020-11-24 21:12:56
% LastEditors  : LYC
% Description  : cal china(CN) and world(UN-90-90) areamean rad and kernel cal rad 时间序列
%                注意辐射闭合检查值只适合晴空的情况
%                时间序列: 2000.03-2014-02
% FilePath     : /code/p3_paperFigIntegrate/Fig1_2_radClosure/ERA5_dRAreaMean.m
%
%%---------------------------------------------------------
clc; clear; tic;
% constant
lon_k = 0:2.5:357.5; nlonk = length(lon_k); % kernel lat lon
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);
var_state = {'d', 'clim_', 'trendm_d', 'trends_d', 'dRclr_sfcKernTemp_d'};
areaName = {'china', 'world'};
areaLabel = {'cn', 'un'};
Label.wave = {'net', 'lw', 'sw'};
Label.radEf = {'wvlw', 'wvsw', 'tslw', 'albsw', 'talw', 'res'};
[readme, level, tLin, vars] = obsParameters('ERA5');

for exmNum = 1:1 % only read 2000.03-2018.02 to cut to 2000.03-2014-02
    varsPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{1}); %rawdata
    dvarsPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{2}); %anomaly
    kernelCalPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{4}); % kernelCal
    radEfectPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{5}); %radEffect
    trend_radEfectPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{6}); %/data1/liuyincheng/cmip6-proces/aimp_2000-2014/MRI-ESM2-0/ensemble/radEffect_trend/
    %% load data
    % 观测值
    load([radEfectPath, 'global_vars.mat'])% 'lon_f', 'lat_f', 'lon_k', 'lat_k', 'plev_k', 'time'
    load([radEfectPath, 'real_dradEfect.mat'])% dR_allsky, dR_clr, l_rad, readme_realradEfect, s_rad (final row: 1sfc cld,2sfc clr,3toa cld,4toa clr)
    % cut to 2000.03-2014-02
    timeStr=string(datestr(datenum(time),'yyyy-mm-dd'));
    cutEnd=find(timeStr=='2014-02-01');
    time=time(1:cutEnd);
    ntime = length(time);
    dRclr_sfc.net = squeeze(dR_clr(:, :, 1:cutEnd, 1));
    dRclr_sfc.sw = squeeze(s_rad(:, :, 1:cutEnd, 2));
    dRclr_sfc.lw = squeeze(l_rad(:, :, 1:cutEnd, 2));
    dRclr_toa.net = squeeze(dR_clr(:, :, 1:cutEnd, 2));
    dRclr_toa.sw = squeeze(s_rad(:, :, 1:cutEnd, 4));
    dRclr_toa.lw = squeeze(l_rad(:, :, 1:cutEnd, 4));

    % toa kernel rad
    load([radEfectPath, 'dradEfect_toa_clr.mat'])%'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', 'taEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect'
    load([radEfectPath, 'dR_residual_clr_toa.mat'])%dR_residual_clr_toa
    dRclr_toaKern.lw = wvlwEffect(:,:,1:cutEnd) + tsEffect(:,:,1:cutEnd) + taEffect(:,:,1:cutEnd);
    dRclr_toaKern.sw = wvswEffect(:,:,1:cutEnd) + albEffect(:,:,1:cutEnd);
    dRclr_toaKern.net = dRclr_toaKern.sw + dRclr_toaKern.lw;

    % sfc kernel rad
    load([radEfectPath, 'dradEfect_sfc_clr.mat'])%'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', 'taEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect'
    load([radEfectPath, 'dR_residual_clr_sfc.mat'])%dR_residual_clr_toa
    dRclr_sfcKern.lw = wvlwEffect(:,:,1:cutEnd) + tsEffect(:,:,1:cutEnd) + taEffect(:,:,1:cutEnd);
    dRclr_sfcKern.sw = wvswEffect(:,:,1:cutEnd) + albEffect(:,:,1:cutEnd);
    dRclr_sfcKern.net = dRclr_sfcKern.sw + dRclr_sfcKern.lw;

    % use one var to calculation
    dRclr_sfcKernTemp = zeros(nlonf, nlatf, ntime, 3);
    dRclr_toaKernTemp = dRclr_sfcKernTemp;
    dRclr_sfcTemp = zeros(nlonf, nlatf, ntime, 3);
    dRclr_toaTemp = dRclr_sfcTemp;

    for waveNum = 1:length(Label.wave)
        eval(['dRclr_sfcTemp(:,:,:,waveNum)=dRclr_sfc.', Label.wave{waveNum}, ';'])
        eval(['dRclr_toaTemp(:,:,:,waveNum)=dRclr_toa.', Label.wave{waveNum}, ';'])
        eval(['dRclr_sfcKernTemp(:,:,:,waveNum)=dRclr_sfcKern.', Label.wave{waveNum}, ';'])
        eval(['dRclr_toaKernTemp(:,:,:,waveNum)=dRclr_toaKern.', Label.wave{waveNum}, ';'])
    end

    dRclr_sfcRes = dRclr_sfcTemp - dRclr_sfcKernTemp; % (144,72,time, 3)
    dRclr_toaRes = dRclr_toaTemp - dRclr_toaKernTemp; % (144,72,time, 3)

    % mask set
    [Lonf, Latf] = ndgrid(lon_f, lat_f);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part1: world mean rad
    maxlat = 90; areaNum = 1;
    dRclr_sfc_un = zeros(nlonf, nlatf, ntime, 3);
    dRclr_toa_un = dRclr_sfc_un; dRclr_sfcKern_un = dRclr_sfc_un; dRclr_toaKern_un = dRclr_sfc_un; dRclr_sfcRes_un = dRclr_sfc_un; dRclr_toaRes_un = dRclr_sfc_un;

    for waveNum = 1:length(Label.wave)
        dRclr_sfc_un(:, :, :, waveNum) = maskLand(squeeze(dRclr_sfcTemp(:, :, :, waveNum)), lat_f, maxlat, -maxlat, areaNum); % world
        dRclr_toa_un(:, :, :, waveNum) = maskLand(squeeze(dRclr_toaTemp(:, :, :, waveNum)), lat_f, maxlat, -maxlat, areaNum);
        dRclr_sfcKern_un(:, :, :, waveNum) = maskLand(squeeze(dRclr_sfcKernTemp(:, :, :, waveNum)), lat_f, maxlat, -maxlat, areaNum);
        dRclr_toaKern_un(:, :, :, waveNum) = maskLand(squeeze(dRclr_toaKernTemp(:, :, :, waveNum)), lat_f, maxlat, -maxlat, areaNum);
        dRclr_sfcRes_un(:, :, :, waveNum) = maskLand(squeeze(dRclr_sfcRes(:, :, :, waveNum)), lat_f, maxlat, -maxlat, areaNum);
        dRclr_toaRes_un(:, :, :, waveNum) = maskLand(squeeze(dRclr_toaRes(:, :, :, waveNum)), lat_f, maxlat, -maxlat, areaNum);
    end

    % 加权平均
    dRclrMean_sfc_un = zeros(ntime, 3);
    dRclrMean_toa_un = dRclrMean_sfc_un; dRclrMean_sfcKern_un = dRclrMean_sfc_un; dRclrMean_toaKern_un = dRclrMean_sfc_un; dRclrMean_sfcRes_un = dRclrMean_sfc_un; dRclrMean_toaRes_un = dRclrMean_sfc_un;
    wgt = cos(Latf ./ 180 .* pi); %纬度加权算区域平均值
    wgt=maskLand(wgt,lat_f, maxlat, -maxlat, areaNum);
    for waveNum = 1:length(Label.wave)

        for waveTime = 1:ntime
            dRclrMean_sfc_un(waveTime, waveNum) = nansum(nansum(squeeze(dRclr_sfc_un(:, :, waveTime, waveNum)) .* wgt)) / nansum(nansum(wgt));
            dRclrMean_toa_un(waveTime, waveNum) = nansum(nansum(squeeze(dRclr_toa_un(:, :, waveTime, waveNum)) .* wgt)) / nansum(nansum(wgt));
            dRclrMean_sfcKern_un(waveTime, waveNum) = nansum(nansum(squeeze(dRclr_sfcKern_un(:, :, waveTime, waveNum)) .* wgt)) / nansum(nansum(wgt));
            dRclrMean_toaKern_un(waveTime, waveNum) = nansum(nansum(squeeze(dRclr_toaKern_un(:, :, waveTime, waveNum)) .* wgt)) / nansum(nansum(wgt));
            dRclrMean_sfcRes_un(waveTime, waveNum) = nansum(nansum(squeeze(dRclr_sfcRes_un(:, :, waveTime, waveNum)) .* wgt)) / nansum(nansum(wgt));
            dRclrMean_toaRes_un(waveTime, waveNum) = nansum(nansum(squeeze(dRclr_toaRes_un(:, :, waveTime, waveNum)) .* wgt)) / nansum(nansum(wgt));
        end

    end

    radEfectPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{3}, 'ERA5', level.standVarPath{5}); %radEffect
    auto_mkdir(radEfectPath);
    unFileName = fullfile(radEfectPath, 'RadMean_world.mat');
    readme='(time,  band({net, lw, sw}))';
    save(unFileName, 'lon_f', 'lat_f', 'time','readme', 'dRclrMean_sfc_un', 'dRclrMean_toa_un', ...
        'dRclrMean_sfcKern_un', 'dRclrMean_toaKern_un', ...
        'dRclrMean_sfcRes_un', 'dRclrMean_toaRes_un')
    clear dRclrMean_sfc_un dRclrMean_toa_un dRclrMean_sfcKern_un dRclrMean_toaKern_un dRclrMean_sfcRes_un dRclrMean_toaRes_un
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part2: china mean rad
    maxlat = 60; areaNum = 2;
    dRclr_sfc_cn = zeros(nlonf, nlatf, ntime, 3);
    dRclr_toa_cn = dRclr_sfc_cn; dRclr_sfcKern_cn = dRclr_sfc_cn; dRclr_toaKern_cn = dRclr_sfc_cn; dRclr_sfcRes_cn = dRclr_sfc_cn; dRclr_toaRes_cn = dRclr_sfc_cn;
    wgt = cos(Latf ./ 180 .* pi); %纬度加权算区域平均值
    wgt=maskLand(wgt,lat_f, maxlat, -maxlat, areaNum);
    wgt = cos(Latf ./ 180 .* pi); %纬度加权算区域平均值
    wgt=maskLand(wgt,lat_f, maxlat, -maxlat, areaNum);
    for waveNum = 1:length(Label.wave)
        dRclr_sfc_cn(:, :, :, waveNum) = maskLand(squeeze(dRclr_sfcTemp(:, :, :, waveNum)), lat_f, maxlat, -maxlat, areaNum); % world
        dRclr_toa_cn(:, :, :, waveNum) = maskLand(squeeze(dRclr_toaTemp(:, :, :, waveNum)), lat_f, maxlat, -maxlat, areaNum);
        dRclr_sfcKern_cn(:, :, :, waveNum) = maskLand(squeeze(dRclr_sfcKernTemp(:, :, :, waveNum)), lat_f, maxlat, -maxlat, areaNum);
        dRclr_toaKern_cn(:, :, :, waveNum) = maskLand(squeeze(dRclr_toaKernTemp(:, :, :, waveNum)), lat_f, maxlat, -maxlat, areaNum);
        dRclr_sfcRes_cn(:, :, :, waveNum) = maskLand(squeeze(dRclr_sfcRes(:, :, :, waveNum)), lat_f, maxlat, -maxlat, areaNum);
        dRclr_toaRes_cn(:, :, :, waveNum) = maskLand(squeeze(dRclr_toaRes(:, :, :, waveNum)), lat_f, maxlat, -maxlat, areaNum);
    end

    % 加权平均
    dRclrMean_sfc_cn = zeros(ntime, 3);
    dRclrMean_toa_cn = dRclrMean_sfc_cn; dRclrMean_sfcKern_cn = dRclrMean_sfc_cn; dRclrMean_toaKern_cn = dRclrMean_sfc_cn; dRclrMean_sfcRes_cn = dRclrMean_sfc_cn; dRclrMean_toaRes_cn = dRclrMean_sfc_cn;

    for waveNum = 1:length(Label.wave)

        for waveTime = 1:ntime
            dRclrMean_sfc_cn(waveTime, waveNum) = nansum(nansum(squeeze(dRclr_sfc_cn(:, :, waveTime, waveNum)) .* wgt)) / nansum(nansum(wgt));
            dRclrMean_toa_cn(waveTime, waveNum) = nansum(nansum(squeeze(dRclr_toa_cn(:, :, waveTime, waveNum)) .* wgt)) / nansum(nansum(wgt));
            dRclrMean_sfcKern_cn(waveTime, waveNum) = nansum(nansum(squeeze(dRclr_sfcKern_cn(:, :, waveTime, waveNum)) .* wgt)) / nansum(nansum(wgt));
            dRclrMean_toaKern_cn(waveTime, waveNum) = nansum(nansum(squeeze(dRclr_toaKern_cn(:, :, waveTime, waveNum)) .* wgt)) / nansum(nansum(wgt));
            dRclrMean_sfcRes_cn(waveTime, waveNum) = nansum(nansum(squeeze(dRclr_sfcRes_cn(:, :, waveTime, waveNum)) .* wgt)) / nansum(nansum(wgt));
            dRclrMean_toaRes_cn(waveTime, waveNum) = nansum(nansum(squeeze(dRclr_toaRes_cn(:, :, waveTime, waveNum)) .* wgt)) / nansum(nansum(wgt));
        end

    end

    cnFileName = fullfile(radEfectPath, 'RadMean_china.mat');
    readme='(time,  band({net, lw, sw}))';
    save(cnFileName, 'lon_f', 'lat_f', 'time', 'readme', 'dRclrMean_sfc_cn', 'dRclrMean_toa_cn', ...
        'dRclrMean_sfcKern_cn', 'dRclrMean_toaKern_cn', ...
        'dRclrMean_sfcRes_cn', 'dRclrMean_toaRes_cn')
    clear dRclrMean_sfc_cn dRclrMean_toa_cn dRclrMean_sfcKern_cn dRclrMean_toaKern_cn dRclrMean_sfcRes_cn dRclrMean_toaRes_cn
    
    clear dRclr_sfc dRclr_toa dRclr_toaKern dRclr_sfcKern

end
t = toc; disp(t)
