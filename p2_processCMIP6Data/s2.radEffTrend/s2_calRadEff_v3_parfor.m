%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-09 15:52:00
% LastEditTime : 2021-03-16 21:19:54
% LastEditors  : Please set LastEditors
% Description  : cal mainly include 1.regrid vars, 2.vars anomly
%                CMIP6 mothly data
%                time:2000.01-2014.12(interval:15*12);1980.01-2014.12(interval:35*12); 2015.01-2099.12(interval:85*12)
%                initial time in amip(432 total): 253 of 432(2000.03);13 of 432(1980.01);
%                initial time in futrue(1032 total): 1 of 1032(2015.01);
%                initial time in amip-hist(1740 total): 1,561 of 1740(2000.03);1,321 of 1740(1980.01);
%                PS: m mean month, s mean season, yr mean year
% inputPath    : /code/p2_processCMIP6Data/s1.modelDataProcess/
% Attention!!!
% check lat_f: model lat_f disagree with kernels lat_f (Opposite direction)
%%---------------------------------------------------------
clear; clc; tic;
% pause(9999)
dbstop if error
nowpath = pwd;

% 开启并行环境
poolobj = gcp('nocreate'); % If no pool,  create new one.

if isempty(poolobj)
    MyPar = parpool(18);
else
    MyPar = gcp('nocreate');
    disp('Already initialized'); %说明并行环境已经启动。
end

%% 预选读取所有的路径
exm1 = 3; exm2 = 4;
mdl1 = 1; mdl2 = 'end';
esm1 = 1; esm2 = 'end';
esmPath_assmble = cell(2, 1);
mdlPath_assmble = esmPath_assmble;
exmNum_assmble = zeros(2, 1);
mdlNum_assmble = exmNum_assmble;
esmNum_assmble = exmNum_assmble;
esmCount = 0;

for exmNum = exm1:exm2%1 mean amip 2000; 2 mean amip 1980; 3 means ssp245, 4 means ssp370, 6 abrupt-4xCO2_150years
    % model parameters
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % experiment path (tLin:1740)
    inputPath = '/data1/liuyincheng/CMIP6-process/';
    exmPath_all = cell(1, length(Experiment));

    for varNum = 1:length(Experiment)
        exmPath_all{varNum} = fullfile(inputPath, level.time1{exmNum});
    end

    exmPath = exmPath_all{exmNum};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model
    if strcmp(mdl2, 'end')
        mdlR = length(level.model2);
    else
        mdlR = mdl2;
    end

    % mdlPath_assmble=cell((mdl2-mdl1+1)*(exm2-exm1+1),1);
    for mdlNum = mdl1:mdlR% model number
        % model path
        mdlPath = fullfile(exmPath, level.model2{mdlNum});
        % ensemble member path
        esmName = getPath_fileName(mdlPath, '.');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensemble member
        if strcmp(esm2, 'end')
            esmR = length(esmName);
        end

        for esmNum = esm1:esmR
            esmCount = esmCount + 1;
            esmPath = fullfile(mdlPath, esmName{esmNum, 1});
            esmPath_assmble{esmCount, 1} = esmPath;
            mdlPath_assmble{esmCount, 1} = mdlPath;
            exmNum_assmble(esmCount, 1) = exmNum;
            mdlNum_assmble(esmCount, 1) = mdlNum;
            esmNum_assmble(esmCount, 1) = esmNum;
        end

    end

end

parfor esmNum = 1:length(esmPath_assmble)
    esmFun(mdlPath_assmble{esmNum}, esmPath_assmble{esmNum}, exmNum_assmble(esmNum), mdlNum_assmble(esmNum), esmNum_assmble(esmNum));
end

delete(MyPar)

eval(['cd ', nowpath]);

t = toc; tickTok(t)

function [] = esmFun(mdlPath, esmPath, exmNum, mdlNum, esmNum)

    %% parameters and path
    Rv = 487.5;
    Lv = 2.5e6;
    %%% kernels in differnt conditions
    kernelsHightLabel = {'sfc', 'toa'}; kernelsLevel2_sky = {'cld', 'clr'};
    kernelsName = cell(1, 4);
    saveradEfectName = cell(1, 4);
    n = 0;

    for i = 1:2

        for j = 1:2
            n = n + 1;
            kernelsName{n} = ['kernels_', kernelsHightLabel{i}, '_', kernelsLevel2_sky{j}, '.mat'];
            saveradEfectName{n} = ['dradEffect_', kernelsHightLabel{i}, '_', kernelsLevel2_sky{j}, '.mat'];
        end

    end

    saveTrend_radEfectName = {'trend_dradEffect_sfc_cld.mat', 'trend_dradEffect_toa_cld.mat'};
    plevel = 24; % wv_lwkernel and wv_swkernel level;
    t_scflevel = 25; % sfc t_lwkernel level;
    load('/data1/liuyincheng/CMIP6-process/z_globalVar/ERF_rec.mat')% 'ERF_rec', 'timeERF_rec'
    lon_k = 0:2.5:357.5; nlonk = length(lon_k); % kernel lat lon
    lat_k = 90:-2.5:-90; nlatk = length(lat_k);
    lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
    lon_f = lon_k; nlonf = length(lon_f);
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    %%% load ERF data
    if strcmp(Experiment{exmNum}(1:3), 'ami')
        timeERF = timeERF_rec.hist(find(timeERF_rec.hist == tLin.startYear{exmNum}):find(timeERF_rec.hist == tLin.endYear{exmNum}));
        ERForcing = ERF_rec.hist(find(timeERF_rec.hist == tLin.startYear{exmNum}):find(timeERF_rec.hist == tLin.endYear{exmNum}));
    elseif strcmp(Experiment{exmNum}(1:6), 'ssp245')
        timeERF = timeERF_rec.ssp(find(timeERF_rec.ssp == tLin.startYear{exmNum}):find(timeERF_rec.ssp == tLin.endYear{exmNum}));
        ERForcing = ERF_rec.ssp245(find(timeERF_rec.ssp == tLin.startYear{exmNum}):find(timeERF_rec.ssp == tLin.endYear{exmNum}));
    elseif strcmp(Experiment{exmNum}(1:6), 'ssp370')
        timeERF = timeERF_rec.ssp(find(timeERF_rec.ssp == tLin.startYear{exmNum}):find(timeERF_rec.ssp == tLin.endYear{exmNum}));
        ERForcing = ERF_rec.ssp370(find(timeERF_rec.ssp == tLin.startYear{exmNum}):find(timeERF_rec.ssp == tLin.endYear{exmNum}));
    else
        ERForcing = 0;
        disp('this experient doesnt need ERF or havent input ERF Data!')
    end

    if length(ERForcing) ~= 1% if ERForcing=0 mean .eg. abrupt experint
        ERForcing = repmat(ERForcing, [12 1]);
        ERForcing = reshape(ERForcing, 1, length(ERForcing(:))); % yearx12, assume every month have the same radiative forcing
    end

    % ensemble member path
    esmName = getPath_fileName(mdlPath, '.');
    disp([level.model2{mdlNum}, ', ', level.time1{exmNum}, ', ', esmName{esmNum, 1}, ' ensemble start!'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ensemble member
    % path
    dvarsPath = fullfile(esmPath, level.process3{2}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/essember/anomaly
    varsPath = fullfile(esmPath, level.process3{1}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/essember/rawdata
    radEfectPath = fullfile(esmPath, level.process3{6}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/essember/radEffect/
    auto_mkdir(radEfectPath)
    %%%%% step1: cal dReffect
    % load dvars
    load([dvarsPath, 'global_vars.mat'])% lat_k lon_k time plev_k readme
    load([varsPath, 'hus.mat'])
    load([dvarsPath, 'dhus.mat'])% dhus and clim_hus
    load([dvarsPath, 'dalb.mat'])%
    load([dvarsPath, 'dta.mat'])%
    load([dvarsPath, 'dts.mat'])%
    dta(isnan(dta)) = 0;
    dalb(isnan(dalb)) = 0; dts(isnan(dts)) = 0;
    ntime = length(time.date);
    % cal the fixed moisture because kernel q's unit is W/m2/K/100mb
    startMonth = 1;
    dhus2 = zeros(144, 73, 24, ntime);

    for monNum = 1:ntime
        dhus2(:, :, :, monNum) = (log(hus(:, :, :, monNum)) - log(clim_hus(:, :, :, mod(monNum + startMonth - 2, 12) + 1))) .* clim_ta(:, :, :, mod(monNum + startMonth - 2, 12) + 1).^2 * Rv / Lv;
    end

    clear hus dhus clim_hus clim_ta clim_ts clim_alb

    dhus2(isnan(dhus2)) = 0;
    save([dvarsPath, 'dhus2.mat'], 'dhus2');
    % loop kernels in differnt conditions( 1sfc cld,2sfc clr,3toa cld,4toa clr)
    kernelPath = fullfile(esmPath, level.process3{5}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/kernelsCal
    % dR_tas=zeros(144,72,ntime,2);dR_taOnly=dR_tas;
    % store vars in one var( 1sfc cld,2sfc clr,3toa cld,4toa clr)
    
    dR_hus = zeros(144, 72, ntime, 4);
    dR_alb = dR_hus; dR_ts = dR_hus; dR_ta = dR_hus; dR_nonCloud = dR_hus; dR_nonCloudAndTs = dR_hus;
    lw_dR_hus = dR_hus; lw_dR_ta = dR_hus; lw_dR_ts = dR_hus; sw_dR_hus = dR_hus; sw_dR_alb = dR_hus;
    lw_nonCloud = dR_hus; sw_nonCloud = dR_hus;
    lw_nonCloudAndTs = dR_hus; sw_nonCloudAndTs = dR_hus;

    for property_of_LevelSKy = 1:4%( 1sfc cld,2sfc clr,3toa cld,4toa clr)
        load([kernelPath, kernelsName{property_of_LevelSKy}])% read kernel
        % prepare zero mixre
        lw_husEffect = zeros(144, 73, plevel, ntime);
        sw_husEffect = zeros(144, 73, plevel, ntime);
        tsEffect = zeros(144, 73, ntime);
        albEffect = zeros(144, 73, ntime);

        if property_of_LevelSKy <= 2% sfc
            taEffect = zeros(144, 73, t_scflevel, ntime);
            % consider 1 levels near sfc
            tasEffect = zeros(144, 73, ntime);
            taOnlyEffect = zeros(144, 73, plevel, ntime);
            % consider 2 levels near sfc
            tasEffect2 = taEffect;
            tasEffect1 = taEffect;
            % taOnlyEffect2 = zeros(144, 73, plevel-2, ntime);
            % t0Effect = zeros(144, 73, ntime);% mean sfc ta radeffect
        else
            taEffect = zeros(144, 73, plevel, ntime);
            % consider 2 levels near sfc
            tasEffect2 = taEffect;
            tasEffect1 = taEffect;
            % taOnlyEffect2 = zeros(144, 73, plevel - 2, ntime);
        end

        % Kernels and mat data are both from bottom to top
        % mod function used to process the data starts from 1 and should be consist with kernel(start from 1)
        startMonth = 1;

        for monNum = 1:ntime
            lw_husEffect(:, :, :, monNum) = wv_lwkernel(:, :, :, mod(monNum + startMonth - 2, 12) + 1) .* dhus2(:, :, :, monNum);
            sw_husEffect(:, :, :, monNum) = wv_swkernel(:, :, :, mod(monNum + startMonth - 2, 12) + 1) .* dhus2(:, :, :, monNum);
            tsEffect(:, :, monNum) = ts_lwkernel(:, :, mod(monNum + startMonth - 2, 12) + 1) .* dts(:, :, monNum);
            albEffect(:, :, monNum) = alb_swkernel(:, :, mod(monNum + startMonth - 2, 12) + 1) .* dalb(:, :, monNum) * 100; % note that alb kernerl unite: W/m2/0.01

            if property_of_LevelSKy <= 2% sfc
                taEffect(:, :, 2:end, monNum) = t_lwkernel(:, :, 2:end, mod(monNum + startMonth - 2, 12) + 1) .* dta(:, :, :, monNum);
                taEffect(:, :, 1, monNum) = squeeze(t_lwkernel(:, :, 1, mod(monNum + startMonth - 2, 12) + 1)) .* dts(:, :, monNum);
                tasEffect(:, :, monNum) = taEffect(:, :, 1, monNum);
                taOnlyEffect(:, :, :, monNum) = taEffect(:, :, 2:end, monNum);

                % consider 1 levels near sfc
                tasEffect1(:, :, 2:end, monNum) = t_level1_lwkernel(:, :, 2:end, mod(monNum + startMonth - 2, 12) + 1) .* dta(:, :, :, monNum);
                tasEffect1(:, :, 1, monNum) = squeeze(t_level1_lwkernel(:, :, 1, mod(monNum + startMonth - 2, 12) + 1)) .* dts(:, :, monNum);

                % consider 2 levels near sfc
                tasEffect2(:, :, 2:end, monNum) = t_level2_lwkernel(:, :, 2:end, mod(monNum + startMonth - 2, 12) + 1) .* dta(:, :, :, monNum);
                tasEffect2(:, :, 1, monNum) = squeeze(t_level2_lwkernel(:, :, 1, mod(monNum + startMonth - 2, 12) + 1)) .* dts(:, :, monNum);
                else% toa

                taEffect(:, :, :, monNum) = t_lwkernel(:, :, :, mod(monNum + startMonth - 2, 12) + 1) .* dta(:, :, :, monNum);
                % consider 1 levels near sfc
                tasEffect1(:, :, :, monNum) = t_level1_lwkernel(:, :, :, mod(monNum + startMonth - 2, 12) + 1) .* dta(:, :, :, monNum);
                % consider 2 levels near sfc
                tasEffect2(:, :, :, monNum) = t_level2_lwkernel(:, :, :, mod(monNum + startMonth - 2, 12) + 1) .* dta(:, :, :, monNum);
            end

        end


        lw_husEffect = squeeze(nansum(lw_husEffect(:, :, :, :), 3));
        sw_husEffect = squeeze(nansum(sw_husEffect(:, :, :, :), 3));
        taEffect = squeeze(nansum(taEffect(:, :, :, :), 3));
        tasEffect2 = squeeze(nansum(tasEffect2(:, :, :, :), 3));
        tasEffect1 = squeeze(nansum(tasEffect1(:, :, :, :), 3));

        % regrid to figure lat and lon 144x72
        lw_husEffect = autoRegrid3(lon_k, lat_k, time.date, lw_husEffect, lon_f, lat_f, time.date);
        sw_husEffect = autoRegrid3(lon_k, lat_k, time.date, sw_husEffect, lon_f, lat_f, time.date);
        taEffect = autoRegrid3(lon_k, lat_k, time.date, taEffect, lon_f, lat_f, time.date);
        tsEffect = autoRegrid3(lon_k, lat_k, time.date, tsEffect, lon_f, lat_f, time.date);
        albEffect = autoRegrid3(lon_k, lat_k, time.date, albEffect, lon_f, lat_f, time.date);
        tasEffect2 = autoRegrid3(lon_k, lat_k, time.date, tasEffect2, lon_f, lat_f, time.date);
        tasEffect1 = autoRegrid3(lon_k, lat_k, time.date, tasEffect1, lon_f, lat_f, time.date);

        % integration(taEffect,tsEffect,albEffect,nonCloudEffect,husEffect)
        husEffect = lw_husEffect + sw_husEffect;
        nonCloudEffect = husEffect + taEffect + tsEffect + albEffect;
        nonCloudAndTsEffect = husEffect + taEffect + albEffect; % prepare for cal ta+alb+q+clod effect
        % sw and lw radEffect
        lw_taEffect = taEffect;
        lw_tsEffect = tsEffect;
        sw_albEffect = albEffect;
        lw_nonCloudEffect = lw_husEffect + lw_taEffect + lw_tsEffect;
        sw_nonCloudEffect = sw_husEffect + sw_albEffect;
        lw_nonCloudAndTsEffect = lw_husEffect + lw_taEffect;
        sw_nonCloudAndTsEffect = sw_husEffect + sw_albEffect;

        if property_of_LevelSKy <= 2
            taOnlyEffect = squeeze(nansum(taOnlyEffect(:, :, :, :), 3));
            % regrid special values to figure lat and lon 144x72
            taOnlyEffect = autoRegrid3(lon_k, lat_k, time.date, taOnlyEffect, lon_f, lat_f, time.date);
            tasEffect = autoRegrid3(lon_k, lat_k, time.date, tasEffect, lon_f, lat_f, time.date);
            % save the radEffect: dradEffect_sfc_cld.mat
            save([radEfectPath, 'global_vars.mat'], 'lon_f', 'lat_f', 'time', 'plev_k', 'readme', 'timeseries', 'modelname')
            save([radEfectPath, saveradEfectName{property_of_LevelSKy}], 'lw_husEffect', 'lw_taEffect', 'lw_tsEffect', 'sw_husEffect', 'sw_albEffect', ...
                'tsEffect', 'albEffect', 'husEffect', ...
                'taEffect', 'taOnlyEffect', 'tasEffect', 'tasEffect2', 'tasEffect1', ...
                'nonCloudEffect', 'lw_nonCloudEffect', 'sw_nonCloudEffect', ...
                'nonCloudAndTsEffect', 'lw_nonCloudAndTsEffect', 'sw_nonCloudAndTsEffect');
        else
            % save the radEffect: dradEffect_sfc_cld.mat
            save([radEfectPath, 'global_vars.mat'], 'lon_f', 'lat_f', 'time', 'plev_k', 'readme', 'timeseries', 'modelname')
            save([radEfectPath, saveradEfectName{property_of_LevelSKy}], 'lw_husEffect', 'lw_taEffect', 'lw_tsEffect', 'sw_husEffect', 'sw_albEffect', ...
                'tsEffect', 'albEffect', 'husEffect', ...
                'taEffect', 'tasEffect2', 'tasEffect1', ...
                'nonCloudEffect', 'lw_nonCloudEffect', 'sw_nonCloudEffect', ...
                'nonCloudAndTsEffect', 'lw_nonCloudAndTsEffect', 'sw_nonCloudAndTsEffect');
        end

        clear tasEffect2 tasEffect1

        dR_hus(:, :, :, property_of_LevelSKy) = husEffect; % q
        lw_dR_hus(:, :, :, property_of_LevelSKy) = lw_husEffect; % lw_q
        sw_dR_hus(:, :, :, property_of_LevelSKy) = sw_husEffect; % sw_q
        clear husEffect lw_husEffect sw_husEffect

        dR_alb(:, :, :, property_of_LevelSKy) = albEffect; %albedo
        sw_dR_alb(:, :, :, property_of_LevelSKy) = sw_albEffect; %sw_albedo
        clear albEffect sw_albEffect

        dR_ts(:, :, :, property_of_LevelSKy) = tsEffect; % ts
        lw_dR_ts(:, :, :, property_of_LevelSKy) = lw_tsEffect; % lw_ts
        clear tsEffect lw_tsEffect

        dR_ta(:, :, :, property_of_LevelSKy) = taEffect; % ta
        lw_dR_ta(:, :, :, property_of_LevelSKy) = lw_taEffect; % lw_ta
        clear taEffect lw_taEffect

        dR_nonCloud(:, :, :, property_of_LevelSKy) = nonCloudEffect; % nonCloudEffect
        lw_nonCloud(:, :, :, property_of_LevelSKy) = lw_nonCloudEffect; % lw_nonCloudEffect
        sw_nonCloud(:, :, :, property_of_LevelSKy) = sw_nonCloudEffect; % sw_nonCloudEffect
        dR_nonCloudAndTs(:, :, :, property_of_LevelSKy) = nonCloudAndTsEffect; % prepare for cal ta+alb+q+clod effect
        lw_nonCloudAndTs(:, :, :, property_of_LevelSKy) = lw_nonCloudEffect; % lw_nonCloudAndTsEffect
        sw_nonCloudAndTs(:, :, :, property_of_LevelSKy) = sw_nonCloudEffect; % sw_nonCloudAndTsEffect
        clear nonCloudEffect lw_nonCloudEffect sw_nonCloudEffect nonCloudAndTsEffect lw_nonCloudAndTsEffect sw_nonCloudAndTsEffect

    end

    clear dta dts dalb dhus2

    readme_radEfect = 'final column: 1sfc cld,2sfc clr,3toa cld,4toa clr, note that dR_tas and dR_taOnly shouldnt have value in toa';
    save([radEfectPath, 'dradEffect_union.mat'], 'dR_hus', 'dR_alb', 'dR_ts', 'dR_ta', 'dR_nonCloud', 'readme_radEfect', ...
        'lw_dR_hus', 'lw_dR_ta', 'lw_dR_ts', 'sw_dR_hus', 'sw_dR_alb', ...
        'lw_nonCloud', 'sw_nonCloud', 'lw_nonCloudAndTs', 'sw_nonCloudAndTs');

    %%%%% step2: cal cloud effect(definen down is postive)
    % load origan model rad dvars
    % 1.sfc dvars(includ cld and clr sky)
    load([dvarsPath, 'drlus.mat'])% surface_upwelling_longwave_flux_in_air
    load([dvarsPath, 'drsus.mat']); load([dvarsPath, 'drsuscs.mat'])% surface_upwelling_shortwave_flux_in_air
    load([dvarsPath, 'drsds.mat']); load([dvarsPath, 'drsdscs.mat'])% surface_downwelling_shortwave_flux_in_air
    load([dvarsPath, 'drlds.mat']); load([dvarsPath, 'drldscs.mat'])% surface_downwelling_longwave_flux_in_air
    % 2.toa dvars(includ cld and clr sky)
    load([dvarsPath, 'drsdt.mat'])% toa_incoming_shortwave_flux
    load([dvarsPath, 'drsut.mat']); load([dvarsPath, 'drsutcs.mat'])% toa_outgoing_shortwave_flux
    load([dvarsPath, 'drlut.mat']); load([dvarsPath, 'drlutcs.mat'])% toa_outgoing_longwave_flux

    % model output rad( 1sfc cld,2sfc clr,3toa cld,4toa clr)
    lw_net14473(:, :, :, 1) = drlds - drlus; % Sfc_cld net thermal radiation,
    lw_net14473(:, :, :, 2) = drldscs - drlus; % Sfc_clr net thermal radiation,
    lw_net14473(:, :, :, 3) = -drlut; % toa_cld net thermal radiation
    lw_net14473(:, :, :, 4) = -drlutcs; % toa_clr net thermal radiation

    sw_net14473(:, :, :, 1) = drsds - drsus; % Sfc_cld net solar radiation,
    sw_net14473(:, :, :, 2) = drsdscs - drsuscs; % Sfc_clr net solar radiation,
    sw_net14473(:, :, :, 3) = drsdt - drsut; % toa_cld net solar radiation,
    sw_net14473(:, :, :, 4) = drsdt - drsutcs; % toa_clr net solar radiation,
    lw_net14473(isnan(lw_net14473)) = 0;
    sw_net14473(isnan(sw_net14473)) = 0;

    %regrid 144x72(unite grids)
    lw_net = zeros(144, 72, ntime, 4); sw_net = lw_net;

    for property_of_LevelSKy = 1:4
        lw_net(:, :, :, property_of_LevelSKy) = autoRegrid3(lon_k, lat_k, time.date, squeeze(lw_net14473(:, :, :, property_of_LevelSKy)), lon_f, lat_f, time.date);
        sw_net(:, :, :, property_of_LevelSKy) = autoRegrid3(lon_k, lat_k, time.date, squeeze(sw_net14473(:, :, :, property_of_LevelSKy)), lon_f, lat_f, time.date);
    end

    clear lw_net14473 sw_net14473

    % 1.sfc, 2.toa( 1sfc cld,2sfc clr,3toa cld,4toa clr)
    dR_net_cld = lw_net(:, :, :, [1 3]) + sw_net(:, :, :, [1 3]); % real net rad
    dR_net_clr = lw_net(:, :, :, [2 4]) + sw_net(:, :, :, [2 4]);
    % save real radEffect: eg:real_dradEffect.mat
    readme_realradEfect = {'final row: 1sfc cld,2sfc clr,3toa cld,4toa clr', 'dR_net_cld,dR_net_clr:1.sfc, 2.toa'};
    save([radEfectPath, 'model_dradEffect.mat'], 'lw_net', 'sw_net', 'dR_net_cld', 'dR_net_clr', 'readme_realradEfect');

    %% cal the dRcloud and dR_co2(144,72,time,2)(final row means contain sfc and toa)
    % note that residual includ( 1sfc cld,2sfc clr,3toa cld,4toa clr)
    % dCRF=dR_net_cld-dR_net_clr;%(1sfc 2toa)
    % method 1 of cal cloud effect
    dvarsFeedback = dR_nonCloud(:, :, :, [2 4]) - dR_nonCloud(:, :, :, [1 3]);
    lw_dvarsFeedback = lw_nonCloud(:, :, :, [2 4]) - lw_nonCloud(:, :, :, [1 3]);
    sw_dvarsFeedback = sw_nonCloud(:, :, :, [2 4]) - sw_nonCloud(:, :, :, [1 3]);

    dCRF = dR_net_cld - dR_net_clr; % cloud radiative forcing
    lw_dCRF = lw_net(:, :, :, [1 3]) - lw_net(:, :, :, [2 4]);
    sw_dCRF = sw_net(:, :, :, [1 3]) - sw_net(:, :, :, [2 4]);

    % dCRF = dR_net_clr - dR_net_cld; % cloud radiative forcing
    dR_cloud = dCRF + dvarsFeedback;
    lw_cloud = lw_dCRF + lw_dvarsFeedback;
    sw_cloud = sw_dCRF + sw_dvarsFeedback;

    % add ERF (G0-G)/G~0.16 according soden paper
    if length(ERForcing) ~= 1

        for timeLine = 1:ntime
            dR_cloud(:, :, timeLine, 2) = dR_cloud(:, :, timeLine, 2) + ERForcing(timeLine) .* 0.16;
            lw_cloud(:, :, timeLine, 2) = lw_cloud(:, :, timeLine, 2) + ERForcing(timeLine) .* 0.16;
        end

    end

    dR_residual_clr = dR_net_clr - dR_nonCloud(:, :, :, [2 4]);
    lw_residual_clr = lw_net(:, :, :, [2 4]) - lw_nonCloud(:, :, :, [2 4]);
    sw_residual_clr = sw_net(:, :, :, [2 4]) - sw_nonCloud(:, :, :, [2 4]);

    dR_residual_cld = dR_net_cld - dR_nonCloud(:, :, :, [1 3]) - dR_cloud; % r terms, indicate co2/aerosol forcing
    lw_residual_cld = lw_net(:, :, :, [1 3]) - lw_nonCloud(:, :, :, [1 3]) - lw_cloud;
    sw_residual_cld = sw_net(:, :, :, [1 3]) - sw_nonCloud(:, :, :, [1 3]) - sw_cloud;

    dR_residual_cld1 = dR_net_cld - dR_nonCloud(:, :, :, [1 3]); % r terms, indicate cloud feedback and co2/aerosol forcing

    % method 2 of cal cloud effect
    % dR_cloud = dR_residual_cld1 - dR_residual_clr;

    dR_cloud_sfc = dR_cloud(:, :, :, 1); dR_cloud_toa = dR_cloud(:, :, :, 2);
    lw_cloud_sfc = lw_cloud(:, :, :, 1); lw_cloud_toa = lw_cloud(:, :, :, 2);
    sw_cloud_sfc = sw_cloud(:, :, :, 1); sw_cloud_toa = sw_cloud(:, :, :, 2);

    dCRF_sfc = dCRF(:, :, :, 1); dCRF_toa = dCRF(:, :, :, 2);
    lw_dCRF_sfc = dCRF(:, :, :, 1); lw_dCRF_toa = lw_dCRF(:, :, :, 2);
    sw_dCRF_sfc = dCRF(:, :, :, 1); sw_dCRF_toa = sw_dCRF(:, :, :, 2);

    dR_residual_cld_sfc = dR_residual_cld(:, :, :, 1); dR_residual_cld_toa = dR_residual_cld(:, :, :, 2);
    dR_residual_clr_sfc = dR_residual_clr(:, :, :, 1); dR_residual_clr_toa = dR_residual_clr(:, :, :, 2);

    lw_residual_cld_sfc = lw_residual_cld(:, :, :, 1); lw_residual_cld_toa = lw_residual_cld(:, :, :, 2);
    lw_residual_clr_sfc = lw_residual_clr(:, :, :, 1); lw_residual_clr_toa = lw_residual_clr(:, :, :, 2);

    sw_residual_cld_sfc = sw_residual_cld(:, :, :, 1); sw_residual_cld_toa = sw_residual_cld(:, :, :, 2);
    sw_residual_clr_sfc = sw_residual_clr(:, :, :, 1); sw_residual_clr_toa = sw_residual_clr(:, :, :, 2);

    % equal to real_dradEffect.mat
    dR_net_cld_sfc = dR_net_cld(:, :, :, 1); dR_net_cld_toa = dR_net_cld(:, :, :, 2);
    dR_net_clr_sfc = dR_net_clr(:, :, :, 1); dR_net_clr_toa = dR_net_clr(:, :, :, 2);
    % save dR_cloud and dR_residual_cld and observe total effect
    save([radEfectPath, 'dCRF.mat'], 'dCRF', 'lw_dCRF', 'sw_dCRF', ...%CRF_sfc
    'dCRF_sfc', 'lw_dCRF_sfc', 'sw_dCRF_sfc', ...% CRF_sfc
    'dCRF_toa', 'lw_dCRF_toa', 'sw_dCRF_toa'); % CRF_toa
    save([radEfectPath, 'dvarsFeedback.mat'], 'dvarsFeedback', 'lw_dvarsFeedback', 'sw_dvarsFeedback');
    clear dvarsFeedback lw_dvarsFeedback sw_dvarsFeedback

    save([radEfectPath, 'dR_cloud.mat'], 'dR_cloud_sfc', 'lw_cloud_sfc', 'sw_cloud_sfc', ...
        'dR_cloud_toa', 'lw_cloud_toa', 'sw_cloud_toa');

    save([radEfectPath, 'dR_residual_cld_sfc.mat'], 'dR_residual_cld_sfc', 'lw_residual_cld_sfc', 'sw_residual_cld_sfc');
    save([radEfectPath, 'dR_residual_cld_toa.mat'], 'dR_residual_cld_toa', 'lw_residual_cld_toa', 'sw_residual_cld_toa');
    save([radEfectPath, 'dR_residual_clr_sfc.mat'], 'dR_residual_clr_sfc', 'lw_residual_clr_sfc', 'sw_residual_clr_sfc');
    save([radEfectPath, 'dR_residual_clr_toa.mat'], 'dR_residual_clr_toa', 'lw_residual_clr_toa', 'sw_residual_clr_toa');
    save([radEfectPath, 'dR_net.mat'], 'dR_net_cld_sfc', 'dR_net_cld_toa', ...
        'dR_net_clr_sfc', 'dR_net_clr_toa');

    %%%%% step3: cal The effect of Ts on the atmosphere(dR_tstoa-dR_tssfc) (note effect on atoms is positive, but effect on earth is nagetive)
    % and nonTs Effect (ta+alb+q+clod effect)
    %( 1sfc cld,2sfc clr,3toa cld,4toa clr)
    dR_tsAtom_cld = dR_ts(:, :, :, 3) - dR_ts(:, :, :, 1);
    dR_tsAtom_clr = dR_ts(:, :, :, 4) - dR_ts(:, :, :, 2);
    save([radEfectPath, 'dR_tsAtom_cld.mat'], 'dR_tsAtom_cld');
    save([radEfectPath, 'dR_tsAtom_clr.mat'], 'dR_tsAtom_clr');
    dR_nonTs_sfc = dR_cloud_sfc + dR_nonCloudAndTs(:, :, :, 1);
    dR_nonTs_toa = dR_cloud_toa + dR_nonCloudAndTs(:, :, :, 3);
    save([radEfectPath, 'dR_nonTs_sfc.mat'], 'dR_nonTs_sfc');
    save([radEfectPath, 'dR_nonTs_toa.mat'], 'dR_nonTs_toa');
    % 20/8/9 补充 cal total effect
    dR_total_sfc = dR_cloud_sfc + dR_nonCloud(:, :, :, 1);
    dR_total_toa = dR_cloud_toa + dR_nonCloud(:, :, :, 3);
    save([radEfectPath, 'dR_total_sfc.mat'], 'dR_total_sfc');
    save([radEfectPath, 'dR_total_toa.mat'], 'dR_total_toa');

    disp([level.model2{mdlNum}, ', ', level.time1{exmNum}, ', ', esmName{esmNum, 1}, ' ensemble is done!'])
    disp(' ')

end
