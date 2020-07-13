%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-07-06 08:55:44
% LastEditTime : 2020-07-13 09:31:46
% LastEditors  : LYC
% Description  :
% FilePath     : /code/p1_processObserveData/ERAi/ERAi_s3_calRadEff_dp.m
%
%%---------------------------------------------------------
clear; clc; tic;
% pause(9999)
%%% constant
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
        saveradEfectName{n} = ['dradEfect_', kernelsHightLabel{i}, '_', kernelsLevel2_sky{j}, '.mat'];
    end

end

saveTrend_radEfectName = {'trend_dradEfect_sfc_cld.mat', 'trend_dradEfect_toa_cld.mat'};
plevel = 24; % wv_lwkernel and wv_swkernel level;
t_scflevel = 25; % sfc t_lwkernel level;

[readme, level, tLin, vars] = obsParameters('ERAi');
for p_1 = 1:2
    varsPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi', level.standVarPath{1}); %rawdata
    dvarsPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi', level.standVarPath{2}); %anomaly
    kernelCalPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi', level.standVarPath{4}); % kernelCal
    radEfectPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi', level.standVarPath{5}); %radEffect
    auto_mkdir(radEfectPath)
    % load dvars
    load([dvarsPath, 'global_vars.mat'])% lat_k lon_k time plev_k
    load([dvarsPath, 'meto_dvars.mat'])% dhus and clim_hus...
    load([varsPath, 'meto_vars.mat'], 'hus')%
    dta(isnan(dta)) = 0;
    dalb(isnan(dalb)) = 0; dts(isnan(dts)) = 0;
    ntime = length(time);
    lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f);
    lon_f = lon_k; nlonf = length(lon_f);
    startMonth = tLin.startMonth{p_1};

    %%%%% step1: cal dReffect
    % cal the fixed moisture because kernel q's unit is W/m2/K/100mb
    dhus2 = zeros(144, 73, 24, ntime);

    for monNum = 1:ntime
        dhus2(:, :, :, monNum) = (log(hus(:, :, :, monNum)) - log(clim_hus(:, :, :, mod(monNum + startMonth - 2, 12) + 1))) .* clim_ta(:, :, :, mod(monNum + startMonth - 2, 12) + 1).^2 * Rv / Lv;
    end

    dhus2(isnan(dhus2)) = 0;
    save([dvarsPath, 'dhus2.mat'], 'dhus2');

    % loop kernels in differnt conditions( 1sfc cld,2sfc clr,3toa cld,4toa clr)
    dR_hus = zeros(144, 72, ntime, 4);
    dR_alb = dR_hus; dR_ts = dR_hus; dR_ta = dR_hus; dR_total = dR_hus; dR_main = dR_hus;
    % dR_tas=zeros(144,72,ntime,2);dR_taOnly=dR_tas;

    for ii = 1:4%( 1sfc cld,2sfc clr,3toa cld,4toa clr)
        load([kernelCalPath, kernelsName{ii}])% read kernel
        % prepare zero mixre
        wvlwEffect = zeros(144, 73, plevel, ntime);
        wvswEffect = zeros(144, 73, plevel, ntime);
        tsEffect = zeros(144, 73, ntime);
        albEffect = zeros(144, 73, ntime);

        if ii <= 2% sfc
            taEffect = zeros(144, 73, t_scflevel, ntime);
            % consider 1 levels near sfc
            tasEffect = zeros(144, 73, ntime);
            taOnlyEffect = zeros(144, 73, plevel, ntime);
            % consider 2 levels near sfc
            tasEffect2 = taEffect;
        else
            taEffect = zeros(144, 73, plevel, ntime);
            % consider 2 levels near sfc
            tasEffect2 = taEffect;
        end

        % Kernels and mat data are both from bottom to top
        % mod function used to process the data starts from 1 and should be consist with kernel(start from 1)

        for monNum = 1:ntime
            wvlwEffect(:, :, :, monNum) = wv_lwkernel(:, :, :, mod(monNum + startMonth - 2, 12) + 1) .* dhus2(:, :, :, monNum);
            wvswEffect(:, :, :, monNum) = wv_swkernel(:, :, :, mod(monNum + startMonth - 2, 12) + 1) .* dhus2(:, :, :, monNum);
            tsEffect(:, :, monNum) = ts_lwkernel(:, :, mod(monNum + startMonth - 2, 12) + 1) .* dts(:, :, monNum);
            albEffect(:, :, monNum) = alb_swkernel(:, :, mod(monNum + startMonth - 2, 12) + 1) .* dalb(:, :, monNum) * 100; % note that alb kernerl unite: W/m2/0.01

            if ii <= 2% sfc
                taEffect(:, :, 2:end, monNum) = t_lwkernel(:, :, 2:end, mod(monNum + startMonth - 2, 12) + 1) .* dta(:, :, :, monNum);
                taEffect(:, :, 1, monNum) = squeeze(t_lwkernel(:, :, 1, mod(monNum + startMonth - 2, 12) + 1)) .* dts(:, :, monNum);
                tasEffect(:, :, monNum) = taEffect(:, :, 1, monNum);
                taOnlyEffect(:, :, :, monNum) = taEffect(:, :, 2:end, monNum);
                % consider 2 levels near sfc(have problem)
                tasEffect2(:, :, 2:end, monNum) = t_level2_lwkernel(:, :, 2:end, mod(monNum + startMonth - 2, 12) + 1) .* dta(:, :, :, monNum);
                tasEffect2(:, :, 1, monNum) = squeeze(t_level2_lwkernel(:, :, 1, mod(monNum + startMonth - 2, 12) + 1)) .* dts(:, :, monNum);
                else% toa
                taEffect(:, :, :, monNum) = t_lwkernel(:, :, :, mod(monNum + startMonth - 2, 12) + 1) .* dta(:, :, :, monNum);
                % consider 2 levels near sfc
                tasEffect2(:, :, :, monNum) = t_level2_lwkernel(:, :, :, mod(monNum + startMonth - 2, 12) + 1) .* dta(:, :, :, monNum);
            end

        end

        taOnlyEffect2 = taEffect - tasEffect2;

        wvlwEffect = squeeze(nansum(wvlwEffect(:, :, :, :), 3));
        wvswEffect = squeeze(nansum(wvswEffect(:, :, :, :), 3));
        taEffect = squeeze(nansum(taEffect(:, :, :, :), 3));
        tasEffect2 = squeeze(nansum(tasEffect2(:, :, :, :), 3));

        if ii <= 2
            taOnlyEffect = squeeze(nansum(taOnlyEffect(:, :, :, :), 3));
            taOnlyEffect2 = squeeze(nansum(taOnlyEffect2(:, :, :, :), 3));
        else
            taOnlyEffect2 = squeeze(nansum(taOnlyEffect2(:, :, :, :), 3));
        end

        % regrid 144x72(unite grids)
        wvlwEffect = autoRegrid3(lon_k, lat_k, time, wvlwEffect, lon_f, lat_f, time);
        wvswEffect = autoRegrid3(lon_k, lat_k, time, wvswEffect, lon_f, lat_f, time);
        taEffect = autoRegrid3(lon_k, lat_k, time, taEffect, lon_f, lat_f, time);
        tsEffect = autoRegrid3(lon_k, lat_k, time, tsEffect, lon_f, lat_f, time);
        albEffect = autoRegrid3(lon_k, lat_k, time, albEffect, lon_f, lat_f, time);

        % integration(taEffect,tsEffect,albEffect,totalEffect,husEffect)
        husEffect = wvlwEffect + wvswEffect;
        totalEffect = husEffect + taEffect + tsEffect + albEffect;
        mainEffect = husEffect + taEffect + albEffect; % prepare for cal ta+alb+q+clod effect

        if ii <= 2
            taOnlyEffect = autoRegrid3(lon_k, lat_k, time, taOnlyEffect, lon_f, lat_f, time);
            taOnlyEffect2 = autoRegrid3(lon_k, lat_k, time, taOnlyEffect2, lon_f, lat_f, time);
            tasEffect = autoRegrid3(lon_k, lat_k, time, tasEffect, lon_f, lat_f, time);
            tasEffect2 = autoRegrid3(lon_k, lat_k, time, tasEffect2, lon_f, lat_f, time);

            % save the radEffect: dradEfect_sfc_cld.mat
            save([radEfectPath, 'global_vars.mat'], 'lon_f', 'lat_f', 'lon_k', 'lat_k', 'plev_k', 'time')
            save([radEfectPath, saveradEfectName{ii}], 'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', ...
                'taEffect', 'taOnlyEffect', 'tasEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect', 'mainEffect');
        else
            taOnlyEffect2 = autoRegrid3(lon_k, lat_k, time, taOnlyEffect2, lon_f, lat_f, time);
            tasEffect2 = autoRegrid3(lon_k, lat_k, time, tasEffect2, lon_f, lat_f, time);
            % save the radEffect: dradEfect_sfc_cld.mat
            save([radEfectPath, 'global_vars.mat'], 'lon_f', 'lat_f', 'lon_k', 'lat_k', 'plev_k', 'time')
            save([radEfectPath, saveradEfectName{ii}], 'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', ...
                'taEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect', 'mainEffect');
        end

        % store vars in one var( 1sfc cld,2sfc clr,3toa cld,4toa clr)
        dR_hus(:, :, :, ii) = husEffect; % q
        dR_alb(:, :, :, ii) = albEffect; %albedo
        dR_ts(:, :, :, ii) = tsEffect; % ts
        dR_ta(:, :, :, ii) = taEffect; % t
        dR_total(:, :, :, ii) = totalEffect; % total
        dR_main(:, :, :, ii) = mainEffect; % prepare for cal ta+alb+q+clod effect
    end

    readme_radEfect = 'final column: 1sfc cld,2sfc clr,3toa cld,4toa clr, note that dR_tas and dR_taOnly shouldnt have value in toa';
    save([radEfectPath, 'dradEfect_union.mat'], 'dR_hus', 'dR_alb', 'dR_ts', 'dR_ta', 'dR_total', 'readme_radEfect');

    %%%%% step2: cal cloud effect(definen down is postive)
    % load origan model rad dvars
    % 1.sfc dvars(includ cld and clr sky)
    load([dvarsPath, 'dstr.mat'])
    load([dvarsPath, 'dstrc.mat'])
    load([dvarsPath, 'dttr.mat'])
    load([dvarsPath, 'dttrc.mat'])
    ll_rad(:, :, :, 1) = dstr; % Sfc_cld net thermal radiation,
    ll_rad(:, :, :, 2) = dstrc; % Sfc_clr net thermal radiation,
    ll_rad(:, :, :, 3) = dttr; % toa_cld net thermal radiation
    ll_rad(:, :, :, 4) = dttrc; % toa_clr net thermal radiation

    load([dvarsPath, 'dssr.mat'])
    load([dvarsPath, 'dssrc.mat'])
    load([dvarsPath, 'dtsr.mat'])
    load([dvarsPath, 'dtsrc.mat'])
    ss_rad(:, :, :, 1) = dssr; % Sfc_cld net solar radiation,
    ss_rad(:, :, :, 2) = dssrc; % Sfc_clr net solar radiation,
    ss_rad(:, :, :, 3) = dtsr; % toa_cld net solar radiation,
    ss_rad(:, :, :, 4) = dtsrc; % toa_clr net solar radiation,

    ll_rad(isnan(ll_rad)) = 0;
    ss_rad(isnan(ss_rad)) = 0;

    %regrid 144x72(unite grids)
    for ii = 1:4
        l_rad(:, :, :, ii) = autoRegrid3(lon_k, lat_k, time, squeeze(ll_rad(:, :, :, ii)), lon_f, lat_f, time);
        s_rad(:, :, :, ii) = autoRegrid3(lon_k, lat_k, time, squeeze(ss_rad(:, :, :, ii)), lon_f, lat_f, time);
    end
    
    % 1.sfc, 2.toa( 1sfc cld,2sfc clr,3toa cld,4toa clr)
    dR_allsky = l_rad(:, :, :, [1 3]) + s_rad(:, :, :, [1 3]); % real net rad
    dR_clr = l_rad(:, :, :, [2 4]) + s_rad(:, :, :, [2 4]);
    % save real radEffect: eg:real_dradEfect.mat
    readme_realradEfect = {'final row: 1sfc cld,2sfc clr,3toa cld,4toa clr', 'dR_allsky,dR_clr:1.sfc, 2.toa'};
    save([radEfectPath, 'real_dradEfect.mat'], 'l_rad', 's_rad', 'dR_allsky', 'dR_clr', 'readme_realradEfect');
    clear ll_rad ss_rad l_rad s_rad
    %% cal the dRcloud and dR_co2(144,72,time,2)(final row means contain sfc and toa)
    % note that residual includ( 1sfc cld,2sfc clr,3toa cld,4toa clr)
    % dCRF=dR_allsky-dR_clr;%(1sfc 2toa)
    %method 1
    dvarsFeedback = dR_total(:, :, :, [2 4]) - dR_total(:, :, :, [1 3]);
    dCRF = dR_allsky - dR_clr; % cloud radiative forcing
    % dCRF = dR_clr - dR_allsky; % cloud radiative forcing
    dR_cloud = dCRF + dvarsFeedback;
    %method 2
    dR_residual_clr = dR_clr - dR_total(:, :, :, [2 4]);
    dR_residual_cld = dR_allsky - dR_total(:, :, :, [1 3]) - dR_cloud; % r terms, indicate co2/aerosol forcing
    dR_residual_cld1 = dR_allsky - dR_total(:, :, :, [1 3]); % r terms, indicate cloud feedback and co2/aerosol forcing
    % dR_cloud = dR_residual_cld1 - dR_residual_clr;

    dR_cloud_sfc = dR_cloud(:, :, :, 1); dR_cloud_toa = dR_cloud(:, :, :, 2);
    dCRF_sfc = dCRF(:, :, :, 1); dCRF_toa = dCRF(:, :, :, 2);
    dR_residual_cld_sfc = dR_residual_cld(:, :, :, 1); dR_residual_cld_toa = dR_residual_cld(:, :, :, 2);
    dR_residual_clr_sfc = dR_residual_clr(:, :, :, 1); dR_residual_clr_toa = dR_residual_clr(:, :, :, 2);
    % equal to real_dradEfect.mat
    dR_ObsTotal_cld_sfc = dR_allsky(:, :, :, 1); dR_ObsTotal_cld_toa = dR_allsky(:, :, :, 2);
    dR_ObsTotal_clr_sfc = dR_clr(:, :, :, 1); dR_ObsTotal_clr_toa = dR_clr(:, :, :, 2);
    % save dR_cloud and dR_residual_cld and observe total effect
    save([radEfectPath, 'dCRF.mat'], 'dCRF');
    save([radEfectPath, 'dCRF_sfc.mat'], 'dCRF_sfc');
    save([radEfectPath, 'dCRF_toa.mat'], 'dCRF_toa');
    save([radEfectPath, 'dR_cloud_sfc.mat'], 'dR_cloud_sfc');
    save([radEfectPath, 'dR_cloud_toa.mat'], 'dR_cloud_toa');
    save([radEfectPath, 'dR_residual_cld_sfc.mat'], 'dR_residual_cld_sfc'); % note that residual cld=clr(in sfc or toa);
    save([radEfectPath, 'dR_residual_cld_toa.mat'], 'dR_residual_cld_toa');
    save([radEfectPath, 'dR_residual_clr_sfc.mat'], 'dR_residual_clr_sfc');
    save([radEfectPath, 'dR_residual_clr_toa.mat'], 'dR_residual_clr_toa');
    save([radEfectPath, 'dR_ObsTotal.mat'], 'dR_ObsTotal_cld_sfc', 'dR_ObsTotal_cld_toa', ...
        'dR_ObsTotal_clr_sfc', 'dR_ObsTotal_clr_toa');
    clear dCRF dCRF_sfc dCRF_toa
        
    %%%%% step3: cal The effect of Ts on the atmosphere(dR_tstoa-dR_tssfc) (note effect on atoms is positive, but effect on earth is nagetive)
    % and main Effect (ta+alb+q+clod effect)
    %( 1sfc cld,2sfc clr,3toa cld,4toa clr)
    dR_tsAtom_cld = dR_ts(:, :, :, 3) - dR_ts(:, :, :, 1);
    dR_tsAtom_clr = dR_ts(:, :, :, 4) - dR_ts(:, :, :, 2);
    save([radEfectPath, 'dR_tsAtom_cld.mat'], 'dR_tsAtom_cld');
    save([radEfectPath, 'dR_tsAtom_clr.mat'], 'dR_tsAtom_clr');
    dR_mainEffect_sfc = dR_cloud_sfc + dR_main(:, :, :, 1);
    dR_mainEffect_toa = dR_cloud_toa + dR_main(:, :, :, 3);
    save([radEfectPath, 'dR_mainEffect_sfc.mat'], 'dR_mainEffect_sfc');
    save([radEfectPath, 'dR_mainEffect_toa.mat'], 'dR_mainEffect_toa');

end
t = toc; disp(t)