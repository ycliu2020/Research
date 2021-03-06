%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-07-06 08:55:44
% LastEditTime : 2021-04-13 20:22:49
% LastEditors  : Please set LastEditors
% Description  : compare the hus computation of Angeline G and now applied
% FilePath     : /code/p1_processObserveData/ERA5/checkHusAndRcloud/ERA5_s3_calRadEff_dp_v2_compareHus.m
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
saveradEffectName = cell(1, 4);
n = 0;

for i = 1:2

    for j = 1:2
        n = n + 1;
        kernelsName{n} = ['kernels_', kernelsHightLabel{i}, '_', kernelsLevel2_sky{j}, '.mat'];
        saveradEffectName{n} = ['dradEffect_', kernelsHightLabel{i}, '_', kernelsLevel2_sky{j}, '.mat'];
    end

end

saveTrend_radEffectName = {'trend_dradEffect_sfc_cld.mat', 'trend_dradEffect_toa_cld.mat'};
plevel = 24; % wv_lwkernel and wv_swkernel level;
t_scflevel = 25; % sfc t_lwkernel level;

[readme, level, tLin, vars] = obsParameters('ERA5');

for p_1 = 4:4
    varsPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{1}); %rawdata
    dvarsPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{2}); %anomaly
    kernelCalPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{4}); % kernelCal
    radEffectPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{5}); %radEffect
    auto_mkdir(radEffectPath)
    % load dvars
    load([dvarsPath, 'global_vars.mat']) % lat_k lon_k time plev_k
    load([dvarsPath, 'meto_dvars.mat']) % dhus and clim_hus...
    load([varsPath, 'meto_vars.mat'], 'hus') %
    dta(isnan(dta)) = 0;
    dalb(isnan(dalb)) = 0; dts(isnan(dts)) = 0;
    ntime = length(time);
    lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f);
    lon_f = lon_k; nlonf = length(lon_f);
    startMonth = tLin.startMonth{p_1};

    %%%%% step1: cal dReffect
    % cal the fixed moisture because kernel q's unit is W/m2/K/100mb
    dhus2 = zeros(144, 73, 24, ntime);

    %% Method1:
    % dln(q)/dT = Lv/(Rv*T**2)->dT=dln(q)*T^2*Rv/Lv
    logHus = log(hus);
    [dlogHus, logHus_clim] = monthlyAnomaly(144, 73, plevel, time, logHus, startMonth); %[nlongitude,nlatitude,time,var,startmonth]

    for monNum = 1:ntime
        dhus2(:, :, :, monNum) = dlogHus(:, :, :, monNum) .* clim_ta(:, :, :, mod(monNum + startMonth - 2, 12) + 1).^2 * Rv / Lv;
    end

    dhus2(isnan(dhus2)) = 0;
    dhus2Infor = 'Unite: K, according to dln(q)/dT = Lv/(Rv*T**2) equation';

    %% Method2: Angeline G method
    p=zeros(144,73,24);
    for i = 1:144
        for j = 1:73
            p(i,j,:)=plev_k;
        end
    end
    p=repmat(p,1,1,1,12);
    qs1 = calcSatSpecHum(clim_ta, p); %g/kg p(lon, lat, plev, time) hpa
    rh = 1000 * clim_hus ./ qs1; %

    dlogqdt=zeros(144,73,24,ntime);
    for yearNum = 1:ntime/12
        qs2 = calcSatSpecHum(clim_ta + dta(:,:,:,(yearNum-1)*12+1:yearNum*12), p); % g/kg
        dqsdt = (qs2 - qs1) ./ dta(:,:,:,(yearNum-1)*12+1:yearNum*12); % g/kg/K
        dqdt = rh .* dqsdt;
        dlogqdt(:,:,:,(yearNum-1)*12+1:yearNum*12) = dqdt ./ (1000 * clim_hus);
    end
    
    %% Method3: Yi huang method
    dhus3=zeros(144,73,24,ntime);
    dws=dlnws(clim_ta);
    for yearNum = 1:ntime/12
        dhus3(:,:,:,(yearNum-1)*12+1:yearNum*12)=dhus(:,:,:,(yearNum-1)*12+1:yearNum*12).*dws.*clim_ta.^2 * Rv / Lv;
    end
    
    
    % save([dvarsPath, 'dhus2.mat'], 'dhus2', 'dhus2Infor');

    % loop kernels in differnt conditions(1sfc cld,2sfc clr,3toa cld,4toa clr)
    dR_hus = zeros(144, 72, ntime, 4);
    dR_alb = dR_hus; dR_ts = dR_hus; dR_ta = dR_hus; dR_nonCloud = dR_hus; dR_nonCloudAndTs = dR_hus;
    % dR_tas=zeros(144,72,ntime,2);dR_taOnly=dR_tas;

    for ii = 1:4 %( 1sfc cld,2sfc clr,3toa cld,4toa clr)
        load([kernelCalPath, kernelsName{ii}]) % read kernel
        % prepare zero mixre
        wvlwEffect = zeros(144, 73, plevel, ntime);
        wvswEffect = zeros(144, 73, plevel, ntime);
        tsEffect = zeros(144, 73, ntime);
        albEffect = zeros(144, 73, ntime);

        if ii <= 2 % sfc
            taEffect = zeros(144, 73, t_scflevel, ntime);
            % consider 1 levels near sfc
            tasEffect = zeros(144, 73, ntime);
            taOnlyEffect = zeros(144, 73, plevel, ntime);
            % consider 2 levels near sfc
            tasEffect2 = taEffect;
            tasEffect1 = taEffect;
        else
            taEffect = zeros(144, 73, plevel, ntime);
            % consider 2 levels near sfc
            tasEffect2 = taEffect;
            tasEffect1 = taEffect;
        end

        % Kernels and mat data are both from bottom to top
        % mod function used to process the data starts from 1 and should be consist with kernel(start from 1)

        for monNum = 1:ntime
            wvlwEffect(:, :, :, monNum) = wv_lwkernel(:, :, :, mod(monNum + startMonth - 2, 12) + 1) .* dhus2(:, :, :, monNum);
            wvswEffect(:, :, :, monNum) = wv_swkernel(:, :, :, mod(monNum + startMonth - 2, 12) + 1) .* dhus2(:, :, :, monNum);
            tsEffect(:, :, monNum) = ts_lwkernel(:, :, mod(monNum + startMonth - 2, 12) + 1) .* dts(:, :, monNum);
            albEffect(:, :, monNum) = alb_swkernel(:, :, mod(monNum + startMonth - 2, 12) + 1) .* dalb(:, :, monNum) * 100; % note that alb kernerl unite: W/m2/0.01

            if ii <= 2 % sfc
                taEffect(:, :, 2:end, monNum) = t_lwkernel(:, :, 2:end, mod(monNum + startMonth - 2, 12) + 1) .* dta(:, :, :, monNum);
                taEffect(:, :, 1, monNum) = squeeze(t_lwkernel(:, :, 1, mod(monNum + startMonth - 2, 12) + 1)) .* dts(:, :, monNum);
                tasEffect(:, :, monNum) = taEffect(:, :, 1, monNum);
                taOnlyEffect(:, :, :, monNum) = taEffect(:, :, 2:end, monNum);

                % consider 1 levels near sfc
                tasEffect1(:, :, 2:end, monNum) = t_level1_lwkernel(:, :, 2:end, mod(monNum + startMonth - 2, 12) + 1) .* dta(:, :, :, monNum);
                tasEffect1(:, :, 1, monNum) = squeeze(t_level1_lwkernel(:, :, 1, mod(monNum + startMonth - 2, 12) + 1)) .* dts(:, :, monNum);

                % consider 2 levels near sfc(have problem)
                tasEffect2(:, :, 2:end, monNum) = t_level2_lwkernel(:, :, 2:end, mod(monNum + startMonth - 2, 12) + 1) .* dta(:, :, :, monNum);
                tasEffect2(:, :, 1, monNum) = squeeze(t_level2_lwkernel(:, :, 1, mod(monNum + startMonth - 2, 12) + 1)) .* dts(:, :, monNum);
                else % toa
                taEffect(:, :, :, monNum) = t_lwkernel(:, :, :, mod(monNum + startMonth - 2, 12) + 1) .* dta(:, :, :, monNum);

                % consider 1 levels near sfc
                tasEffect1(:, :, :, monNum) = t_level1_lwkernel(:, :, :, mod(monNum + startMonth - 2, 12) + 1) .* dta(:, :, :, monNum);
                % consider 2 levels near sfc
                tasEffect2(:, :, :, monNum) = t_level2_lwkernel(:, :, :, mod(monNum + startMonth - 2, 12) + 1) .* dta(:, :, :, monNum);
            end

        end

        wvlwEffect = squeeze(nansum(wvlwEffect(:, :, :, :), 3));
        wvswEffect = squeeze(nansum(wvswEffect(:, :, :, :), 3));
        taEffect = squeeze(nansum(taEffect(:, :, :, :), 3));
        tasEffect2 = squeeze(nansum(tasEffect2(:, :, :, :), 3));
        tasEffect1 = squeeze(nansum(tasEffect1(:, :, :, :), 3));

        if ii <= 2
            taOnlyEffect = squeeze(nansum(taOnlyEffect(:, :, :, :), 3));
        end

        % regrid 144x72(unite grids)
        wvlwEffect = autoRegrid3(lon_k, lat_k, time, wvlwEffect, lon_f, lat_f, time);
        wvswEffect = autoRegrid3(lon_k, lat_k, time, wvswEffect, lon_f, lat_f, time);
        taEffect = autoRegrid3(lon_k, lat_k, time, taEffect, lon_f, lat_f, time);
        tsEffect = autoRegrid3(lon_k, lat_k, time, tsEffect, lon_f, lat_f, time);
        albEffect = autoRegrid3(lon_k, lat_k, time, albEffect, lon_f, lat_f, time);

        % integration(taEffect,tsEffect,albEffect,nonCloudEffect,husEffect)
        husEffect = wvlwEffect + wvswEffect;
        nonCloudEffect = husEffect + taEffect + tsEffect + albEffect;
        nonCloudAndTsEffect = husEffect + taEffect + albEffect; % prepare for cal ta+alb+q+clod effect

        if ii <= 2
            taOnlyEffect = autoRegrid3(lon_k, lat_k, time, taOnlyEffect, lon_f, lat_f, time);
            tasEffect = autoRegrid3(lon_k, lat_k, time, tasEffect, lon_f, lat_f, time);
            tasEffect2 = autoRegrid3(lon_k, lat_k, time, tasEffect2, lon_f, lat_f, time);
            tasEffect1 = autoRegrid3(lon_k, lat_k, time, tasEffect1, lon_f, lat_f, time);

            % save the radEffect: dradEffect_sfc_cld.mat
            % save([radEffectPath, 'global_vars.mat'], 'lon_f', 'lat_f', 'lon_k', 'lat_k', 'plev_k', 'time')
            % save([radEffectPath, saveradEffectName{ii}], 'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', ...
                %     'taEffect', 'taOnlyEffect', 'tasEffect', 'tasEffect2','tasEffect1', 'nonCloudEffect', 'nonCloudAndTsEffect');
        else
            tasEffect2 = autoRegrid3(lon_k, lat_k, time, tasEffect2, lon_f, lat_f, time);
            tasEffect1 = autoRegrid3(lon_k, lat_k, time, tasEffect1, lon_f, lat_f, time);
            % save the radEffect: dradEffect_sfc_cld.mat
            % save([radEffectPath, 'global_vars.mat'], 'lon_f', 'lat_f', 'lon_k', 'lat_k', 'plev_k', 'time')
            % save([radEffectPath, saveradEffectName{ii}], 'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', ...
                %     'taEffect', 'tasEffect2', 'tasEffect1', 'nonCloudEffect', 'nonCloudAndTsEffect');
        end

        % store vars in one var( 1sfc cld,2sfc clr,3toa cld,4toa clr)
        dR_hus(:, :, :, ii) = husEffect; % q
        dR_alb(:, :, :, ii) = albEffect; %albedo
        dR_ts(:, :, :, ii) = tsEffect; % ts
        dR_ta(:, :, :, ii) = taEffect; % t
        dR_nonCloud(:, :, :, ii) = nonCloudEffect; % total
        dR_nonCloudAndTs(:, :, :, ii) = nonCloudAndTsEffect; % prepare for cal ta+alb+q+clod effect
    end

    readme_radEffect = 'final column: 1sfc cld,2sfc clr,3toa cld,4toa clr, note that dR_tas and dR_taOnly shouldnt have value in toa';
    % save([radEffectPath, 'dradEffect_union.mat'], 'dR_hus', 'dR_alb', 'dR_ts', 'dR_ta', 'dR_nonCloud', 'readme_radEffect');

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
    dR_net_cld = l_rad(:, :, :, [1 3]) + s_rad(:, :, :, [1 3]); % real net rad
    dR_net_clr = l_rad(:, :, :, [2 4]) + s_rad(:, :, :, [2 4]);
    % save real radEffect: eg:real_dradEffect.mat
    readme_realradEffect = {'final row: 1sfc cld,2sfc clr,3toa cld,4toa clr', 'dR_net_cld,dR_net_clr:1.sfc, 2.toa'};
    % save([radEffectPath, 'real_dradEffect.mat'], 'l_rad', 's_rad', 'dR_net_cld', 'dR_net_clr', 'readme_realradEffect');
    clear ll_rad ss_rad l_rad s_rad
    %% cal the dRcloud and dR_co2(144,72,time,2)(final row means contain sfc and toa)
    % note that residual includ( 1sfc cld,2sfc clr,3toa cld,4toa clr)
    % dCRF=dR_net_cld-dR_net_clr;%(1sfc 2toa)
    %method 1
    dvarsFeedback = dR_nonCloud(:, :, :, [2 4]) - dR_nonCloud(:, :, :, [1 3]);
    dCRF = dR_net_cld - dR_net_clr; % cloud radiative forcing
    % dCRF = dR_net_clr - dR_net_cld; % cloud radiative forcing
    dR_cloud = dCRF + dvarsFeedback;
    
    dR_residual_clr = dR_net_clr - dR_nonCloud(:, :, :, [2 4]);
    dR_cloud = dR_net_cld - dR_nonCloud(:,:,:,[1 3])- dR_residual_clr=dR_net_cld - dR_nonCloud(:,:,:,[1 3])-(dR_net_clr - dR_nonCloud(:, :, :, [2 4]))

    if strcmp(cldClr_char,'clr')
        dR_res = dR - dR_noncloud;
        dR_lw_res = dR_lw - dR_lw_noncloud;
        dR_sw_res = dR_sw - dR_sw_noncloud;
    elseif strcmp(cldClr_char,'all')
        dR_c = dR - dR_noncloud - dR_res;            
        dR_lw_c = dR_lw - dR_lw_noncloud - dR_lw_res;
        dR_sw_c = dR_sw - dR_sw_noncloud - dR_sw_res;
    end 

    
    %method 2
    dR_residual_clr = dR_net_clr - dR_nonCloud(:, :, :, [2 4]);
    dR_residual_cld = dR_net_cld - dR_nonCloud(:, :, :, [1 3]) - dR_cloud; % r terms, indicate co2/aerosol forcing
    dR_residual_cld1 = dR_net_cld - dR_nonCloud(:, :, :, [1 3]); % r terms, indicate cloud feedback and co2/aerosol forcing
    % dR_cloud = dR_residual_cld1 - dR_residual_clr;

    dR_cloud_sfc = dR_cloud(:, :, :, 1); dR_cloud_toa = dR_cloud(:, :, :, 2);
    dCRF_sfc = dCRF(:, :, :, 1); dCRF_toa = dCRF(:, :, :, 2);
    dR_residual_cld_sfc = dR_residual_cld(:, :, :, 1); dR_residual_cld_toa = dR_residual_cld(:, :, :, 2);
    dR_residual_clr_sfc = dR_residual_clr(:, :, :, 1); dR_residual_clr_toa = dR_residual_clr(:, :, :, 2);
    % equal to real_dradEffect.mat
    dR_net_cld_sfc = dR_net_cld(:, :, :, 1); dR_net_cld_toa = dR_net_cld(:, :, :, 2);
    dR_net_clr_sfc = dR_net_clr(:, :, :, 1); dR_net_clr_toa = dR_net_clr(:, :, :, 2);
    % % save dR_cloud and dR_residual_cld and observe total effect
    % save([radEffectPath, 'dCRF.mat'], 'dCRF');
    % save([radEffectPath, 'dCRF_sfc.mat'], 'dCRF_sfc');
    % save([radEffectPath, 'dCRF_toa.mat'], 'dCRF_toa');
    % save([radEffectPath, 'dR_cloud.mat'], 'dR_cloud_sfc', 'dR_cloud_toa');
    % save([radEffectPath, 'dR_residual_cld_sfc.mat'], 'dR_residual_cld_sfc');
    % save([radEffectPath, 'dR_residual_cld_toa.mat'], 'dR_residual_cld_toa');
    % save([radEffectPath, 'dR_residual_clr_sfc.mat'], 'dR_residual_clr_sfc');
    % save([radEffectPath, 'dR_residual_clr_toa.mat'], 'dR_residual_clr_toa');
    % save([radEffectPath, 'dR_net.mat'], 'dR_net_cld_sfc', 'dR_net_cld_toa', ...
        %     'dR_net_clr_sfc', 'dR_net_clr_toa');
    % clear dCRF dCRF_sfc dCRF_toa

    %%%%% step3: cal The effect of Ts on the atmosphere(dR_tstoa-dR_tssfc) (note effect on atoms is positive, but effect on earth is nagetive)
    % and main Effect (ta+alb+q+clod effect)
    %( 1sfc cld,2sfc clr,3toa cld,4toa clr)3
    
    dR_tsAtom_cld = dR_ts(:, :, :, 3) - dR_ts(:, :, :, 1);
    dR_tsAtom_clr = dR_ts(:, :, :, 4) - dR_ts(:, :, :, 2);
    % save([radEffectPath, 'dR_tsAtom_cld.mat'], 'dR_tsAtom_cld');
    % save([radEffectPath, 'dR_tsAtom_clr.mat'], 'dR_tsAtom_clr');
    dR_nonTs_sfc = dR_cloud_sfc + dR_nonCloudAndTs(:, :, :, 1);
    dR_nonTs_toa = dR_cloud_toa + dR_nonCloudAndTs(:, :, :, 3);
    % save([radEffectPath, 'dR_nonTs_sfc.mat'], 'dR_nonTs_sfc');
    % save([radEffectPath, 'dR_nonTs_toa.mat'], 'dR_nonTs_toa');

end

t = toc; disp(t)

function dws=dlnws(T)
    
    pliq=0.01.*exp(54.842763- 6763.22./T-4.21.*log(T) + 0.000367.*T+tanh(0.0415.*(T-218.8)).*...
                    (53.878 -1331.22/T-9.44523*log(T) + 0.014025*T));
    pice=exp(9.550426-5723.265/T+3.53068*log(T)-0.00728332*T)/100.;
    T1=T+1.;
    pliq1=0.01*exp(54.842763- 6763.22 / T1-4.21*log(T1) + 0.000367*T1+tanh(0.0415*(T1 - 218.8)).*...
                    (53.878 -1331.22/T1-9.44523*log(T1) + 0.014025*T1));
    pice1=exp(9.550426-5723.265/T1+3.53068*log(T1)-0.00728332*T1)/100.;
    
    ws=T;
    ws1=ws;
    ws(T>=273)=pliq(T>=273);
    ws(T<273)=pice(T<273);

    ws1(T>=273)=pliq1(T>=273);
    ws1(T<273)=pice1(T<273);

    dws=ws./(ws1-ws);

end
    