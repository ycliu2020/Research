%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-09 15:52:00
% LastEditTime : 2020-07-14 09:07:01
% LastEditors  : LYC
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
nowpath = pwd;
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
load('/data1/liuyincheng/cmip6-process/z_globalVar/ERF_rec.mat')% 'ERF_rec', 'timeERF_rec'
lon_k = 0:2.5:357.5; nlonk = length(lon_k); % kernel lat lon
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment
for p_1 = 3:4%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    %CAMS-CSM1-0 didn't have sfc clear sky radiation, delete it
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(p_1);
    % exmPath
    exmPath = ['/data1/liuyincheng/cmip6-process/', level.time1{p_1}]; %/data1/liuyincheng/cmip6-process/2000-2014/
    %%% load ERF data
    if strcmp(Experiment{p_1}(1:3), 'ami')
        timeERF = timeERF_rec.hist(find(timeERF_rec.hist == tLin.startYear{p_1}):find(timeERF_rec.hist == tLin.endYear{p_1}));
        ERForcing = ERF_rec.hist(find(timeERF_rec.hist == tLin.startYear{p_1}):find(timeERF_rec.hist == tLin.endYear{p_1}));
    elseif strcmp(Experiment{p_1}(1:6), 'ssp245')
        timeERF = timeERF_rec.ssp(find(timeERF_rec.ssp == tLin.startYear{p_1}):find(timeERF_rec.ssp == tLin.endYear{p_1}));
        ERForcing = ERF_rec.ssp245(find(timeERF_rec.ssp == tLin.startYear{p_1}):find(timeERF_rec.ssp == tLin.endYear{p_1}));
    elseif strcmp(Experiment{p_1}(1:6), 'ssp370')
        timeERF = timeERF_rec.ssp(find(timeERF_rec.ssp == tLin.startYear{p_1}):find(timeERF_rec.ssp == tLin.endYear{p_1}));
        ERForcing = ERF_rec.ssp370(find(timeERF_rec.ssp == tLin.startYear{p_1}):find(timeERF_rec.ssp == tLin.endYear{p_1}));
    else
        ERForcing = 0;
        disp('this experient doesnt need ERF or havent input ERF Data!')
    end

    if length(ERForcing) ~= 1% if ERForcing=0 mean .eg. abrupt experint
        ERForcing = repmat(ERForcing, [12 1]);
        ERForcing = reshape(ERForcing, 1, length(ERForcing(:))); % yearx12, assume every month have the same radiative forcing
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model
    for level1 = 1:length(level.model2)
        % model path
        mdlPath = fullfile(exmPath, level.model2{level1});
        eval(['cd ', mdlPath]);
        disp(' ')
        disp([level.model2{level1}, ' model start!'])

        % ensemble member path
        esmName = getPath_fileName(mdlPath, '.');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensemble member
        for esmNum = 1:length(esmName)
            esmPath = fullfile(mdlPath, esmName{esmNum, 1});
            % path
            dvarsPath = fullfile(esmPath, level.process3{2}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/essember/anomaly
            varsPath = fullfile(esmPath, level.process3{1}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/essember/rawdata
            radEfectPath = fullfile(esmPath, level.process3{6}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/essember/radEffect/
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
            startMonth=1;
            dhus2 = zeros(144, 73, 24, ntime);
            for monNum = 1:ntime
                dhus2(:, :, :, monNum) = (log(hus(:, :, :, monNum)) - log(clim_hus(:, :, :, mod(monNum + startMonth- 2, 12) + 1))) .* clim_ta(:, :, :, mod(monNum + startMonth- 2, 12) + 1).^2 * Rv / Lv;
            end
            dhus2(isnan(dhus2)) = 0;
            save([dvarsPath, 'dhus2.mat'], 'dhus2');
            % loop kernels in differnt conditions( 1sfc cld,2sfc clr,3toa cld,4toa clr)
            kernelPath = fullfile(esmPath,level.process3{5}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/kernelsCal
            dR_hus = zeros(144, 72, ntime, 4);
            dR_alb = dR_hus; dR_ts = dR_hus; dR_ta = dR_hus; dR_total = dR_hus; dR_main = dR_hus;
            % dR_tas=zeros(144,72,ntime,2);dR_taOnly=dR_tas;

            for ii = 1:4%( 1sfc cld,2sfc clr,3toa cld,4toa clr)
                load([kernelPath, kernelsName{ii}])% read kernel
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
                    % taOnlyEffect2 = zeros(144, 73, plevel-2, ntime);
                    % t0Effect = zeros(144, 73, ntime);% mean sfc ta radeffect
                else
                    taEffect = zeros(144, 73, plevel, ntime);
                    % consider 2 levels near sfc
                    tasEffect2 = taEffect;
                    % taOnlyEffect2 = zeros(144, 73, plevel - 2, ntime);
                end

                % Kernels and mat data are both from bottom to top
                % mod function used to process the data starts from 1 and should be consist with kernel(start from 1)
                startMonth=1;
                for monNum = 1:ntime
                    wvlwEffect(:, :, :, monNum) = wv_lwkernel(:, :, :, mod(monNum + startMonth- 2, 12) + 1) .* dhus2(:, :, :, monNum);
                    wvswEffect(:, :, :, monNum) = wv_swkernel(:, :, :, mod(monNum + startMonth- 2, 12) + 1) .* dhus2(:, :, :, monNum);
                    tsEffect(:, :, monNum) = ts_lwkernel(:, :, mod(monNum + startMonth- 2, 12) + 1) .* dts(:, :, monNum);
                    albEffect(:, :, monNum) = alb_swkernel(:, :, mod(monNum + startMonth- 2, 12) + 1) .* dalb(:, :, monNum) * 100; % note that alb kernerl unite: W/m2/0.01

                    if ii <= 2% sfc
                        taEffect(:, :, 2:end, monNum) = t_lwkernel(:, :, 2:end, mod(monNum + startMonth- 2, 12) + 1) .* dta(:, :, :, monNum);
                        taEffect(:, :, 1, monNum) = squeeze(t_lwkernel(:, :, 1, mod(monNum + startMonth- 2, 12) + 1)) .* dts(:, :, monNum);
                        tasEffect(:, :, monNum) = taEffect(:, :, 1, monNum);
                        taOnlyEffect(:, :, :, monNum) = taEffect(:, :, 2:end, monNum);
                        % consider 2 levels near sfc(have problem)
                        tasEffect2(:, :, 2:end, monNum) = t_level2_lwkernel(:, :, 2:end, mod(monNum + startMonth- 2, 12) + 1) .* dta(:, :, :, monNum);
                        tasEffect2(:, :, 1, monNum) = squeeze(t_level2_lwkernel(:, :, 1, mod(monNum + startMonth- 2, 12) + 1)) .* dts(:, :, monNum);
                        else% toa
                        taEffect(:, :, :, monNum) = t_lwkernel(:, :, :, mod(monNum + startMonth- 2, 12) + 1) .* dta(:, :, :, monNum);
                        % consider 2 levels near sfc
                        tasEffect2(:, :, :, monNum) = t_level2_lwkernel(:, :, :, mod(monNum + startMonth- 2, 12) + 1) .* dta(:, :, :, monNum);
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

                %regrid 144x72(unite grids)
                wvlwEffect = autoRegrid3(lon_k, lat_k, time.date, wvlwEffect, lon_f, lat_f, time.date);
                wvswEffect = autoRegrid3(lon_k, lat_k, time.date, wvswEffect, lon_f, lat_f, time.date);
                taEffect = autoRegrid3(lon_k, lat_k, time.date, taEffect, lon_f, lat_f, time.date);
                tsEffect = autoRegrid3(lon_k, lat_k, time.date, tsEffect, lon_f, lat_f, time.date);
                albEffect = autoRegrid3(lon_k, lat_k, time.date, albEffect, lon_f, lat_f, time.date);
                % integration(taEffect,tsEffect,albEffect,totalEffect,husEffect)
                husEffect = wvlwEffect + wvswEffect;
                totalEffect = husEffect + taEffect + tsEffect + albEffect;
                mainEffect = husEffect + taEffect + albEffect; % prepare for cal ta+alb+q+clod effect

                if ii <= 2
                    taOnlyEffect = autoRegrid3(lon_k, lat_k, time.date, taOnlyEffect, lon_f, lat_f, time.date);
                    taOnlyEffect2 = autoRegrid3(lon_k, lat_k, time.date, taOnlyEffect2, lon_f, lat_f, time.date);
                    tasEffect = autoRegrid3(lon_k, lat_k, time.date, tasEffect, lon_f, lat_f, time.date);
                    tasEffect2 = autoRegrid3(lon_k, lat_k, time.date, tasEffect2, lon_f, lat_f, time.date);
                    % save the radEffect: dradEfect_sfc_cld.mat
                    save([radEfectPath, 'global_vars.mat'], 'lon_f', 'lat_f', 'time', 'plev_k', 'readme', 'timeseries', 'modelname')
                    save([radEfectPath, saveradEfectName{ii}], 'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', ...
                        'taEffect', 'taOnlyEffect', 'tasEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect', 'mainEffect');
                else
                    taOnlyEffect2 = autoRegrid3(lon_k, lat_k, time.date, taOnlyEffect2, lon_f, lat_f, time.date);
                    tasEffect2 = autoRegrid3(lon_k, lat_k, time.date, tasEffect2, lon_f, lat_f, time.date);
                    % save the radEffect: dradEfect_sfc_cld.mat
                    save([radEfectPath, 'global_vars.mat'], 'lon_f', 'lat_f', 'time', 'plev_k', 'readme', 'timeseries', 'modelname')
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
            load([dvarsPath, 'drlus.mat'])% surface_upwelling_longwave_flux_in_air
            load([dvarsPath, 'drsus.mat']); load([dvarsPath, 'drsuscs.mat'])% surface_upwelling_shortwave_flux_in_air
            load([dvarsPath, 'drsds.mat']); load([dvarsPath, 'drsdscs.mat'])% surface_downwelling_shortwave_flux_in_air
            load([dvarsPath, 'drlds.mat']); load([dvarsPath, 'drldscs.mat'])% surface_downwelling_longwave_flux_in_air
            % 2.toa dvars(includ cld and clr sky)
            load([dvarsPath, 'drsdt.mat'])% toa_incoming_shortwave_flux
            load([dvarsPath, 'drsut.mat']); load([dvarsPath, 'drsutcs.mat'])% toa_outgoing_shortwave_flux
            load([dvarsPath, 'drlut.mat']); load([dvarsPath, 'drlutcs.mat'])% toa_outgoing_longwave_flux

            % model output rad( 1sfc cld,2sfc clr,3toa cld,4toa clr)
            ll_rad(:, :, :, 1) = drlds - drlus; % Sfc_cld net thermal radiation,
            ll_rad(:, :, :, 2) = drldscs - drlus; % Sfc_clr net thermal radiation,
            ll_rad(:, :, :, 3) = -drlut; % toa_cld net thermal radiation
            ll_rad(:, :, :, 4) = -drlutcs; % toa_clr net thermal radiation

            ss_rad(:, :, :, 1) = drsds - drsus; % Sfc_cld net solar radiation,
            ss_rad(:, :, :, 2) = drsdscs - drsuscs; % Sfc_clr net solar radiation,
            ss_rad(:, :, :, 3) = drsdt - drsut; % toa_cld net solar radiation,
            ss_rad(:, :, :, 4) = drsdt - drsutcs; % toa_clr net solar radiation,
            ll_rad(isnan(ll_rad)) = 0;
            ss_rad(isnan(ss_rad)) = 0;

            %regrid 144x72(unite grids)
            for ii = 1:4
                l_rad(:, :, :, ii) = autoRegrid3(lon_k, lat_k, time.date, squeeze(ll_rad(:, :, :, ii)), lon_f, lat_f, time.date);
                s_rad(:, :, :, ii) = autoRegrid3(lon_k, lat_k, time.date, squeeze(ss_rad(:, :, :, ii)), lon_f, lat_f, time.date);
            end

            % 1.sfc, 2.toa( 1sfc cld,2sfc clr,3toa cld,4toa clr)
            dR_allsky = l_rad(:, :, :, [1 3]) + s_rad(:, :, :, [1 3]); % real net rad
            dR_clr = l_rad(:, :, :, [2 4]) + s_rad(:, :, :, [2 4]);
            % save real radEffect: eg:real_dradEfect.mat
            readme_realradEfect = {'final row: 1sfc cld,2sfc clr,3toa cld,4toa clr', 'dR_allsky,dR_clr:1.sfc, 2.toa'};
            save([radEfectPath, 'real_dradEfect.mat'], 'l_rad', 's_rad', 'dR_allsky', 'dR_clr', 'readme_realradEfect');

            %% cal the dRcloud and dR_co2(144,72,time,2)(final row means contain sfc and toa)
            % note that residual includ( 1sfc cld,2sfc clr,3toa cld,4toa clr)
            % dCRF=dR_allsky-dR_clr;%(1sfc 2toa)
            %method 1
            dvarsFeedback = dR_total(:, :, :, [2 4]) - dR_total(:, :, :, [1 3]);
            dCRF = dR_allsky - dR_clr; % cloud radiative forcing
            % dCRF = dR_clr - dR_allsky; % cloud radiative forcing
            dR_cloud = dCRF + dvarsFeedback;
            % add ERF (G0-G)/G~0.16 according soden paper
            if length(ERForcing) ~= 1

                for timeLine = 1:ntime
                    dR_cloud(:, :, timeLine, 2) = dR_cloud(:, :, timeLine, 2) + ERForcing(timeLine) .* 0.16;
                end

            end

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
            save([radEfectPath, 'dR_residual_cld_sfc.mat'], 'dR_residual_cld_sfc');
            save([radEfectPath, 'dR_residual_cld_toa.mat'], 'dR_residual_cld_toa');
            save([radEfectPath, 'dR_residual_clr_sfc.mat'], 'dR_residual_clr_sfc');
            save([radEfectPath, 'dR_residual_clr_toa.mat'], 'dR_residual_clr_toa');
            save([radEfectPath, 'dR_ObsTotal.mat'], 'dR_ObsTotal_cld_sfc', 'dR_ObsTotal_cld_toa', ...
                'dR_ObsTotal_clr_sfc', 'dR_ObsTotal_clr_toa');

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
            disp([esmName{esmNum, 1}, ' ensemble is done!'])
        end

        disp([level.model2{level1}, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{p_1}, ' era is done!'])
    clear dCRF_sfc dR_cloud_sfc dR_residual_cld_sfc dR_residual_clr_sfc dR_residual_cld_toa dR_residual_clr_toa dR_ObsTotal_cld_sfc dR_ObsTotal_clr_sfc dR_ObsTotal_cld_toa dR_ObsTotal_clr_toa ss_rad ll_rad l_rad s_rad
end

t = toc; disp(t)
