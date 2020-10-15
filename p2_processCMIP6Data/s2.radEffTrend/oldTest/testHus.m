%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% mothly data
% cal the 1.dradEfect 2.cloudradEffect 3.their trend,4.year trend mask and  yr_cc
% this script mainly include:
% PS: m mean month, s mean season, yr mean year
%
% hope to use an unite scrip to read and process all the varies we need
% Compare the pre rad and ts data from 2000-03 to 2018-02(18 years)
%
% time:2000.01-2014.12(interval:15*12);1980.01-2014.12(interval:35*12); 2015.01-2099.12(interval:85*12)
% initial time in hist(1740 total): 1,561 of 1740(2000.03);1,321 of 1740(1980.01)
% initial time in futrue(1032 total): 1 of 1032(2015.01);
%
% cal mainly include dts and dRHS trend
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear; clc; tic;

% model settings
for p_1 = 1:1%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    nowpath = pwd;
    % constant
    Rv = 487.5;
    Lv = 2.5e6;
    %CAMS-CSM1-0 didn't have sfc clear sky radiation, delete it
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(p_1);

    % inputPath
    inputPath = ['/data1/liuyincheng/cmip6-process/', level.time1{p_1}]; %/data1/liuyincheng/cmip6-process/2000-2014/
    % kernels in differnt conditions
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
    % % load maskword file
    % load mask_cp144.mat

    % model loop
    for level1 = 9:9%1:length(level.model2)
        % load dvars
        dvarsPath = [inputPath, level.model2{level1}, '/', level.process3{2}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
        varsPath = [inputPath, level.model2{level1}, '/', level.process3{1}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/rawdata
        load([varsPath, 'hus.mat'])
        load([dvarsPath, 'global_vars.mat'])% latf lonf time plevf readme
        load([dvarsPath, 'dhus.mat'])% dhus and clim_hus
        load([dvarsPath, 'dalb.mat'])%
        load([dvarsPath, 'dta.mat'])%
        load([dvarsPath, 'dts.mat'])%
        dta(isnan(dta)) = 0;
        dalb(isnan(dalb)) = 0; dts(isnan(dts)) = 0;

        % load real rad dvars
        % 1.sfc dvars(includ cld and clr sky)
        load([dvarsPath, 'drlus.mat'])% surface_upwelling_longwave_flux_in_air
        load([dvarsPath, 'drsus.mat']); load([dvarsPath, 'drsuscs.mat'])% surface_upwelling_shortwave_flux_in_air
        load([dvarsPath, 'drsds.mat']); load([dvarsPath, 'drsdscs.mat'])% surface_downwelling_shortwave_flux_in_air
        load([dvarsPath, 'drlds.mat']); load([dvarsPath, 'drldscs.mat'])% surface_downwelling_longwave_flux_in_air
        % 2.toa dvars(includ cld and clr sky)
        load([dvarsPath, 'drsdt.mat'])% toa_incoming_shortwave_flux
        load([dvarsPath, 'drsut.mat']); load([dvarsPath, 'drsutcs.mat'])% toa_outgoing_shortwave_flux
        load([dvarsPath, 'drlut.mat']); load([dvarsPath, 'drlutcs.mat'])% toa_outgoing_longwave_flux

        % cal the fixed moisture because kernel q's unit is W/m2/K/100mb
        dhus2 = zeros(144, 73, 24, length(time));

        for kk = 1:length(time)
            dhus2(:, :, :, kk) = (log(hus(:, :, :, kk)) - log(clim_hus(:, :, :, mod(kk - 1, 12) + 1))) .* clim_ta(:, :, :, mod(kk - 1, 12) + 1).^2 * Rv / Lv;
        end

        dhus2(isnan(dhus2)) = 0;
        % loop kernels in differnt conditions( 1sfc cld,2sfc clr,3toa cld,4toa clr)
        kernelPath = [inputPath, level.model2{level1}, '/', level.process3{5}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/kernelsCal
        dR_hus = zeros(144, 72, length(time), 4);
        dR_alb = dR_hus; dR_ts = dR_hus; dR_ta = dR_hus; dR_total = dR_hus; dR_main = dR_hus;
        % dR_tas=zeros(144,72,length(time),2);dR_taOnly=dR_tas;

        for ii = 1:4%( 1sfc cld,2sfc clr,3toa cld,4toa clr)
            load([kernelPath, kernelsName{ii}])% read kernel
            %% step1: cal dReffect
            % prepare zero mixre
            wvlwEffect = zeros(144, 73, plevel, length(time));
            wvswEffect = zeros(144, 73, plevel, length(time));
            if ii <= 2% sfc
                taEffect = zeros(144, 73, t_scflevel, length(time));
                % consider 1 levels near sfc
                tasEffect = zeros(144, 73, length(time));
                taOnlyEffect = zeros(144, 73, plevel, length(time));
                % consider 2 levels near sfc
                tasEffect2 = taEffect;
                % taOnlyEffect2 = zeros(144, 73, plevel-2, length(time));
                % t0Effect = zeros(144, 73, length(time));% mean sfc ta radeffect
            else
                taEffect = zeros(144, 73, plevel, length(time));
                % consider 2 levels near sfc
                tasEffect2 = taEffect;
                % taOnlyEffect2 = zeros(144, 73, plevel - 2, length(time));
            end

            tsEffect = zeros(144, 73, length(time));
            albEffect = zeros(144, 73, length(time));
            % Kernels and mat data are both from bottom to top
            % mod function used to process the data starts from 1 and should be consist with kernel(start from 1)

            for kk = 1:length(time)
                wvlwEffect(:, :, :, kk) = wv_lwkernel(:, :, :, mod(kk - 1, 12) + 1) .* dhus2(:, :, :, kk);
                wvswEffect(:, :, :, kk) = wv_swkernel(:, :, :, mod(kk - 1, 12) + 1) .* dhus2(:, :, :, kk);
                tsEffect(:, :, kk) = ts_lwkernel(:, :, mod(kk - 1, 12) + 1) .* dts(:, :, kk);
                albEffect(:, :, kk) = alb_swkernel(:, :, mod(kk - 1, 12) + 1) .* dalb(:, :, kk) * 100;

                if ii <= 2% sfc
                    taEffect(:, :, 2:end, kk) = t_lwkernel(:, :, 2:end, mod(kk - 1, 12) + 1) .* dta(:, :, :, kk);
                    taEffect(:, :, 1, kk) = squeeze(t_lwkernel(:, :, 1, mod(kk - 1, 12) + 1)) .* dts(:, :, kk);
                    tasEffect(:, :, kk) = taEffect(:, :, 1, kk);
                    taOnlyEffect(:, :, :, kk) = taEffect(:, :, 2:end, kk);
                    % consider 2 levels near sfc(have problem)
                    tasEffect2(:, :, 2:end, kk) = t_level2_lwkernel(:, :, 2:end, mod(kk - 1, 12) + 1) .* dta(:, :, :, kk);
                    tasEffect2(:, :, 1, kk) = squeeze(t_level2_lwkernel(:, :, 1, mod(kk - 1, 12) + 1)) .* dts(:, :, kk);
                    else% toa
                    taEffect(:, :, :, kk) = t_lwkernel(:, :, :, mod(kk - 1, 12) + 1) .* dta(:, :, :, kk);
                    % consider 2 levels near sfc
                    tasEffect2(:, :, :, kk) = t_level2_lwkernel(:, :, :, mod(kk - 1, 12) + 1) .* dta(:, :, :, kk);
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

            % %regrid 144x72(unite grids)
            % lat = 88.75:-2.5:-88.75; nlat = length(lat);
            % lon = lonf; nlon = length(lon);
            % nlonf = length(lonf); nlatf = length(latf);
            % [Xlon, Ylat, Ttime] = meshgrid(lat, lon, time);
            % [Xlonf, Ylatf, Ttimef] = meshgrid(latf, lonf, time);
            % wvlwEffect = interp3(Xlonf, Ylatf, Ttimef, wvlwEffect, Xlon, Ylat, Ttime);
            % wvswEffect = interp3(Xlonf, Ylatf, Ttimef, wvswEffect, Xlon, Ylat, Ttime);
            % taEffect = interp3(Xlonf, Ylatf, Ttimef, taEffect, Xlon, Ylat, Ttime);
            % tsEffect = interp3(Xlonf, Ylatf, Ttimef, tsEffect, Xlon, Ylat, Ttime);
            % albEffect = interp3(Xlonf, Ylatf, Ttimef, albEffect, Xlon, Ylat, Ttime);
            % % integration(taEffect,tsEffect,albEffect,totalEffect,husEffect)
            % husEffect = wvlwEffect + wvswEffect;
            % totalEffect = husEffect + taEffect + tsEffect + albEffect;
            % mainEffect = husEffect + taEffect + albEffect; % prepare for cal ta+alb+q+clod effect

            % if ii <= 2
            %     taOnlyEffect = interp3(Xlonf, Ylatf, Ttimef, taOnlyEffect, Xlon, Ylat, Ttime);
            %     taOnlyEffect2 = interp3(Xlonf, Ylatf, Ttimef, taOnlyEffect2, Xlon, Ylat, Ttime);
            %     tasEffect = interp3(Xlonf, Ylatf, Ttimef, tasEffect, Xlon, Ylat, Ttime);
            %     tasEffect2 = interp3(Xlonf, Ylatf, Ttimef, tasEffect2, Xlon, Ylat, Ttime);
            %     % save the radEffect: dradEfect_sfc_cld.mat
            %     radefectOutPath = [inputPath, level.model2{level1}, '/', level.process3{6}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
            %     auto_mkdir(radefectOutPath)
            %     save([radefectOutPath, 'global_vars.mat'], 'lon', 'lat', 'time', 'plevf', 'readme', 'timeseries', 'modelname')
            %     save([radefectOutPath, saveradEfectName{ii}], 'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', ...
            %         'taEffect', 'taOnlyEffect', 'tasEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect', 'mainEffect');
            % else
            %     taOnlyEffect2 = interp3(Xlonf, Ylatf, Ttimef, taOnlyEffect2, Xlon, Ylat, Ttime);
            %     tasEffect2 = interp3(Xlonf, Ylatf, Ttimef, tasEffect2, Xlon, Ylat, Ttime);
            %     % save the radEffect: dradEfect_sfc_cld.mat
            %     radefectOutPath = [inputPath, level.model2{level1}, '/', level.process3{6}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
            %     auto_mkdir(radefectOutPath)
            %     save([radefectOutPath, 'global_vars.mat'], 'lon', 'lat', 'time', 'plevf', 'readme', 'timeseries', 'modelname')
            %     save([radefectOutPath, saveradEfectName{ii}], 'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', ...
            %         'taEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect', 'mainEffect');
            % end

            % % store vars in one var( 1sfc cld,2sfc clr,3toa cld,4toa clr)
            % dR_hus(:, :, :, ii) = husEffect; % q
            % dR_alb(:, :, :, ii) = albEffect; %albedo
            % dR_ts(:, :, :, ii) = tsEffect; % ts
            % dR_ta(:, :, :, ii) = taEffect; % t
            % dR_total(:, :, :, ii) = totalEffect; % total
            % dR_main(:, :, :, ii) = mainEffect; % prepare for cal ta+alb+q+clod effect
            % readme_radEfect = 'final column: 1sfc cld,2sfc clr,3toa cld,4toa clr, note that dR_tas and dR_taOnly shouldnt have value in toa';
            % save([radefectOutPath, 'dradEfect_union.mat'], 'dR_hus', 'dR_alb', 'dR_ts', 'dR_ta', 'dR_total', 'readme_radEfect');
        end

        % %% step2: cal cloud effect(definen down is postive)
        % % real rad( 1sfc cld,2sfc clr,3toa cld,4toa clr)
        % ll_rad(:, :, :, 1) = drlds - drlus; % Sfc_cld net thermal radiation,
        % ll_rad(:, :, :, 2) = drldscs - drlus; % Sfc_clr net thermal radiation,
        % ll_rad(:, :, :, 3) = -drlut; % toa_cld net thermal radiation
        % ll_rad(:, :, :, 4) = -drlutcs; % toa_clr net thermal radiation

        % ss_rad(:, :, :, 1) = drsds - drsus; % Sfc_cld net solar radiation,
        % ss_rad(:, :, :, 2) = drsdscs - drsuscs; % Sfc_clr net solar radiation,
        % ss_rad(:, :, :, 3) = drsdt - drsut; % toa_cld net solar radiation,
        % ss_rad(:, :, :, 4) = drsdt - drsutcs; % toa_clr net solar radiation,
        % ll_rad(isnan(ll_rad)) = 0;
        % ss_rad(isnan(ss_rad)) = 0;

        % %regrid 144x72(unite grids)
        % lat = 88.75:-2.5:-88.75; nlat = length(lat);
        % lon = lonf; nlon = length(lon);
        % nlonf = length(lonf); nlatf = length(latf);
        % [Xlon, Ylat, Ttime] = meshgrid(lat, lon, time);
        % [Xlonf, Ylatf, Ttimef] = meshgrid(latf, lonf, time);

        % for ii = 1:4
        %     l_rad(:, :, :, ii) = interp3(Xlonf, Ylatf, Ttimef, squeeze(ll_rad(:, :, :, ii)), Xlon, Ylat, Ttime);
        %     s_rad(:, :, :, ii) = interp3(Xlonf, Ylatf, Ttimef, squeeze(ss_rad(:, :, :, ii)), Xlon, Ylat, Ttime);
        % end

        % % 1.sfc, 2.toa( 1sfc cld,2sfc clr,3toa cld,4toa clr)
        % dR_allsky = l_rad(:, :, :, [1 3]) + s_rad(:, :, :, [1 3]); % real net rad
        % dR_clr = l_rad(:, :, :, [2 4]) + s_rad(:, :, :, [2 4]);
        % % save real radEffect: eg:real_dradEfect.mat
        % readme_realradEfect = {'final row: 1sfc cld,2sfc clr,3toa cld,4toa clr', 'dR_allsky,dR_clr:1.sfc, 2.toa'};
        % save([radefectOutPath, 'real_dradEfect.mat'], 'l_rad', 's_rad', 'dR_allsky', 'dR_clr', 'readme_realradEfect');

        % %% cal the dRcloud and dR_co2(144,72,time,2)(final row means contain sfc and toa)
        % % note that residual includ( 1sfc cld,2sfc clr,3toa cld,4toa clr)
        % % dCRF=dR_clr-dR_allsky;%(1sfc 2toa)
        % % dR_comp=
        % %method 1
        % dKernel = dR_total(:, :, :, [2 4]) - dR_total(:, :, :, [1 3]);
        % dCRF = dR_clr - dR_allsky; % cloud radiative forcing
        % dR_cloud = dCRF + dKernel;
        % %method 2
        % dR_residual_clr = dR_clr - dR_total(:, :, :, [2 4]);
        % dR_residual_cld1 = dR_allsky - dR_total(:, :, :, [1 3]); % r terms, indicate cloud feedback and co2/aerosol forcing
        % % dR_cloud = dR_residual_cld1 - dR_residual_clr;
        % dR_residual_cld = dR_allsky - dR_total(:, :, :, [1 3]) - dR_cloud; % r terms, indicate co2/aerosol forcing

        % dR_cloud_sfc = dR_cloud(:, :, :, 1); dR_cloud_toa = dR_cloud(:, :, :, 2);
        % dR_residual_cld_sfc = dR_residual_cld(:, :, :, 1); dR_residual_cld_toa = dR_residual_cld(:, :, :, 2);
        % dR_residual_clr_sfc = dR_residual_clr(:, :, :, 1); dR_residual_clr_toa = dR_residual_clr(:, :, :, 2);
        % % equal to real_dradEfect.mat
        % dR_ObsTotal_cld_sfc = dR_allsky(:, :, :, 1); dR_ObsTotal_cld_toa = dR_allsky(:, :, :, 2);
        % dR_ObsTotal_clr_sfc = dR_clr(:, :, :, 1); dR_ObsTotal_clr_toa = dR_clr(:, :, :, 2);
        % % save dR_cloud and dR_residual_cld and observe total effect
        % save([radefectOutPath, 'dR_cloud_sfc.mat'], 'dR_cloud_sfc');
        % save([radefectOutPath, 'dR_cloud_toa.mat'], 'dR_cloud_toa');
        % save([radefectOutPath, 'dR_residual_cld_sfc.mat'], 'dR_residual_cld_sfc');
        % save([radefectOutPath, 'dR_residual_cld_toa.mat'], 'dR_residual_cld_toa');
        % save([radefectOutPath, 'dR_residual_clr_sfc.mat'], 'dR_residual_clr_sfc');
        % save([radefectOutPath, 'dR_residual_clr_toa.mat'], 'dR_residual_clr_toa');
        % save([radefectOutPath, 'dR_ObsTotal.mat'], 'dR_ObsTotal_cld_sfc', 'dR_ObsTotal_cld_toa', ...
        %     'dR_ObsTotal_clr_sfc', 'dR_ObsTotal_clr_toa');

        % %% step3: cal The effect of Ts on the atmosphere(dR_tstoa-dR_tssfc) (note effect on atoms is positive, but effect on earth is nagetive)
        % % and main Effect (ta+alb+q+clod effect)
        % %( 1sfc cld,2sfc clr,3toa cld,4toa clr)
        % dR_tsAtom_cld = dR_ts(:, :, :, 3) - dR_ts(:, :, :, 1);
        % dR_tsAtom_clr = dR_ts(:, :, :, 4) - dR_ts(:, :, :, 2);
        % save([radefectOutPath, 'dR_tsAtom_cld.mat'], 'dR_tsAtom_cld');
        % save([radefectOutPath, 'dR_tsAtom_clr.mat'], 'dR_tsAtom_clr');
        % dR_mainEffect_sfc = dR_cloud_sfc + dR_main(:, :, :, 1);
        % dR_mainEffect_toa = dR_cloud_toa + dR_main(:, :, :, 3);
        % save([radefectOutPath, 'dR_mainEffect_sfc.mat'], 'dR_mainEffect_sfc');
        % save([radefectOutPath, 'dR_mainEffect_toa.mat'], 'dR_mainEffect_toa');
        % disp([level.model2{level1}, ' model is done!'])
        % disp(' ')
    end

    disp([level.time1{p_1}, ' era is done!'])
    % clear
end

t = toc; disp(t)
