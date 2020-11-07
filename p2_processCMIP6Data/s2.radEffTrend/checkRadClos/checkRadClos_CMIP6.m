%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-10-30 18:26:55
% LastEditTime : 2020-11-01 20:01:46
% LastEditors  : LYC
% Description  : test rad closure of CMIP6 model on sfc and toa
% FilePath     : /code/p2_processCMIP6Data/s1.modelDataProcess/checkRadClos/checkRadClos_CMIP6.m
% symbol_custom_string_obkoro1:
%%---------------------------------------------------------
clc; clear;
nowpath = pwd;

%% calculate mask file
% load mask map
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/中国干湿以及青藏高原分区地图/main_use/mask720.mat')
% 国境线
bou_china = shaperead('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/map_data/中国及各省shp/国界线/bou1_4p.shp'); % load china  boundary
bou_chinaX = [bou_china(:).X]; bou_chinaY = [bou_china(:).Y];
% 国境线(带南海线)
bou_china_line = shaperead('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/map_data/中国及各省shp/国界线/bou1_4l.shp'); % load china  boundary
bou_china_lineX = [bou_china_line(:).X]; bou_china_lineY = [bou_china_line(:).Y];
% 省界线
bou_chinaProvince = shaperead('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/map_data/中国及各省shp/省界线/bou2_4p.shp'); % load china  boundary
bou_chinaProvinceX = [bou_chinaProvince(:).X]; bou_chinaProvinceY = [bou_chinaProvince(:).Y];

%save maskFile to area of china east
maskchina_east = maskchina_cp;
maskchina_east(lonw < 112, :) = 0;
maskchina_east(:, latw > 38) = 0;
readme_mask_CN_east = 'area of east China, Range: lonw>112E,latw<38N; Based on maskchina_cp(lon, lat same as it)';
save('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472_east.mat', 'maskchina_east', 'lat14472', 'lon14472', 'readme_mask_CN_east')

%% loop and process
[mlabels, areaNum] = cmipPlotParameters('atm', 'land', 'radEffect'); % plot parameters

esm = 'r1i1p1f1';
exmStart = 1; exmEnd = 1;

for exmNum = exmStart:exmEnd
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % mPath.input:E:/data/cmip6-process/2000-2014/
    exmPath = fullfile('/data1/liuyincheng/cmip6-process/', level.time1{exmNum});
    % mPath.output:a_research/P02.Ts_change_research/figure/04.cmip6Result/2000-2014/ 根据输出的格式进行修改
    % mPath.uniOutput = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/02.cmip6Result/1.6/', lower(mlabels.level), level.time1{exmNum});
    % mPath.Output = fullfile(mPath.uniOutput);
    % auto_mkdir(mPath.Output)

    % model loop
    mdlStart = 1; mdlEnd = length(level.model2); % differnt models%length(level.model2)

    for mdlNum = mdlStart:mdlEnd
        % model path
        mdlName = level.model2{mdlNum};
        mdlPath = fullfile(exmPath, level.model2{mdlNum});
        eval(['cd ', mdlPath]);
        disp(' ')
        disp([level.model2{mdlNum}, ' model start!'])

        % ensemble member path
        esmName = getPath_fileName(mdlPath, '.');
        eval(['cd ', nowpath]);
        %% 暂时只看esm实验
        esm = 'r1i1p1f1';

        if sum(strcmp(esmName, esm)) == 0
            disp(['the ', esm, ' ensemble of ', mdlName, ' didnt exist']);
            continue
        end

        specificNum = find(strcmp(esmName, 'r1i1p1f1') == 1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensemble member
        for esmNum = specificNum:specificNum%1:1%length(esmName)% note that r1i1p1 sometime not the first folder
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% load and read
            esmPath = fullfile(mdlPath, esmName{esmNum, 1});
            dradTrendPath = fullfile(esmPath, level.process3{7}); %/data1/liuyincheng/cmip6-process/amip_1980-2014/CESM2/Effect_trend
            dvarsPath = fullfile(esmPath, level.process3{2}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
            vsTsEffectPath = fullfile(esmPath, level.process3{9}); %/data1/liuyincheng/cmip6-process/amip_1980-2014/CESM2/vsTsEffect
            load([dradTrendPath, 'global_vars.mat'])% lat_f lon_f time plevf readme
            load([dvarsPath, 'global_vars.mat'])% lat_k lon_k time plevk readme
            load([dvarsPath, 'dts.mat'])% dts clim_ts
            load([vsTsEffectPath, 'dTs_x_sfc.mat'])% dTs_alb, dTs_cloud, dTs_hus, dTs_residual, dTs_ta
            nlonf = length(lon_f); nlatf = length(lat_f); ntime = length(time.date);
            dts = autoRegrid3(lon_k, lat_k, time.date, dts, lon_f, lat_f, time.date);

            % read model out put
            load([dradEffectPath, 'real_dradEfect.mat'])% 'l_rad', 's_rad', 'dR_allsky', 'dR_clr', 'readme_realradEfect'
            % cal rhs(RHeating)
            load([dvarsPath, 'drlds.mat'])% surface_downwelling_longwave_flux_in_air
            load([dvarsPath, 'drsds.mat'])% surface_downwelling_shortwave_flux_in_air
            load([dvarsPath, 'drsus.mat'])% surface_upwelling_shortwave_flux_in_air
            dR_swnet = drsds - drsus; % sfc net shortwave flux
            drhs = drlds + dR_swnet; % equilibrium equation's RHS, nearly equal to sfc upward rad
            drhs = autoRegrid3(lon_k, lat_k, time.date, drhs, lon_f, lat_f, time.date);
            % read sfc values
            dradEffectPath = fullfile(esmPath, level.process3{6}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
            load([dradEffectPath, 'dR_cloud_sfc.mat'])% dR_cloud_sfc
            load([dradEffectPath, 'dradEfect_sfc_cld.mat'])% albEffect, husEffect, mainEffect, taEffect, taOnlyEffect, taOnlyEffect2, tasEffect, tasEffect2, totalEffect, tsEffect, wvlwEffect, wvswEffect
            load([dradEffectPath, 'dR_residual_cld_sfc.mat'])% dR_resiual_cld_sfc
            albEffect_sfc = albEffect;
            husEffect_sfc = husEffect;
            taEffect_sfc = taEffect;
            tsEffect_sfc = tsEffect;
            mainEffect_sfc = mainEffect;
            % read toa values
            load([dradEffectPath, 'dR_cloud_toa.mat'])% dR_cloud_toa
            load([dradEffectPath, 'dradEfect_toa_cld.mat'])% albEffect, husEffect, mainEffect, taEffect, taOnlyEffect, taOnlyEffect2, tasEffect, tasEffect2, totalEffect, tsEffect, wvlwEffect, wvswEffect
            load([dradEffectPath, 'dR_residual_cld_toa.mat'])% dR_resiual_cld_toa
            albEffect_toa = albEffect;
            husEffect_toa = husEffect;
            taEffect_toa = taEffect;
            tsEffect_toa = tsEffect;
            mainEffect_toa = mainEffect;


            varUsed = zeros(nlonf, nlatf, ntime, 2);
            varUsed(:, :, :, 1) = squeeze(dR_allsky(:, :, :, 1)); % sfc net
            varUsed(:, :, :, 2) = mainEffect + dR_cloud_sfc + tsEffect; % kernel cal sfc net
            varUsed(:, :, :, 3) = dR_residual_cld_sfc; % residual
            size_varUsed = size(varUsed);
            nValue = size_varUsed(4);

            % mask
            latRange = 88.75;
            areaNum = 'world';

            for x_ind = 1:nValue
                [varUsed(:, :, :, x_ind), ~, ~] = maskArea(varUsed(:, :, :, x_ind), lat_f, latRange, -latRange, areaNum);
            end

            % cal the weight average value of east china
            glb_varUsed = zeros(ntime, nValue);

            for timeNum = 1:ntime

                for x_ind = 1:nValue
                    glb_varUsed(timeNum, x_ind) = areaMeanLatWeight(varUsed(:, :, timeNum, x_ind), lat_f);
                end

            end

            % cal the coef

            disp([esmName{esmNum, 1}, ' ensemble is done!'])

        end

        disp([level.model2{mdlNum}, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{exmNum}, ' era is done!'])
    disp(' ')
end

eval(['cd ', nowpath]);
