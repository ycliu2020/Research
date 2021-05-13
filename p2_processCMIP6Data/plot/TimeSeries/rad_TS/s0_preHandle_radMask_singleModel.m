%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-08-31 17:00:15
% LastEditTime : 2021-05-10 16:36:04
% LastEditors  : Please set LastEditors
% Description  : MME result of 时间序列图
% FilePath     : /code/p2_processCMIP6Data/plot/TimeSeries/rad_TS/s0_preHandle_radMask_singleModel.m
% note :
%%---------------------------------------------------------
clc; clear; tic;
nowpath = pwd;
% load mask map
run '/home/liuyc/lib/tools/matlab/myTools/autoScript/preLoadVar.m'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CMIP6 part
% experiment
for exmNum = [1 2 4] %1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    %CAMS-CSM1-0 didn't have sfc clear sky radiation, delete it
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % exmPath
    exmPath = ['/data1/liuyincheng/CMIP6-process/', level.time1{exmNum}]; %/data1/liuyincheng/CMIP6-process/2000-2014/

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model
    for mdlNum = 1:length(level.model2)
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
        for esmNum = specificNum:specificNum %1:length(esmName)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% load and read
            esmPath = fullfile(mdlPath, esmName{esmNum, 1});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% load and read
            % data path
            dvarsPath = fullfile(esmPath, level.process3{2}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/anomaly
            dradEffectPath = fullfile(esmPath, level.process3{6}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/radEffect/
            load([dradEffectPath, 'global_vars.mat']) % 'lon_f', 'lat_f', 'timeEssmble', 'time', 'plev_k', 'readme', 'timeseries', 'MME_Models'
            ntime = length(time.date);

            % cal rhs(RHeating)
            load([dvarsPath, 'drlds.mat']) % surface_downwelling_longwave_flux_in_air
            load([dvarsPath, 'drsds.mat']) % surface_downwelling_shortwave_flux_in_air
            load([dvarsPath, 'drsus.mat']) % surface_upwelling_shortwave_flux_in_air
            load([dvarsPath, 'drlus.mat']) % surface_upwelling_longwave_flux_in_air
            % LH and SH
            load([dvarsPath, 'dhfss.mat']) % Surface Upward Sensible Heat Flux
            load([dvarsPath, 'dhfls.mat']) % Surface Upward Latent Heat Flux

            dR_swnet = drsds - drsus; % sfc net shortwave flux
            drhs = drlds + dR_swnet; % equilibrium equation's RHS, nearly equal to sfc upward rad
            drhs = autoRegrid3(lon_k, lat_k, time.date, drhs, lon_f, lat_f, time.date);
            dlwUPRawData = autoRegrid3(lon_k, lat_k, time.date, drlus, lon_f, lat_f, time.date);

            dhfss = autoRegrid3(lon_k, lat_k, time.date, dhfss, lon_f, lat_f, time.date);
            dhfls = autoRegrid3(lon_k, lat_k, time.date, dhfls, lon_f, lat_f, time.date);

            % cloud radRffect
            load([dradEffectPath, 'dR_cloud.mat'])

            % radEffect of models
            load([dradEffectPath, 'dradEffect_sfc_cld.mat']) % 'totalEffect', 'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', 'taEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect', 'mainEffect'
            load([dradEffectPath, 'dR_residual_cld_sfc.mat']) % dR_residual_cld_sfc, lw_residual_cld_sfc, sw_residual_cld_sfc
            drhsKern = albEffect + husEffect + taEffect + dR_cloud_sfc;
            dlwUPKern = -tsEffect;

            time = time.date;
            ntime = length(time);
            startMonth = 1;

            run regnMask.m
            
            % save value
            outPutPath = fullfile(esmPath, 'FigData/');
            auto_mkdir(outPutPath)
            save([outPutPath, 'regionalVarsRad_sfc_cld.mat'], 'dR_cloudMask', 'dR_residualMask', 'dR_albMask', 'dR_husMask', 'dR_taMask', 'dR_tsMask', 'dhfssMask', 'dhflsMask', ...
                'drhsMask', 'drhsKernMask', 'dlwUPRawDataMask', 'dlwUPKernMask')
        end

    end

end

eval(['cd ', nowpath]);
t = toc; disp(t)
