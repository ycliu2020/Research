%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-09 15:52:00
% LastEditTime : 2020-06-25 17:25:56
% LastEditors  : LYC
% Description  : cal mainly include 1.regrid vars, 2.vars anomly
%                CMIP6 mothly data
%                PS: m mean month, s mean season, yr mean year
%                time.date:2000.01-2014.12(interval:15*12);1980.01-2014.12(interval:35*12); 2015.01-2099.12(interval:85*12)
%                initial time.date in amip(432 total): 253 of 432(2000.03);13 of 432(1980.01);
%                initial time.date in futrue(1032 total): 1 of 1032(2015.01);
%                initial time.date in amip-hist(1740 total): 1,561 of 1740(2000.03);1,321 of 1740(1980.01);
% exmPath     : /Research/p2_processCMIP6Data/s2.varsTrend/dSFC_trend.m
% Attention!!!
% check lat: model lat disagree with kernels lat (Opposite direction)
%%---------------------------------------------------------
clear; clc; tic;
nowpath = pwd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment
for p_1 = 1:2%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    % model parameters
    [~, modlist_Experiment, level, tLin, mPlev, vars] = modelParameters(p_1);
    % chose right month, very important !!!
    startmonth = 1;
    % input path
    exmPath = ['/data1/liuyincheng/cmip6-process/', level.time1{p_1}]; %~/data/cmip6/2000-2014/
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model
    for level1 = 1:length(level.model2)% model numbers
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
            % input and output Path
            esmPath=fullfile(mdlPath,esmName{esmNum,1});
            inputPath = [esmPath, '/', level.process3{2}]; %~/data/cmip6/2000-2014/MRI-ESM2-0/ensemble member/anomaly
            outpathname = [esmPath, '/', level.process3{3}]; %/home/lyc/data/cmip6/2000-2014/MIROC6/anomaly_trend
            auto_mkdir(outpathname)
            % load
            load([inputPath, 'global_vars.mat'])
            load([inputPath, 'dts.mat'])
            load([inputPath, 'drlds.mat'])% surface_downwelling_longwave_flux_in_air
            load([inputPath, 'drsds.mat'])% surface_downwelling_shortwave_flux_in_air
            load([inputPath, 'drsus.mat'])% surface_upwelling_shortwave_flux_in_air
            load([inputPath, 'dhfls.mat'])% Surface Upward Latent Heat Flux
            load([inputPath, 'dhfss.mat'])% Surface Upward Sensible Heat Flux
            % unite define a vector component which is positive when directed downward
            dhFlux = -dhfls - dhfss; % LH+SH
            dr_swnet = drsds - drsus; % sfc net shortwave flux
            drhs = drlds + dr_swnet; % equilibrium equation's RHS, nearly equal to sfc upward rad
            drhsPlus = drhs + dhFlux; % should more equal to sfc upward rad
    
            %regrid 144x72(unite grids)
            lat = 88.75:-2.5:-88.75; nlat = length(lat);
            lon = lonf; nlon = length(lon);
            nlonf = length(lonf); nlatf = length(latf);

            dhFlux = autoRegrid3(latf, lonf, time.date, dhFlux, lat, lon, time.date);
            dts = autoRegrid3(latf, lonf, time.date, dts, lat, lon, time.date);
            drhs = autoRegrid3(latf, lonf, time.date, drhs, lat, lon, time.date);
            drhsPlus = autoRegrid3(latf, lonf, time.date, drhsPlus, lat, lon, time.date);
    
            % cal the trend
            [trendm_dts, trends_dts, trendyr_dts, p_dts, cons_dts] = autoCalTrend(dts, nlon, nlat, time.date, startmonth);
            [trendm_drhs, trends_drhs, trendyr_drhs, p_drhs, cons_drhs] = autoCalTrend(drhs, nlon, nlat, time.date, startmonth);
            [trendm_drhsPlus, trends_drhsPlus, trendyr_drhsPlus, p_drhsPlus, cons_drhsPlus] = autoCalTrend(drhsPlus, nlon, nlat, time.date, startmonth);
            [trendm_dhFlux, trends_dhFlux, trendyr_dhFlux, p_dhFlux, cons_dhFlux] = autoCalTrend(dhFlux, nlon, nlat, time.date, startmonth);
    
            % now we done all the job, now save and output.
            save([outpathname, 'trend_dts.mat'], 'trendm_dts', 'cons_dts', 'p_dts', 'trends_dts', 'trendyr_dts');
            save([outpathname, 'trend_drhs.mat'], 'trendm_drhs', 'cons_drhs', 'p_drhs', 'trends_drhs', 'trendyr_drhs');
            save([outpathname, 'trend_drhsPlus.mat'], 'trendm_drhsPlus', 'cons_drhsPlus', 'p_drhsPlus', 'trends_drhsPlus', 'trendyr_drhsPlus');
            save([outpathname, 'trend_dhFlux.mat'], 'trendm_dhFlux', 'cons_dhFlux', 'p_dhFlux', 'trends_dhFlux', 'trendyr_dhFlux');
            save([outpathname, 'global_vars.mat'], 'lon', 'lat', 'time', 'plevf', 'readme', 'timeseries', 'modelname')
            
            disp([esmName{esmNum,1}, ' ensemble is done!'])
        end
        disp([level.model2{level1}, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{p_1}, ' era is done!'])
    disp(' ')
end

t = toc; disp(t)
