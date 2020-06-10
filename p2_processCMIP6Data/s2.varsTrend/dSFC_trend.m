%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-09 15:52:00
% LastEditTime : 2020-06-10 15:18:19
% LastEditors  : LYC
% Description  : 
% FilePath     : /Research/p2_processCMIP6Data/s2.varsTrend/dSFC_trend.m
%  
%%---------------------------------------------------------
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% mothly data
% cal the anomly trend about any vars in CMIP6 as showed on script name
% this script mainly include:trendm_dts,trends_dts,trendyr_dts; trendm_drhs,trends_drhs,trendyr_drhs...
% PS: m mean month, s mean season, yr mean year
%
% hope to use an unite scrip to read and process all the varies we need
% Compare the pre rad and ts data from 2000-03 to 2018-02(18 years)
%
% experiment information:
% time:2000.01-2014.12(interval:15*12);1980.01-2014.12(interval:35*12); 2015.01-2099.12(interval:85*12)
% initial time in hist(1740 total): 1,561 of 1740(2000.03);1,321 of 1740(1980.01)
% initial time in futrue(1032 total): 1 of 1032(2015.01);
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear; clc; tic;
nowpath = pwd;

for p_1 = 1:4 %1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
% model parameters
[~,modlist_Experiment,level,tLin, mPlev, vars] = modelParameters(p_1);
% chose right month, very important !!! 
startmonth=1;
% input path
inputPath = ['/data1/liuyincheng/cmip6-process/', level.time1{p_1}];%~/data/cmip6/2000-2014/

for level1 =  1:length(level.model2)% model numbers
    inputPath1 = [inputPath, level.model2{level1}, '/', level.process3{2}]; %~/data/cmip6/2000-2014/MRI-ESM2-0/anomaly
    load([inputPath1, 'global_vars.mat'])
    load([inputPath1, 'dts.mat'])
    load([inputPath1, 'drlds.mat'])% surface_downwelling_longwave_flux_in_air
    load([inputPath1, 'drsds.mat'])% surface_downwelling_shortwave_flux_in_air
    load([inputPath1, 'drsus.mat'])% surface_upwelling_shortwave_flux_in_air
    load([inputPath1, 'dhfls.mat'])% Surface Upward Latent Heat Flux
    load([inputPath1, 'dhfss.mat'])% Surface Upward Sensible Heat Flux
    % unite define a vector component which is positive when directed downward
    dhFlux = -dhfls - dhfss; % LH+SH
    dr_swnet = drsds - drsus;% sfc net shortwave flux 
    drhs = drlds + dr_swnet;% equilibrium equation's RHS, nearly equal to sfc upward rad
    drhsPlus = drhs + dhFlux;% should more equal to sfc upward rad

    %regrid 144x72(unite grids)
    lat = 88.75:-2.5:-88.75; nlat = length(lat);
    lon=lonf;nlon = length(lon);
    nlonf=length(lonf);nlatf=length(latf);
    [Xlon, Ylat, Ttime] = meshgrid(lat, lon, time);
    [Xlonf, Ylatf, Ttimef] = meshgrid(latf, lonf, time);
    dhFlux = interp3(Xlonf, Ylatf, Ttimef, dhFlux, Xlon, Ylat, Ttime);
    dts = interp3(Xlonf, Ylatf, Ttimef, dts, Xlon, Ylat, Ttime);
    drhs = interp3(Xlonf, Ylatf, Ttimef, drhs, Xlon, Ylat, Ttime);
    drhsPlus = interp3(Xlonf, Ylatf, Ttimef, drhsPlus, Xlon, Ylat, Ttime);
    
    % cal the trend
    [trendm_dts, trends_dts,trendyr_dts, p_dts, cons_dts] = autoCalTrend(dts, nlon, nlat, time, startmonth);
    [trendm_drhs, trends_drhs,trendyr_drhs, p_drhs, cons_drhs] = autoCalTrend(drhs, nlon, nlat, time, startmonth);
    [trendm_drhsPlus, trends_drhsPlus,trendyr_drhsPlus, p_drhsPlus, cons_drhsPlus] = autoCalTrend(drhsPlus, nlon, nlat, time, startmonth);
    [trendm_dhFlux, trends_dhFlux,trendyr_dhFlux, p_dhFlux, cons_dhFlux] = autoCalTrend(dhFlux, nlon, nlat, time, startmonth);

    % now we done all the job, now save and output.
    outpathname = [inputPath ,level.model2{level1}, '/', level.process3{3}];%/home/lyc/data/cmip6/2000-2014/MIROC6/anomaly_trend
    auto_mkdir(outpathname)
    save([outpathname,'trend_dts.mat'],'trendm_dts',  'cons_dts',  'p_dts', 'trends_dts', 'trendyr_dts');
    save([outpathname,'trend_drhs.mat'], 'trendm_drhs', 'cons_drhs', 'p_drhs', 'trends_drhs', 'trendyr_drhs');
    save([outpathname,'trend_drhsPlus.mat'], 'trendm_drhsPlus', 'cons_drhsPlus', 'p_drhsPlus', 'trends_drhsPlus', 'trendyr_drhsPlus');
    save([outpathname,'trend_dhFlux.mat'], 'trendm_dhFlux', 'cons_dhFlux', 'p_dhFlux', 'trends_dhFlux', 'trendyr_dhFlux');
    save([outpathname,'global_vars.mat'], 'lon', 'lat', 'time','plevf','readme','timeseries','modelname')
    disp([level.model2{level1},' model is done!'])
    disp(' ')
end
disp([level.time1{p_1},' era is done!'])
disp(' ')
end
t = toc; disp(t)
