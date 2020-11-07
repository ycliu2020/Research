%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-09 15:52:00
% LastEditTime : 2020-11-07 10:35:44
% LastEditors  : LYC
% Description  : includ sfc and toa trend of vars
%                cal mainly include 1.regrid vars, 2.vars anomly
%                CMIP6 mothly data
%                PS: m mean month, s mean season, yr mean year
%                time.date:2000.01-2014.12(interval:15*12);1980.01-2014.12(interval:35*12); 2015.01-2099.12(interval:85*12)
%                initial time.date in amip(432 total): 253 of 432(2000.03);13 of 432(1980.01);
%                initial time.date in futrue(1032 total): 1 of 1032(2015.01);
%                initial time.date in amip-hist(1740 total): 1,561 of 1740(2000.03);1,321 of 1740(1980.01);
% exmPath     :
% Attention!!!
% check lat_f: model lat_f disagree with kernels lat_f (Opposite direction)
%%---------------------------------------------------------
clear; clc; tic;
nowpath = pwd;

% 开启并行环境
poolobj = gcp('nocreate'); % If no pool,  create new one.

if isempty(poolobj)
    MyPar = parpool(16);
else
    MyPar = gcp('nocreate');
    disp('Already initialized'); %说明并行环境已经启动。
end

%% 预选读取所有的路径
exm1 = 1; exm2 = 4;
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
    lon_k = 0:2.5:357.5; nlonk = length(lon_k);
    lat_k = 90:-2.5:-90; nlatk = length(lat_k);
    lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
    lon_f = lon_k; nlonf = length(lon_f);

    % model parameters
    [~, modlist_Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % chose right month, very important !!!
    startmonth = 1;
    % ensemble member path
    esmName = getPath_fileName(mdlPath, '.');
    eval(['cd ', esmPath]);
    disp([level.model2{mdlNum}, ', ', level.time1{exmNum}, ', ', esmName{esmNum, 1}, ' ensemble start!'])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ensemble member
    for esmNum = 1:length(esmName)
        % input and output Path
        esmPath = fullfile(mdlPath, esmName{esmNum, 1});
        inputPath = [esmPath, '/', level.process3{2}]; %~/data/CMIP6/2000-2014/MRI-ESM2-0/ensemble member/anomaly
        anomTrendPath = [esmPath, '/', level.process3{3}]; %/home/lyc/data/CMIP6/2000-2014/MIROC6/anomaly_trend
        auto_mkdir(anomTrendPath)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SFC vars trend
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
        dR_swnet = drsds - drsus; % sfc net shortwave flux
        drhs = drlds + dR_swnet; % equilibrium equation's RHS, nearly equal to sfc upward rad
        drhsPlus = drhs + dhFlux; % should more equal to sfc upward rad
        %regrid 144x72(unite grids)
        dhFlux = autoRegrid3(lon_k, lat_k, time.date, dhFlux, lon_f, lat_f, time.date);
        dts = autoRegrid3(lon_k, lat_k, time.date, dts, lon_f, lat_f, time.date);
        drhs = autoRegrid3(lon_k, lat_k, time.date, drhs, lon_f, lat_f, time.date);
        drhsPlus = autoRegrid3(lon_k, lat_k, time.date, drhsPlus, lon_f, lat_f, time.date);
        % cal the trend
        [trendm_dts, trends_dts, trendyr_dts, p_dts, cons_dts] = autoCalTrend(dts, nlonf, nlatf, time.date, startmonth);
        [trendm_drhs, trends_drhs, trendyr_drhs, p_drhs, cons_drhs] = autoCalTrend(drhs, nlonf, nlatf, time.date, startmonth);
        [trendm_drhsPlus, trends_drhsPlus, trendyr_drhsPlus, p_drhsPlus, cons_drhsPlus] = autoCalTrend(drhsPlus, nlonf, nlatf, time.date, startmonth);
        [trendm_dhFlux, trends_dhFlux, trendyr_dhFlux, p_dhFlux, cons_dhFlux] = autoCalTrend(dhFlux, nlonf, nlatf, time.date, startmonth);
        % now we done all the job, now save and output.
        save([anomTrendPath, 'trend_dts.mat'], 'trendm_dts', 'cons_dts', 'p_dts', 'trends_dts', 'trendyr_dts');
        save([anomTrendPath, 'trend_drhs.mat'], 'trendm_drhs', 'cons_drhs', 'p_drhs', 'trends_drhs', 'trendyr_drhs');
        save([anomTrendPath, 'trend_drhsPlus.mat'], 'trendm_drhsPlus', 'cons_drhsPlus', 'p_drhsPlus', 'trends_drhsPlus', 'trendyr_drhsPlus');
        save([anomTrendPath, 'trend_dhFlux.mat'], 'trendm_dhFlux', 'cons_dhFlux', 'p_dhFlux', 'trends_dhFlux', 'trendyr_dhFlux');
        save([anomTrendPath, 'global_vars.mat'], 'lon_f', 'lat_f', 'time', 'plev_k', 'readme', 'timeseries', 'modelname')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TOA vars trend
        % load
        load([inputPath, 'global_vars.mat'])
        load([inputPath, 'drlut.mat'])% toa_outgoing_longwave_flux
        load([inputPath, 'drsut.mat'])% toa_outgoing_shortwave_flux
        load([inputPath, 'drsdt.mat'])% toa_incoming_shortwave_flux
        % unite define a vector component which is positive when directed downward
        dnetTOA = drsdt - drlut - drsut; % net radiation in toa
        %regrid 144x72(unite grids)
        dnetTOA = autoRegrid3(lon_k, lat_k, time.date, dnetTOA, lon_f, lat_f, time.date);
        % cal the trend
        [trendm_dnetTOA, trends_dnetTOA, trendyr_dnetTOA, p_dnetTOA, cons_dnetTOA] = autoCalTrend(dnetTOA, nlonf, nlatf, time.date, startmonth);
        % now we done all the job, now save and output.
        save([anomTrendPath, 'trend_dnetTOA.mat'], 'trendm_dnetTOA', 'cons_dnetTOA', 'p_dnetTOA', 'trends_dnetTOA', 'trendyr_dnetTOA');
        
        disp([level.model2{mdlNum}, ', ', level.time1{exmNum}, ', ', esmName{esmNum, 1}, ' ensemble is done!'])
        disp(' ')
    end

end
