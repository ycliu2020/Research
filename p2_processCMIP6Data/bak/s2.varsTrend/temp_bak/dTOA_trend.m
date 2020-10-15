%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-09 15:52:00
% LastEditTime : 2020-07-05 17:16:15
% LastEditors  : LYC
% Description  : cal mainly include 1.regrid vars, 2.vars anomly
%                CMIP6 mothly data
%                PS: m mean month, s mean season, yr mean year
%                time.date:2000.01-2014.12(interval:15*12);1980.01-2014.12(interval:35*12); 2015.01-2099.12(interval:85*12)
%                initial time.date in amip(432 total): 253 of 432(2000.03);13 of 432(1980.01);
%                initial time.date in futrue(1032 total): 1 of 1032(2015.01);
%                initial time.date in amip-hist(1740 total): 1,561 of 1740(2000.03);1,321 of 1740(1980.01);
% exmPath     : /Research/p2_processCMIP6Data/s2.varsTrend/dTOA_trend.m
% Attention!!!
% check lat_f: model lat_f disagree with kernels lat_f (Opposite direction)
%%---------------------------------------------------------
% note that TOA use net rad flux
clear; clc; tic;
nowpath = pwd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment
for p_1 = 4:4%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    % model parameters
    [~, modlist_Experiment, level, tLin, mPlev, vars] = cmipParameters(p_1);
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
            anomTrendPath = [esmPath, '/', level.process3{3}]; %/home/lyc/data/cmip6/2000-2014/MIROC6/anomaly_trend
            auto_mkdir(anomTrendPath)
            % load
            load([inputPath, 'global_vars.mat'])
            load([inputPath, 'drlut.mat'])% toa_outgoing_longwave_flux
            load([inputPath, 'drsut.mat'])% toa_outgoing_shortwave_flux
            load([inputPath, 'drsdt.mat'])% toa_incoming_shortwave_flux
            % unite define a vector component which is positive when directed downward
            dnetTOA = drsdt - drlut - drsut; % net radiation in toa
            %regrid 144x72(unite grids)
            lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f);
            lon_f = lon_k; nlonf = length(lon_f);
            nlonk = length(lon_k); nlatk = length(lat_k);
            dnetTOA = autoRegrid3(lon_k, lat_k, time.date, dnetTOA, lon_f, lat_f, time.date);
            % cal the trend
            [trendm_dnetTOA, trends_dnetTOA, trendyr_dnetTOA, p_dnetTOA, cons_dnetTOA] = autoCalTrend(dnetTOA, nlonf, nlatf, time.date, startmonth);
            % now we done all the job, now save and output.
            save([anomTrendPath, 'trend_dnetTOA.mat'], 'trendm_dnetTOA', 'cons_dnetTOA', 'p_dnetTOA', 'trends_dnetTOA', 'trendyr_dnetTOA');
            
            disp([esmName{esmNum,1}, ' ensemble is done!'])
        end

        disp([level.model2{level1}, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{p_1}, ' era is done!'])
    disp(' ')
end

t = toc; disp(t)
