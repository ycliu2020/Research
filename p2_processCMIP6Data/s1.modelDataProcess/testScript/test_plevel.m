%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-09 15:52:00
% LastEditTime : 2020-07-02 21:22:08
% LastEditors  : LYC
% Description  : cal mainly include 1.regrid vars, 2.vars anomly
%                CMIP6 mothly data
%                time:2000.01-2014.12(interval:15*12);1980.01-2014.12(interval:35*12); 2015.01-2099.12(interval:85*12)
%                initial time in amip(432 total): 253 of 432(2000.03);13 of 432(1980.01);
%                initial time in futrue(1032 total): 1 of 1032(2015.01);
%                initial time in amip-hist(1740 total): 1,561 of 1740(2000.03);1,321 of 1740(1980.01);
% inputPath     : /code/p2_processCMIP6Data/s1.modelDataProcess/dvar.m
% Attention!!!
% check lat: model lat disagree with kernels lat (Opposite direction)
%%---------------------------------------------------------
clear; clc; tic;
nowpath = pwd;
outPath = '/data1/liuyincheng/cmip6-process/';
var_state = {'d', 'clim_', 'trendm_d', 'trends_d', 'trendyr_d'}; % m,s,yr indicate that month, season, year
kernels_path = '/data1/liuyincheng/y_kernels/kernels_YiH/toa/dp.nc';
lonf = 0:2.5:357.5; nlonf = length(lonf);
latf = 90:-2.5:-90; nlatf = length(latf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment
for p_1 = 1:4%1 mean amip 2000; 2 mean amip 1980; 3 means ssp245, 4 means ssp370, 6 abrupt-4xCO2_150years
    % model parameters
    [readme, Experiment, level, tLin, mPlev, vars] = modelParameters(p_1);
    % experiment path (tLin:1740)
    inputPath = '/data1/liuyincheng/CMIP6-mirror/';
    exmPath_all = cell(1, length(Experiment));

    for ii = 1:length(Experiment)
        exmPath_all{ii} = fullfile(inputPath, Experiment{ii});
    end

    exmPath = exmPath_all{p_1};
    time2 = 1:tLin.inter{p_1};
    plevf = ncread(kernels_path, 'player'); % pay attention to plevf's range must smaller than plev's
    plevfnum = length(plevf);
    readme.timeseries = tLin.read{p_1};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model
    for level1 = 1:length(level.model2)% model number
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
            esmPath=fullfile(mdlPath,esmName{esmNum,1});
            eval(['cd ', esmPath]);
            % outputpath and make it
            % outPathName{1} = fullfile(outPath, level.time1{p_1}, level.model2{level1}, esmName{esmNum,1}, level.process3{1}); %'model/ensemble member/rawdata_regrid/'
            % outPathName{2} = fullfile(outPath, level.time1{p_1}, level.model2{level1}, esmName{esmNum,1}, level.process3{2}); %'model/ensemble member/anomaly/'
            % auto_mkdir(outPathName{1})
            % auto_mkdir(outPathName{2})

            %% check1:  Time line
            % % locate first year
            % temp = dir(fullfile(esmPath,[vars.D3{1}, '_Amon*.nc']));
            % formatOut = 'yyyy-mm';
            % startT = cdftime2loc(temp, formatOut, tLin.start{p_1});
            % % read time
            % startLoc = startT; count = tLin.inter{p_1}; stride = 1;
            % time = cmipTimeRead(temp.name,startLoc,count,stride); %date, Units, Calendar, length 
            % % test time line consistence
            % startYear = str2num(tLin.time{p_1}(1:4)); endYear = str2num(tLin.time{p_1}(end - 3:end));
            % disp('Time line check:')
            % testTime(time.date, startYear, endYear, 1)

            %% check2:  plev
            temp = dir(fullfile(esmPath,[vars.D4{1}, '_Amon*.nc']));
            % test plev consistence
            testPlev(temp, mPlev, p_1)

     
            disp([esmName{esmNum,1}, ' ensemble is done!'])
        end
        disp([level.model2{level1}, ' model is done!'])
        disp(' ')
    end
    disp([level.time1{p_1}, ' era is done!'])
    disp(' ')
end

eval(['cd ', nowpath]);

t = toc; disp(t)
