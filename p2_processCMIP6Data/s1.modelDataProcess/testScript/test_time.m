%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% mothly data
% cal the anomly
% hope to use an unite scrip to read and process all the varies we need
% Compare the pre rad and ts data from 2000-03 to 2018-02(18 years)
%
% time:2000.01-2014.12(interval:15*12);1980.01-2014.12(interval:35*12); 2015.01-2099.12(interval:85*12)
% initial time in amip(432 total): 253 of 432(2000.03);13 of 432(1980.01);
% initial time in futrue(1032 total): 1 of 1032(2015.01);
% initial time in amip-hist(1740 total): 1,561 of 1740(2000.03);1,321 of 1740(1980.01);
%
% cal mainly include 1.regrid vars, 2.vars anomly
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Attention!!!
% check lat: model lat disagree with kernels lat (Opposite direction)
% check plev drection and unite: model plev disagree with kernels plev
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear; clc; tic;
nowpath = pwd;

% file structure
% set time parameter
for p_1 = 1:4%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    % model parameters
    [readme, Experiment, level, tLin, vars] = modelParameters(p_1);
    % input and output path (tLin:1740)
    inputPath = '/data1/liuyincheng/CMIP6-mirror/';
    filepath_all = cell(1, length(Experiment));

    for ii = 1:length(Experiment)
        filepath_all{ii} = [inputPath, Experiment{ii}, '/'];
    end

    filepath = filepath_all{p_1};
    outPath = '/data1/liuyincheng/cmip6-process/';
    var_state = {'d', 'clim_', 'trendm_d', 'trends_d', 'trendyr_d'}; % m,s,yr indicate that month, season, year
    kernels_path = '/data/pub/kernels_YiH/toa/dp.nc';
    time2 = 1:tLin.inter(p_1);
    plevf = ncread(kernels_path, 'player'); % pay attention to plevf's range must smaller than plev's
    plevfnum = length(plevf);
    lonf = 0:2.5:357.5; nlonf = length(lonf);
    latf = 90:-2.5:-90; nlatf = length(latf);
    readme.timeseries = tLin.read{p_1};

    for level1 = 1:length(level.model2)% model number
        % inputpath
        filepath_t = [filepath, level.model2{level1}, '/'];
        eval(['cd ', filepath_t]);
        % outputpath and make it
        outPathName{1} = [outPath, level.time1{p_1}, level.model2{level1}, '/', level.process3{1}];
        outPathName{2} = [outPath, level.time1{p_1}, level.model2{level1}, '/', level.process3{2}];
        auto_mkdir(outPathName{1})
        auto_mkdir(outPathName{2})

        %% Time line test
        % locate first year
        temp = dir([filepath_t, vars.D3{1}, '_Amon*.nc']);
        formatOut = 'yyyy-mm';
        startT = cdftime2loc(temp, formatOut, tLin.start{p_1});
        % read time
        startLoc = startT; count = tLin.inter(p_1); stride = 1;
        time = ncread(temp.name, 'time', startLoc, count, stride); %unite:day
        timeUnits = ncreadatt(temp.name, 'time', 'units');
        timeCalendar = ncreadatt(temp.name, 'time', 'calendar');
        time = cdfdate2num(timeUnits, timeCalendar, time);
        % test
        startYear = str2num(tLin.time{p_1}(1:4)); endYear = str2num(tLin.time{p_1}(end - 3:end));
        ntime = length(time);
        disp(' ')
        disp([level.model2{level1}, ' model start!'])
        disp('Time line check:')
        testTime(time, startYear, endYear, 1)
    end

    disp([level.time1{p_1}, ' era is done!'])
    disp(' ')
end

eval(['cd ', nowpath]);

t = toc; disp(t)
