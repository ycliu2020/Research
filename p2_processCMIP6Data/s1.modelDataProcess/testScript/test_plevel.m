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
for p_1 = 4:4%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    % model parameters
    [readme, Experiment, level, tLin, mPlev, vars] = modelParameters(p_1);
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
        disp(' ')
        disp([level.model2{level1}, ' model start!'])
        auto_mkdir(outPathName{1})
        auto_mkdir(outPathName{2})

        %% Time line test
        % locate first year
        temp = dir([filepath_t, vars.D3{1}, '_Amon*.nc']);
        % formatOut = 'yyyy-mm';
        % startT = cdftime2loc(temp, formatOut, tLin.start{p_1});
        % % read time
        % startLoc = startT; count = tLin.inter(p_1); stride = 1;
        % time = ncread(temp.name, 'time', startLoc, count, stride); %unite:day
        % timeUnits = ncreadatt(temp.name, 'time', 'units');
        % timeCalendar = ncreadatt(temp.name, 'time', 'calendar');
        % time = cdfdate2num(timeUnits, timeCalendar, time);
        % % test
        % startYear = str2num(tLin.time{p_1}(1:4)); endYear = str2num(tLin.time{p_1}(end - 3:end));
        % ntime = length(time);
        % disp('Time line check:')
        % testTime(time, startYear, endYear, 1)

        %% 2) read data
        % 4D vars
        v4_names = {'hus', 'ta'};
        % define anomaly\mean\trend data
        dv4_names = strcat(var_state{1}, v4_names); clim_v4_names = strcat(var_state{2}, v4_names);
        % read 4-D variables
        temp = dir([filepath_t, v4_names{1}, '_Amon*.nc']);
        lon_v4 = ncread(temp.name, 'lon'); lat_v4 = ncread(temp.name, 'lat');
        lat_v4 = lat_v4(end:-1:1);

        % plev.Units = ncreadatt(temp.name, 'plev', 'units')
        % plev.Wards = ncreadatt(temp.name, 'plev', 'positive')
        % plev.data=ncread(temp.name, 'plev');
        % length(plev.data)
        % max(plev.data)
        % min(plev.data)
        testPlev(temp,mPlev,p_1)


        % if strcmp(level.model2{level1}(1:end), 'CAMS-CSM1-0') == 1
        %     plev = ncread(temp.name, 'plev');
        %     plev = plev(1:17) / 100; % paùùhpa,check the unite and direction
        %     startLoc = [1 1 1 startT]; count = [inf inf 17 tLin.inter(p_1)]; stride = [1 1 1 1];
        % else
        %     plev = ncread(temp.name, 'plev') / 100; %
        %     startLoc = [1 1 1 startT]; count = [inf inf inf tLin.inter(p_1)]; stride = [1 1 1 1];
        % end

        % temp_v4 = zeros(length(v4_names), length(lon_v4), length(lat_v4), length(plev), tLin.inter(p_1));

        % for i1 = 1:length(v4_names)%2vars
        %     temp = dir([filepath_t, v4_names{i1}, '_Amon*.nc']);
        %     temp_v4(i1, :, :, :, :) = ncread(temp.name, v4_names{i1}, startLoc, count, stride);
        % end

        % temp_v4 = temp_v4(:, :, end:-1:1, :, :); % transfor the lat
    end
        

    disp([level.time1{p_1}, ' era is done!'])
    disp(' ')
end

eval(['cd ', nowpath]);

t = toc; disp(t)

% %read 4-D variables
% filepath_t = [filepath, level.model2{level1}];
% eval(['cd ', filepath_t]);
% startLoc = [1 1 1 1515]; count = [inf inf inf 218]; stride = [1 1 1 1];
% dir_t1 = dir([filepath_t, v4{1}, '_Amon*.nc']); %filename=[filepath_t,dir_t1.name];
% dir_t2 = dir([filepath_t, v4{2}, '_Amon*.nc']); %filename=[filepath_t,dir_t1.name];
% hus = ncread(dir_t1.name, v4{1}, startLoc, count, stride);
% ta = ncread(dir_t2.name, v4{2}, startLoc, count, stride);
