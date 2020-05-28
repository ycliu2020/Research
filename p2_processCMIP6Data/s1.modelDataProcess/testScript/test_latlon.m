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

    for level1 = 5:5%1:length(level.model2)% model number
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
        formatOut = 'yyyy-mm';
        startT = cdftime2loc(temp, formatOut, tLin.start{p_1});
        % read time
        startLoc = startT; count = tLin.inter(p_1); stride = 1;
        time = ncread(temp.name, 'time', startLoc, count, stride); %unite:day
        timeUnits = ncreadatt(temp.name, 'time', 'units');
        timeCalendar = ncreadatt(temp.name, 'time', 'calendar');
        time = cdfdate2num(timeUnits, timeCalendar, time);
        % % test time line consistence
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
        plev = ncread(temp.name, 'plev') / 100; %
        startLoc = [1 1 1 startT]; count = [inf inf inf tLin.inter(p_1)]; stride = [1 1 1 1];
        temp_v = ncread(temp.name, v4_names{1}, startLoc, count, stride);
        temp_v = temp_v(:, end:-1:1, :, :); % transfor the lat
        disp('4_D')
        [pp]=testLonlat(temp);
        p_2=4;
        [pp, lon_v4, lat_v4,temp_v] = fixLonlat(pp, p_2, lon_v4, lat_v4, temp_v, temp);

        % disp('3_D')
        % temp = dir([filepath_t, vars.D3{1}, '_Amon*.nc']);
        % [pp]=testLonlat(temp);
        % lon_v.units = ncreadatt(temp.name, 'lon', 'units');
        % lon_v.data = ncread(temp.name, 'lon');
        % lat_v.units = ncreadatt(temp.name, 'lat', 'units');
        % lat_v.data = ncread(temp.name, 'lat');

        % test plev consistence
        % testPlev(temp,mPlev,p_1)

        % if max(lon_v3) < 357.5%strcmp(level.model2{level1}(1:end),'CanESM5') == 1
        %     temp_v3(:, end + 1, :, :) = temp_v3(:, 1, :, :);
        %     lon_v3(end + 1) = 360;
        % elseif max(lon_v4) < 357.5%strcmp(level.model2{level1}(1:end),'CanESM5') == 1
        %     temp_v4(:, end + 1, :, :, :) = temp_v4(:, 1, :, :, :);
        %     lon_v4(end + 1) = 360;
        % elseif strcmp(level.model2{level1}(1:end), 'GISS-E2-1-G') == 1
        %     temp_vv4 = temp_v4; temp_vv3 = temp_v3;
        %     temp_vv4(:, end + 1, :, :, :) = 0; temp_vv3(:, end + 1, :, :) = 0;
        %     temp_vv4(:, 1, :, :, :) = temp_v4(:, end, :, :, :);
        %     temp_vv4(:, 2:end, :, :, :) = temp_v4;
        %     temp_vv3(:, 1, :, :) = temp_v3(:, end, :, :);
        %     temp_vv3(:, 2:end, :, :) = temp_v3;
        %     lon_v3 = [-1.25; lon_v3]; lon_v4 = [-1.25; lon_v4];
        %     temp_v3 = temp_vv3; temp_v4 = temp_vv4;
        %     clear temp_vv3 temp_vv4
        % end


        % plev = ncread(temp.name, 'plev') / 100; %
        % startLoc = [1 1 1 startT]; count = [inf inf inf tLin.inter(p_1)]; stride = [1 1 1 1];
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
% if strcmp(level.model2{level1}(1:end), 'CAMS-CSM1-0') == 1
%     plev = ncread(temp.name, 'plev');
%     plev = plev(1:17) / 100; % pahpa,check the unite and direction
%     startLoc = [1 1 1 startT]; count = [inf inf 17 tLin.inter(p_1)]; stride = [1 1 1 1];
% else
%     plev = ncread(temp.name, 'plev') / 100; %
%     startLoc = [1 1 1 startT]; count = [inf inf inf tLin.inter(p_1)]; stride = [1 1 1 1];
% end


        % % 3D vars
        % % find 3D vars exist or not
        % nn = 0;
        % for i1 = 1:length(vars.D3)% 13vars maybe more
        %     temp = dir([filepath_t, vars.D3{i1}, '_Amon*.nc']);
        %     tt = size(temp);
        %     if tt(1) == 0
        %         continue
        %     end
        %     nn = nn + 1;
        %     v3_names{nn} = vars.D3{i1};
        % end

        % % define anomaly\mean\trend data
        % dv3_names = strcat(var_state{1}, v3_names); clim_v3_names = strcat(var_state{2}, v3_names);
        % % trend_v3.mon_name = strcat(var_state{3}, v3_names);trend_v3.sea_name = strcat(var_state{4}, v3_names);
        % % trend_v3.yr_name = strcat(var_state{5}, v3_names);

        % % read 3-D variables
        % startLoc = [1 1 startT]; count = [inf inf tLin.inter(p_1)]; stride = [1 1 1];
        % temp = dir([filepath_t, v3_names{1}, '_Amon*.nc']);
        % lon_v3 = ncread(temp.name, 'lon'); lat_v3 = ncread(temp.name, 'lat');
        % lat_v3 = lat_v3(end:-1:1);
        % temp_v3 = zeros(length(v3_names), length(lon_v3), length(lat_v3), tLin.inter(p_1));

        % for i1 = 1:length(v3_names)% 13vars or 10vars
        %     temp = dir([filepath_t, v3_names{i1}, '_Amon*.nc']);
        %     temp_v3(i1, :, :, :) = ncread(temp.name, v3_names{i1}, startLoc, count, stride);
        % end

        % temp_v3 = temp_v3(:, :, end:-1:1, :); % transfor the lat

        % %% 3) now we finished reading data, next we regrid to 2.5*2.5
        % if max(lon_v3) < 357.5%strcmp(level.model2{level1}(1:end),'CanESM5') == 1
        %     temp_v3(:, end + 1, :, :) = temp_v3(:, 1, :, :);
        %     lon_v3(end + 1) = 360;
        % elseif max(lon_v4) < 357.5%strcmp(level.model2{level1}(1:end),'CanESM5') == 1
        %     temp_v4(:, end + 1, :, :, :) = temp_v4(:, 1, :, :, :);
        %     lon_v4(end + 1) = 360;
        % elseif strcmp(level.model2{level1}(1:end), 'GISS-E2-1-G') == 1
        %     temp_vv4 = temp_v4; temp_vv3 = temp_v3;
        %     temp_vv4(:, end + 1, :, :, :) = 0; temp_vv3(:, end + 1, :, :) = 0;
        %     temp_vv4(:, 1, :, :, :) = temp_v4(:, end, :, :, :);
        %     temp_vv4(:, 2:end, :, :, :) = temp_v4;
        %     temp_vv3(:, 1, :, :) = temp_v3(:, end, :, :);
        %     temp_vv3(:, 2:end, :, :) = temp_v3;
        %     lon_v3 = [-1.25; lon_v3]; lon_v4 = [-1.25; lon_v4];
        %     temp_v3 = temp_vv3; temp_v4 = temp_vv4;
        %     clear temp_vv3 temp_vv4
        % end