%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-09 15:52:00
% LastEditTime : 2020-06-11 21:23:23
% LastEditors  : LYC
% Description  : 
% FilePath     : /Research/p2_processCMIP6Data/s1.modelDataProcess/dvar.m
%  
%%---------------------------------------------------------
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
for p_1 = 3:4%1 mean amip 2000; 2 mean amip 1980; 4 means ssp245, 5 means ssp370, 6 abrupt-4xCO2_150years
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
    kernels_path = '/data1/liuyincheng/y_kernels/kernels_YiH/toa/dp.nc';
    time2 = 1:tLin.inter{p_1};
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

        %% check Time line
        % locate first year
        temp = dir([filepath_t, vars.D3{1}, '_Amon*.nc']);
        formatOut = 'yyyy-mm';
        startT = cdftime2loc(temp, formatOut, tLin.start{p_1});
        % read time
        startLoc = startT; count = tLin.inter{p_1}; stride = 1;
        time = ncread(temp.name, 'time', startLoc, count, stride); %unite:day
        timeUnits = ncreadatt(temp.name, 'time', 'units');
        timeCalendar = ncreadatt(temp.name, 'time', 'calendar');
        time = cdfdate2num(timeUnits, timeCalendar, time);
        ntime = length(time);
        % test time line consistence
        startYear = str2num(tLin.time{p_1}(1:4)); endYear = str2num(tLin.time{p_1}(end - 3:end));
        disp('Time line check:')
        testTime(time, startYear, endYear, 1)

        %% check plev
        temp = dir([filepath_t, vars.D4{1}, '_Amon*.nc']);
        % test plev consistence
        testPlev(temp, mPlev, p_1)

        %% 2) read data
        % find 4D and 3D vars exist or not
        nn = 0;

        for ii = 1:length(vars.all)% 13vars maybe more
            temp = dir([filepath_t, vars.all{ii}, '_Amon*.nc']);
            tt = size(temp);

            if tt(1) == 0
                continue
            end

            nn = nn + 1;
            v_names{nn} = vars.all{ii};
        end

        % define anomaly\mean\trend data
        dv_names = strcat(var_state{1}, v_names); clim_v_names = strcat(var_state{2}, v_names);
        % one time only cal once
        for ii = 1:length(v_names)% 13vars maybe more
            % read variables
            temp = dir([filepath_t, v_names{ii}, '_Amon*.nc']);
            lon_v = ncread(temp.name, 'lon'); lat_v = ncread(temp.name, 'lat');
            lat_v = lat_v(end:-1:1);
            disp([v_names{ii}, ': '])

            if ismember(1, strcmp(v_names{ii}, vars.D4))% 4D vars
                p_2 = 4; % mean 4D vars
                plev = ncread(temp.name, 'plev') / 100; %
                startLoc = [1 1 1 startT]; count = [inf inf inf tLin.inter{p_1}]; stride = [1 1 1 1];
                temp_v = ncread(temp.name, v_names{ii}, startLoc, count, stride);
                temp_v = temp_v(:, end:-1:1, :, :); % transfor the lat

                if ismember(1, strcmp(v_names{ii}, 'hus'))
                    temp_v(temp_v < 0) = nan;
                end

                % check and fix lon/lat to contain scale
                [pp] = testLonlat(temp);
                [pp, lon_v, lat_v, temp_v] = fixLonlat(pp, p_2, lon_v, lat_v, temp_v, temp);
                % regrid to unite axis
                [Xlonf, Ylatf, Zplevf, Ttimef] = ndgrid(lonf, latf, plevf, time2);
                [Xlon, Ylat, Zplev, Ttime] = ndgrid(lon_v, lat_v, plev, time2);
                v_regrid = interpn(Xlon, Ylat, Zplev, Ttime, temp_v, Xlonf, Ylatf, Zplevf, Ttimef);
                clear temp_v
            elseif ismember(1, strcmp(v_names{ii}, vars.D3))% 3D vars
                p_2 = 3; % mean 3D vars
                startLoc = [1 1 startT]; count = [inf inf tLin.inter{p_1}]; stride = [1 1 1];
                temp_v = ncread(temp.name, v_names{ii}, startLoc, count, stride);
                temp_v = temp_v(:, end:-1:1, :); % transfor the lat
                % check and fix lon/lat to contain scale
                [pp] = testLonlat(temp);
                [pp, lon_v, lat_v, temp_v] = fixLonlat(pp, p_2, lon_v, lat_v, temp_v, temp);
                % regrid to unite axis
                [Xlon, Ylat, Ttime] = meshgrid(lat_v, lon_v, time2);
                [Xlonf, Ylatf, Ttimef] = meshgrid(latf, lonf, time2);
                v_regrid = interp3(Xlon, Ylat, Ttime, temp_v, Xlonf, Ylatf, Ttimef);
                clear temp_v
            end

            % save and release space
            eval([v_names{ii}, '= v_regrid;']); %
            save([outPathName{1}, v_names{ii}], v_names{ii});
            eval(['clear ', v_names{ii}])
            %% now we finished, next we cal the anomaly.
            [v_anom, v_clim] = monthlyAnomaly(144, 73, plevfnum, time, v_regrid, 1); %[nlongitude,nlatitude,time,var,startmonth]
            clear v_regrid
            % save and release space
            eval([dv_names{ii}, '= v_anom;']); %
            eval([clim_v_names{ii}, '= v_clim;']); %
            save([outPathName{2}, dv_names{ii}], dv_names{ii}, clim_v_names{ii});
            eval(['clear ', dv_names{ii}, ' ', clim_v_names{ii}])
            clear v_anom v_clim

        end

        % save grobal vars
        timeseries = tLin.time{p_1};
        modelname = level.model2{level1}(1:end);
        save([outPathName{1}, 'global_vars.mat'], 'lonf', 'latf', 'time', 'plevf', 'readme', 'timeseries', 'modelname')
        save([outPathName{2}, 'global_vars.mat'], 'lonf', 'latf', 'time', 'plevf', 'readme', 'timeseries', 'modelname')
        disp([level.model2{level1}, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{p_1}, ' era is done!'])
    disp(' ')
end

eval(['cd ', nowpath]);

t = toc; disp(t)
