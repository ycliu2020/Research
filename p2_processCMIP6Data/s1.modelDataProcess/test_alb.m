%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% monthly data
% cal the anomly
% hope to use an unite scrip to read and process all the varies we need
% Compare the pre rad and ts data from 2000-03 to 2018-02(18 years)
%
% time:2000.01-2014.12(interval:15*12);1980.01-2014.12(interval:35*12); 2015.01-2099.12(interval:85*12)
% initial time in amip(432 total): 253 of 432(2000.03);13 of 432(1980.01);
% initial time in futrue(1032 total): 1 of 1032(2015.01);
% initial time in amip-hist(1740 total): 1,561 of 1740(2000.03);1,321 of 1740(1980.01);
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Attention!!!
% albeo mean surface albeo
% check lat: model lat disagree with kernels lat (Opposite direction)
% check plev drection and unite: model plev disagree with kernels plev
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear; clc; tic;

for p_1 = 1:1%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
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
    alb_names = {'alb'};
    var_state = {'d', 'clim_', 'trendm_d', 'trends_d', 'trendyr_d'}; % m,s,yr indicate that month, season, year
    kernels_path = '/data/pub/kernels_YiH/toa/dp.nc';
    time2 = 1:tLin.inter{p_1};
    plevf = ncread(kernels_path, 'player'); % pay attention to plevf's range must smaller than plev's
    plevfnum = length(plevf);
    lonf = 0:2.5:357.5; nlonf = length(lonf);
    latf = 90:-2.5:-90; nlatf = length(latf);
    readme.timeseries = tLin.read{p_1};
    readme.sfcAlbeo = 'sfc albeo';

    for level1 = 3:3% BCC-CSM2-MR
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
        % test time line consistence
        startYear = str2num(tLin.time{p_1}(1:4)); endYear = str2num(tLin.time{p_1}(end - 3:end));
        ntime = length(time);
        disp('Time line check:')
        testTime(time, startYear, endYear, 1)

        %% check plev
        temp = dir([filepath_t, vars.D4{1}, '_Amon*.nc']);
        % test plev consistence
        testPlev(temp, mPlev, p_1)

        %% 2) read data
        % find only rsds and rsus.4D and 3D vars exist or not
        nn = 0;
        var1Num = find(strcmp(vars.all, 'rsds'));
        var2Num = find(strcmp(vars.all, 'rsus'));

        for ii = [var1Num var2Num]% only 'rsds', 'rsus',
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
        dalb_names = strcat(var_state{1}, alb_names); clim_alb_names = strcat(var_state{2}, alb_names);

        %% read 'rsds' and 'rsus',
        startLoc = [1 1 startT]; count = [inf inf tLin.inter{p_1}]; stride = [1 1 1];
        temp = dir([filepath_t, v_names{1}, '_Amon*.nc']);
        lon_v = ncread(temp.name, 'lon'); lat_v = ncread(temp.name, 'lat');
        lat_v = lat_v(end:-1:1);
        temp_v3 = zeros(length(v_names), length(lon_v), length(lat_v), tLin.inter{p_1});

        for i1 = 1:length(v_names)%1.'rsds', 2.'rsus',
            temp = dir([filepath_t, v_names{i1}, '_Amon*.nc']);
            temp_v3(i1, :, :, :) = ncread(temp.name, v_names{i1}, startLoc, count, stride);
        end

        temp_v3 = temp_v3(:, :, end:-1:1, :); % transfor the lat

        %% check and fix cal alb 
        %1.'rsds', 2.'rsus',
        if ~isempty(find(temp_v3(2, :, :, :) < 0, 1))
            disp([level.model2{level1}, ' model var rsus<0, now replace <0 with 0;'])
        end

        if ~isempty(find(temp_v3(1, :, :, :) < 0, 1))
            disp([level.model2{level1}, ' model var rsds<0, now replace <0 with 0!!!'])
        end

        fix_1 = squeeze(temp_v3(1, :, :, :));
        fix_1(fix_1 < 0) = 0;
        temp_v3(1, :, :, :) = fix_1;
        fix_2 = squeeze(temp_v3(2, :, :, :));
        fix_2(fix_2 < 0) = 0;
        temp_v3(2, :, :, :) = fix_2;

        temp_v3(temp_v3 == 0) = eps(0);
        temp_alb = squeeze(temp_v3(2, :, :, :) ./ temp_v3(1, :, :, :));
        temp_alb(temp_v3(1, :, :, :) == eps(0)) = 0;

        if ~isempty(find(temp_alb > 1, 1))
            disp([level.model2{level1}, ' model have some area up greater than down radiation, now replace alb >1 or inf with 1;'])
        end
        fix_3 = temp_alb;
        fix_3(fix_3 > 1) = 1;
        temp_alb = fix_3;

        %% check and fix lon/lat to contain scale
        [pp] = testLonlat(temp);
        [pp, lon_v, lat_v, temp_alb] = fixLonlat(pp, p_2, lon_v, lat_v, temp_alb, temp);

        [Xlon, Ylat, Ttime] = meshgrid(lat_v, lon_v, time2);
        [Xlonf, Ylatf, Ttimef] = meshgrid(latf, lonf, time2);
        albeo_regrid = interp3(Xlon, Ylat, Ttime, temp_alb, Xlonf, Ylatf, Ttimef);

        % now we finished, next we cal the anomaly.
        [alb_anom, alb_clim] = monthlyAnomaly3D(144, 73, time, albeo_regrid, 1);

        %% save sfc albeo
        eval([alb_names{1}, '= albeo_regrid;']); %
        save([outPathName{1}, alb_names{1}], alb_names{1});
        eval([dalb_names{1}, '= alb_anom;']); %
        eval([clim_alb_names{1}, '= alb_clim;']); %
        save([outPathName{2}, dalb_names{1}], dalb_names{1}, clim_alb_names{1});
        % % save grobal vars
        % timeseries = tLin.time{p_1};
        % modelname = level.model2{level1}(1:end);
        % save([outPathName{1}, 'global_vars.mat'], 'lonf', 'latf', 'time', 'plevf', 'readme', 'timeseries', 'modelname')
        % save([outPathName{2}, 'global_vars.mat'], 'lonf', 'latf', 'time', 'plevf', 'readme', 'timeseries', 'modelname')
        % disp([level.model2{level1}, ' model is done!'])
        % disp(' ')

    end

    % disp([level.time1{p_1}, ' era is done!'])
    % disp(' ')
end

% eval(['cd ', nowpath]);

t = toc; disp(t)
