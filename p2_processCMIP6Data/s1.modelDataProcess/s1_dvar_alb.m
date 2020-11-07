%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-09 15:52:00
% LastEditTime : 2020-11-03 17:14:39
% LastEditors  : LYC
% Description  : cal mainly alb include 1.regrid vars, 2.vars anomly
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
outPath = '/data1/liuyincheng/CMIP6-process/';
alb_names = {'alb'};
var_state = {'d', 'clim_', 'trendm_d', 'trends_d', 'trendyr_d'}; % m,s,yr indicate that month, season, year
kernels_path = '/data1/liuyincheng/y_kernels/YiH/kernels_YiH/toa/dp.nc';
lon_k = 0:2.5:357.5; nlonk = length(lon_k);
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment
for exmNum = 1:1%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    % model parameters
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % input and output path (tLin:1740)
    inputPath = '/data1/liuyincheng/CMIP6-mirror/';
    exmPath_all = cell(1, length(Experiment));

    for ii = 1:length(Experiment)
        exmPath_all{ii} = [inputPath, Experiment{ii}, '/'];
    end

    exmPath = exmPath_all{exmNum};
    time2 = 1:tLin.inter{exmNum};
    plev_k = ncread(kernels_path, 'player'); % pay attention to plev_k's range must smaller than plev's
    nplevf = length(plev_k);
    readme.timeseries = tLin.read{exmNum};
    readme.sfcAlbeo = 'sfc albeo';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model
    for mdlNum = 1:length(level.model2)% model number
        % inputpath
        mdlPath = fullfile(exmPath, level.model2{mdlNum});
        eval(['cd ', mdlPath]);
        disp(' ')
        disp([level.model2{mdlNum}, ' model start!'])
        % ensemble member path
        esmName = getPath_fileName(mdlPath, '.');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensemble member
        for esmNum = 1:length(esmName)
            esmPath = fullfile(mdlPath, esmName{esmNum, 1});
            eval(['cd ', esmPath]);
            % outputpath and make it
            outPathName{1} = fullfile(outPath, level.time1{exmNum}, level.model2{mdlNum}, esmName{esmNum,1}, level.process3{1}); %'model/ensemble member/rawdata_regrid/'
            outPathName{2} = fullfile(outPath, level.time1{exmNum}, level.model2{mdlNum}, esmName{esmNum,1}, level.process3{2}); %'model/ensemble member/anomaly/'
            auto_mkdir(outPathName{1})
            auto_mkdir(outPathName{2})

            %% check1:  Time line
            % locate first year
            temp = dir(fullfile(esmPath,[vars.D3{1}, '_Amon*.nc']));
            formatOut = 'yyyy-mm';
            startT = cdftime2loc(temp, formatOut, tLin.start{exmNum});
            % read time
            startLoc = startT; count = tLin.inter{exmNum}; stride = 1;
            time = cmipTimeRead(temp.name,startLoc,count,stride); %date, Units, Calendar, length 
            % test time line consistence
            startYear = str2double(tLin.time{exmNum}(1:4)); endYear = str2double(tLin.time{exmNum}(end - 3:end));
            disp('Time line check:')
            testTime(time.date, startYear, endYear, 1)

            %% check2:  plev
            temp = dir(fullfile(esmPath,[vars.D4{1}, '_Amon*.nc']));
            % test plev consistence
            testPlev(temp, mPlev, exmNum)

            %% 2) read data
            % find only rsds and rsus.4D and 3D vars exist or not
            nn = 0;
            var1Num = find(strcmp(vars.all, 'rsdscs'));
            var2Num = find(strcmp(vars.all, 'rsuscs'));

            for ii = [var1Num var2Num]  % only 'rsds', 'rsus',
                temp = dir([esmPath,'/', vars.all{ii}, '_Amon*.nc']);
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
            startLoc = [1 1 startT]; count = [inf inf tLin.inter{exmNum}]; stride = [1 1 1];
            temp = dir([esmPath,'/', v_names{1}, '_Amon*.nc']);
            lon_v = ncread(temp.name, 'lon'); lat_v = ncread(temp.name, 'lat');
            lat_v = lat_v(end:-1:1);
            temp_v3 = zeros(length(v_names), length(lon_v), length(lat_v), tLin.inter{exmNum});

            for i1 = 1:length(v_names) % 1.'rsds', 2.'rsus',
                temp = dir([esmPath,'/', v_names{i1}, '_Amon*.nc']);
                temp_v3(i1, :, :, :) = ncread(temp.name, v_names{i1}, startLoc, count, stride);
            end

            temp_v3 = temp_v3(:, :, end:-1:1, :); % transfor the lat

            %% check and fix cal alb
            %1.'rsds', 2.'rsus',
            if ~isempty(find(temp_v3(2, :, :, :) < 0, 1))
                disp([level.model2{mdlNum}, ' model var rsus<0, now replace <0 with 0;'])
            end

            if ~isempty(find(temp_v3(1, :, :, :) < 0, 1))
                disp([level.model2{mdlNum}, ' model var rsds<0, now replace <0 with 0!!!'])
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
                disp([level.model2{mdlNum}, ' model have some area up greater than down radiation, now replace alb >1 or inf with 1;'])
            end

            fix_3 = temp_alb;
            fix_3(fix_3 > 1) = 1;
            temp_alb = fix_3;

            p_2 = 3; % 3D vars
            %% check and fix lon/lat to contain scale
            [pp] = testLonlat(temp);
            [pp, lon_v, lat_v, temp_alb] = fixLonlat(pp, p_2, lon_v, lat_v, temp_alb, temp);

            albeo_regrid = autoRegrid3(lon_v, lat_v, time2, temp_alb, lon_k, lat_k, time2);
            % now we finished, next we cal the anomaly.
            [alb_anom, alb_clim] = monthlyAnomaly3D(144, 73, time.date, albeo_regrid, 1);

            %% save sfc albeo
            eval([alb_names{1}, '= albeo_regrid;']); %
            save([outPathName{1}, alb_names{1}], alb_names{1});
            eval([dalb_names{1}, '= alb_anom;']); %
            eval([clim_alb_names{1}, '= alb_clim;']); %
            save([outPathName{2}, dalb_names{1}], dalb_names{1}, clim_alb_names{1});
            % save grobal vars
            timeseries = tLin.time{exmNum};
            modelname = level.model2{mdlNum}(1:end);
            save([outPathName{1}, 'global_vars.mat'], 'lon_k', 'lat_k', 'time', 'plev_k', 'readme', 'timeseries', 'modelname')
            save([outPathName{2}, 'global_vars.mat'], 'lon_k', 'lat_k', 'time', 'plev_k', 'readme', 'timeseries', 'modelname')
            
            disp([esmName{esmNum,1}, ' ensemble is done!'])
        end

        disp([level.model2{mdlNum}, ' model is done!'])
        disp(' ')

    end

    disp([level.time1{exmNum}, ' era is done!'])
    disp(' ')
end

eval(['cd ', nowpath]);

t = toc; disp(t)
