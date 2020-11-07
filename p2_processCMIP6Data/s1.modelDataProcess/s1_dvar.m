%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-09 15:52:00
% LastEditTime : 2020-11-05 22:08:14
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
outPath = '/data1/liuyincheng/CMIP6-process/';
var_state = {'d', 'clim_', 'trendm_d', 'trends_d', 'trendyr_d'}; % m,s,yr indicate that month, season, year
kernels_path = '/data1/liuyincheng/y_kernels/YiH/kernels_YiH/toa/dp.nc';
lon_k = 0:2.5:357.5; nlonk = length(lon_k);
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment
for exmNum = 4:4%1 mean amip 2000; 2 mean amip 1980; 3 means ssp245, 4 means ssp370, 6 abrupt-4xCO2_150years
    % model parameters
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % experiment path (tLin:1740)
    inputPath = '/data1/liuyincheng/CMIP6-mirror/';
    exmPath_all = cell(1, length(Experiment));

    for varNum = 1:length(Experiment)
        exmPath_all{varNum} = fullfile(inputPath, Experiment{varNum});
    end

    exmPath = exmPath_all{exmNum};
    time2 = 1:tLin.inter{exmNum};
    plev_k = ncread(kernels_path, 'player'); % pay attention to plev_k's range must smaller than plev's
    nplevk = length(plev_k);
    readme.timeseries = tLin.read{exmNum};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model
    for mdlNum = 6:length(level.model2)% model number
        % model path
        mdlPath = fullfile(exmPath, level.model2{mdlNum});
        eval(['cd ', mdlPath]);
        disp(' ')
        disp([level.model2{mdlNum}, ' model start!'])
        % ensemble member path
        esmName = getPath_fileName(mdlPath, '.');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensemble member
        for esmNum = 1:length(esmName)
            esmPath=fullfile(mdlPath,esmName{esmNum,1});
            eval(['cd ', esmPath]);
            % outputpath and make it
            varsPath = fullfile(outPath, level.time1{exmNum}, level.model2{mdlNum}, esmName{esmNum,1}, level.process3{1}); %'model/ensemble member/rawdata_regrid/'
            dvarsPath = fullfile(outPath, level.time1{exmNum}, level.model2{mdlNum}, esmName{esmNum,1}, level.process3{2}); %'model/ensemble member/anomaly/'
            auto_mkdir(varsPath)
            auto_mkdir(dvarsPath)

            %% check1:  Time line
            % locate first year
            temp = dir(fullfile(esmPath,[vars.D3{1}, '_Amon*.nc']));
            formatOut = 'yyyy-mm';
            startT = cdftime2loc(temp, formatOut, tLin.start{exmNum});
            % read time
            startLoc = startT; count = tLin.inter{exmNum}; stride = 1;
            time = cmipTimeRead(temp.name,startLoc,count,stride); %date, Units, Calendar, length 
            % test time line consistence
            startYear = str2num(tLin.time{exmNum}(1:4)); endYear = str2num(tLin.time{exmNum}(end - 3:end));
            disp('Time line check:')
            testTime(time.date, startYear, endYear, 1)

            %% check2:  plev
            temp = dir(fullfile(esmPath,[vars.D4{1}, '_Amon*.nc']));
            % test plev consistence
            testPlev(temp, mPlev, exmNum)

            %% 2) read data
            % find 4D and 3D vars exist or not
            nn = 0;
            for varNum = 1:length(vars.all)% 13vars maybe more
                temp = dir([esmPath,'/', vars.all{varNum}, '_Amon*.nc']);
                tt = size(temp);
                if tt(1) == 0
                    continue
                end
                nn = nn + 1;
                v_names{nn} = vars.all{varNum};
            end

            % define anomaly\mean\trend data
            dv_names = strcat(var_state{1}, v_names); clim_v_names = strcat(var_state{2}, v_names);
            % one time only cal once
            for varNum = 1:length(v_names)% 13vars maybe more
                % read variables
                temp = dir([esmPath,'/', v_names{varNum}, '_Amon*.nc']);
                lon_v = ncread(temp.name, 'lon'); lat_v = ncread(temp.name, 'lat');
                lat_v = lat_v(end:-1:1);
                disp([v_names{varNum}, '(',level.model2{mdlNum},', ',esmName{esmNum,1},'): '])

                if ismember(1, strcmp(v_names{varNum}, vars.D4))% 4D vars
                    p_2 = 4; % mean 4D vars
                    plev = ncread(temp.name, 'plev') / 100; %
                    startLoc = [1 1 1 startT]; count = [inf inf inf tLin.inter{exmNum}]; stride = [1 1 1 1];
                    temp_v = ncread(temp.name, v_names{varNum}, startLoc, count, stride);
                    temp_v = temp_v(:, end:-1:1, :, :); % transfor the lat
                    if ismember(1, strcmp(v_names{varNum}, 'hus'))
                        temp_v(temp_v < 0) = nan;
                    end
                    % check3 and fix lon/lat to contain scale
                    [pp] = testLonlat(temp);
                    [pp, lon_v, lat_v, temp_v] = fixLonlat(pp, p_2, lon_v, lat_v, temp_v, temp);
                    % regrid to unite axis
                    v_regrid = autoRegrid4(lon_v, lat_v, plev, time2, temp_v, lon_k, lat_k, plev_k, time2);
                    clear temp_v
                elseif ismember(1, strcmp(v_names{varNum}, vars.D3))% 3D vars
                    p_2 = 3; % mean 3D vars
                    startLoc = [1 1 startT]; count = [inf inf tLin.inter{exmNum}]; stride = [1 1 1];
                    temp_v = ncread(temp.name, v_names{varNum}, startLoc, count, stride);
                    temp_v = temp_v(:, end:-1:1, :); % transfor the lat
                    % check3 and fix lon/lat to contain scale
                    [pp] = testLonlat(temp);
                    [pp, lon_v, lat_v, temp_v] = fixLonlat(pp, p_2, lon_v, lat_v, temp_v, temp);
                    % regrid to unite axis
                    v_regrid = autoRegrid3(lon_v, lat_v, time2, temp_v, lon_k, lat_k, time2);
                    clear temp_v
                end

                % save and release space
                eval([v_names{varNum}, '= v_regrid;']); %
                save([varsPath, v_names{varNum}], v_names{varNum});
                eval(['clear ', v_names{varNum}])
                %% now we finished, next we cal the anomaly.
                [v_anom, v_clim] = monthlyAnomaly(144, 73, nplevk, time.date, v_regrid, 1); %[nlongitude,nlatitude,time,var,startmonth]
                clear v_regrid
                % save and release space
                eval([dv_names{varNum}, '= v_anom;']); %
                eval([clim_v_names{varNum}, '= v_clim;']); %
                save([dvarsPath, dv_names{varNum}], dv_names{varNum}, clim_v_names{varNum});
                eval(['clear ', dv_names{varNum}, ' ', clim_v_names{varNum}])
                clear v_anom v_clim

            end

            % save grobal vars
            timeseries = tLin.time{exmNum};
            modelname = level.model2{mdlNum}(1:end);
            save([varsPath, 'global_vars.mat'], 'lon_k', 'lat_k', 'time', 'plev_k', 'readme', 'timeseries', 'modelname')
            save([dvarsPath, 'global_vars.mat'], 'lon_k', 'lat_k', 'time', 'plev_k', 'readme', 'timeseries', 'modelname')
            
            disp([esmName{esmNum,1}, ' ensemble is done!'])
        end
        disp([level.model2{mdlNum}, ' model is done!'])
        disp(' ')
    end
    disp([level.time1{exmNum}, ' era is done!'])
    disp(' ')
end

eval(['cd ', nowpath]);

t = toc; tickTok(t)
