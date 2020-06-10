%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-09 15:52:00
% LastEditTime : 2020-06-10 15:16:50
% LastEditors  : LYC
% Description  : 
% FilePath     : /Research/p2_processCMIP6Data/s3.radEffTrend/s1_calPsmonthMean.m
%  
%%---------------------------------------------------------
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% mothly data
% cal month mean ps and thickness dps, dp
% this script mainly include:ps and dps, dp
% PS:
% ps_m(144x73x12) denotes the month mean ps
% dps denotes the thickness (in hPa) of surface layer(month mean): dps((144x73x12))
% dp denotes the thickness (in hPa) of each layer(month mean): dp((144x73x12))
%
% hope to use an unite scrip to read and process all the varies we need
% Compare the data from 2000-03 to 2018-02(18 years)
%
% experiment information:
% time:2000.01-2014.12(interval:15*12);1980.01-2014.12(interval:35*12); 2015.01-2099.12(interval:85*12)
% initial time in hist(1740 total): 1,561 of 1740(2000.03);1,321 of 1740(1980.01)
% initial time in futrue(1032 total): 1 of 1032(2015.01);
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear; clc; tic;
nowpath = pwd;
% model settings
for p_1 = 1:4%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    [readme, Experiment, level, tLin, mPlev, vars] = modelParameters(p_1);
    % inputPath
    inputPath = ['/data1/liuyincheng/cmip6-process/', level.time1{p_1}]; %~/data/cmip6/2000-2014/

    % read kernels data
    kernelPaths = '/data1/liuyincheng/y_kernels/kernels_YiH/surface/';
    plevel = ncread([kernelPaths, 'dp.nc'], 'plevel');
    player = ncread([kernelPaths, 'dp.nc'], 'player');
    dp_raw = ncread([kernelPaths, 'dp.nc'], 'dp');

    for level1 = 1:length(level.model2)% model numbers
        % Part1 month mean ps (ps_m(144x73x12))
        inputPath1 = [inputPath, level.model2{level1}, '/', level.process3{1}]; %~/data/cmip6/2000-2014/MRI-ESM2-0/rawdata_regrids
        load([inputPath1, 'global_vars.mat'])% latf lonf time plevf readme
        load([inputPath1, 'ps.mat'])
        ntime = length(time); nlatf = length(latf); nlonf = length(lonf);
        ps = ps ./ 100; % transform hpa unite
        ps_m = zeros(nlonf, nlatf, 12);
        [yy, mm, dd] = datevec(time);

        for im = 1:12
            index_mm = mm == im;
            ps_m(:, :, im) = nanmean(ps(:, :, index_mm), 3); %this is the average of each month
        end

        month = 1:12;
        % Part2 cal dps,dp
        dp = zeros(nlonf, nlatf, 24, 12); % revise dp in all layers, under the lowest layer are zeros
        dp_level2 = zeros(nlonf, nlatf, 24, 12); % only contain the near sfc
        dps = zeros(nlonf, nlatf, 12);

        for i = 1:nlonf

            for j = 1:nlatf

                for nt = 1:12

                    if isnan(ps_m(i, j, nt)) == 1
                        continue
                    end

                    temp = find(player < ps_m(i, j, nt), 1, 'first');
                    dps(i, j, nt) = ps_m(i, j, nt) - plevel(temp + 1);
                    % dps(i, j, nt) = dp_bottom(i,j)123;
                    dp(i, j, temp, nt) = dps(i, j, nt);
                    dp(i, j, temp + 1:24, nt) = dp_raw(temp + 1:24);
                    % near surface level2
                    dp_level2(i, j, temp, nt) = dps(i, j, nt);
                    dp_level2(i, j, temp + 1, nt) = dp_raw(temp + 1);
                end

            end

        end

        kernelsOutPath = [inputPath, level.model2{level1}, '/', level.process3{5}]; % /home/lyc/data/cmip6/2000-2014/MRI-ESM2-0/kernelsCal/
        auto_mkdir(kernelsOutPath)
        save([kernelsOutPath, 'global_vars.mat'], 'lonf', 'latf', 'time', 'plevf', 'readme', 'timeseries', 'modelname')
        save([kernelsOutPath, 'kernel_dp.mat'], 'dps', 'dp', 'dp_level2');
        save([kernelsOutPath, 'kernel_ps_m.mat'], 'ps_m');
        disp([level.model2{level1}, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{p_1}, ' era is done!'])
    disp(' ')
end

t = toc; disp(t)
