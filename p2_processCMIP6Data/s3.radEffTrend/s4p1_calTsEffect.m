%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-15 10:44:34
% LastEditTime : 2020-06-15 11:23:40
% LastEditors  : LYC
% Description  : cal dRvars effect on Ts : according Ts=dRx/Kts
% FilePath     : /Research/p2_processCMIP6Data/s3.radEffTrend/s4p1_calTsEffect.m
%  
%%---------------------------------------------------------
clear; clc; tic;
p_3 = 88.75; % Latitude range
lon1 = [2.5 357.5]; lat1 = [-p_3 + 1 p_3 - 1]; % world area
toaSfc = {'toa', 'sfc'};

for p_1 = 1:4%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    nowpath = pwd;
    [readme, Experiment, level, tLin, mPlev, vars] = modelParameters(p_1);
    % inputPath
    inputPath = ['/data1/liuyincheng/cmip6-process/', level.time1{p_1}]; %/data1/liuyincheng/cmip6-process/amip_2000-2014/
    dvarsPath = [inputPath, level.model2{1}, '/', level.process3{6}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
    load([dvarsPath, 'global_vars.mat'])% lat lon time plevf readme
    nlat = length(lat); nlon = length(lon); ntime = length(time);


    for level1 = 1:length(level.model2)


        disp([level.model2{level1}, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{p_1}, ' era is done!'])
    disp(' ')

end

t = toc; disp(t)
