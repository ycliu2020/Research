%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-14 19:53:54
% LastEditTime : 2020-07-11 10:41:13
% LastEditors  : LYC
% Description  : 计算非局地云的辐射效应和温度贡献(用lamd_c和dtsg的线性关系求)
% FilePath     : /code/p2_processCMIP6Data/s3.nonLocalCld/s1p2_cal1_nonLocalCld.m
% Attention    : only ssp caled
%%---------------------------------------------------------
clear; clc; tic;
p_3 = 88.75; % Latitude range
lon1 = [2.5 357.5]; lat1 = [-p_3 + 1 p_3 - 1]; % world area
toaSfc = {'toa', 'sfc'};
% load lamda_cloud
load('/data1/liuyincheng/cmip6-process/z_globalVar/lamda_cloud.mat')

for p_1 = 3:4%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    nowpath = pwd;
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(p_1);
    % inputPath
    inputPath = ['/data1/liuyincheng/cmip6-process/', level.time1{p_1}]; %/data1/liuyincheng/cmip6-process/amip_2000-2014/
    dradEffectPath = [inputPath, level.model2{1}, '/', level.process3{6}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
    load([dradEffectPath, 'global_vars.mat'])% lat_f lon_f time plev_k readme
    nlatf = length(lat_f); nlonf = length(lon_f); ntime = length(time.date);

    % load k1_cld, k2_cld, existModelName
    load([inputPath, 'z_assembleData/non_localCld/k_cld.mat'])

    for level1 = 1:length(existModelName)
        % load data
        dvarsPath = [inputPath, existModelName{level1}, '/', level.process3{2}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
        dvarsTrendPath = [inputPath, existModelName{level1}, '/', level.process3{3}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly_trend
        kernelPath = [inputPath, existModelName{level1}, '/', level.process3{5}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/kernelsCal
        varsPath = [inputPath, existModelName{level1}, '/', level.process3{1}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/rawdata
        dradEffectPath = [inputPath, existModelName{level1}, '/', level.process3{6}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
        dnonLocalCldPath = [inputPath, existModelName{level1}, '/', level.process3{8}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/non_localCld/
        % load dtsg
        load([dvarsPath, 'dtsg_nomask.mat'])

        lamda_c = lamda_cloud.modelVar(strcmp(lamda_cloud.modelName, existModelName{level1}) == 1);
        varUsed = zeros(nlonf, nlatf, ntime, 3); % dR
        varKerUsed=varUsed; % dR/kernels
        for skyLevel = 1:2% 1 mean toa, 2 mean sfc
            load([dradEffectPath, ['dradEfect_', toaSfc{skyLevel}, '_cld.mat']])%albEffect, husEffect, taEffect, mainEffect, totalEffect, tsEffect, talwEffect, taswEffect
            varUsed(:, :, :, 1) = husEffect;
            varUsed(:, :, :, 2) = taEffect;
            varUsed(:, :, :, 3) = albEffect;
            varParitial = zeros(nlonf, nlatf, 3);
            varKerParitial=varParitial;
            %% Part1: cal trendyr_dRnonLocalCld_toa/sfc and trendyr_dRk2_toa/sfc(unite: per day)
            for varNum = 1:3
                tempTemp = squeeze(varUsed(:, :, :, varNum));

                for latNum = 1:nlatf

                    for lonNum = 1:nlonf
                        tempUsed = squeeze(squeeze(tempTemp(lonNum, latNum, :)));
                        k = polyfit(dtsg, tempUsed, 1); % 一元一次拟合
                        varParitial(lonNum, latNum, varNum) = k(1);
                    end

                end

            end

            dRnonLocalCld_hus = squeeze(varParitial(:, :, 1)) * k1_cld * lamda_c;
            dRnonLocalCld_ta = squeeze(varParitial(:, :, 2)) * k1_cld * lamda_c;
            dRnonLocalCld_alb = squeeze(varParitial(:, :, 3)) * k1_cld * lamda_c;
            dRk2_hus = squeeze(varParitial(:, :, 1)) * k2_cld;
            dRk2_ta = squeeze(varParitial(:, :, 2)) * k2_cld;
            dRk2_alb = squeeze(varParitial(:, :, 3)) * k2_cld;
            dRnonLocalCld = dRnonLocalCld_hus + dRnonLocalCld_ta + dRnonLocalCld_alb;
            dRk2 = dRk2_hus + dRk2_ta + dRk2_alb;

            trendyr_dRnonLocalCld = dRnonLocalCld ./ (time.date(end) - time.date(1));
            trendyr_dRk2 = dRk2 ./ (tim.datee(end) - time.date(1));
            save([dnonLocalCldPath, 'dRnonLocalCld_', toaSfc{skyLevel}, '.mat'], 'dRnonLocalCld_hus', 'dRnonLocalCld_ta', 'dRnonLocalCld_alb', 'dRnonLocalCld')
            save([dnonLocalCldPath, 'dRk2_', toaSfc{skyLevel}, '.mat'], 'dRk2_hus', 'dRk2_ta', 'dRk2_alb', 'dRk2')
            save([dnonLocalCldPath, 'trendyr_dRnonLocalCld_', toaSfc{skyLevel}, '.mat'], 'trendyr_dRnonLocalCld')
            save([dnonLocalCldPath, 'trendyr_dRk2_', toaSfc{skyLevel}, '.mat'], 'trendyr_dRk2')
            
            %% Part2: cal trendyr_dTsnonLocalCld_toa/sfc and trendyr_dTsk2_toa/sfc(unite: per day)
            % load ts kernel
            load([kernelPath, '/kernels_', toaSfc{skyLevel}, '_cld'], 'ts_lwkernel');
            load([kernelPath, '/global_vars.mat']); % lat_k,lon_k,time
            % regrid ts_lwkernel to 144x72(unite grids)
            kernelTime = 1:12;
            ts_lwkernel = autoRegrid3(lon_k, lat_k, kernelTime, ts_lwkernel, lon_f, lat_f, kernelTime);
            % extend to the whole time series
            startmonth = 1;
            if startmonth ~= 1
                tempK=zeros(nlonf,nlatf,12);
                tempK(1:12-startmonth+1) = ts_lwkernel(:,:,startmonth:12);
                tempK(12-startmonth+2:12) =ts_lwkernel(:,:,1:startmonth-1);
                ts_lwkernel=tempK;
            end
            nyear=ntime/12;
            ts_lwkernelSeries=repmat(ts_lwkernel,[1 1 nyear]);
            
            for varNum = 1:3
                varKerUsed(:,:,:,varNum)=squeeze(varUsed(:,:,:,varNum))./-ts_lwkernelSeries;% note - didnt run 
            end
            
            for varNum = 1:3
                tempTemp = squeeze(varKerUsed(:, :, :, varNum));
                for latNum = 1:nlatf
                    for lonNum = 1:nlonf
                        tempUsed = squeeze(squeeze(tempTemp(lonNum, latNum, :)));
                        k = polyfit(dtsg, tempUsed, 1); % 一元一次拟合
                        varKerParitial(lonNum, latNum, varNum) = k(1);
                    end
                end
            end
            dTsnonLocalCld_hus = squeeze(varKerParitial(:, :, 1)) * k1_cld * lamda_c;
            dTsnonLocalCld_ta = squeeze(varKerParitial(:, :, 2)) * k1_cld * lamda_c;
            dTsnonLocalCld_alb = squeeze(varKerParitial(:, :, 3)) * k1_cld * lamda_c;
            dTsk2_hus = squeeze(varKerParitial(:, :, 1)) * k2_cld;
            dTsk2_ta = squeeze(varKerParitial(:, :, 2)) * k2_cld;
            dTsk2_alb = squeeze(varKerParitial(:, :, 3)) * k2_cld;
            dTsnonLocalCld = dTsnonLocalCld_hus + dTsnonLocalCld_ta + dTsnonLocalCld_alb;
            dTsk2 = dTsk2_hus + dTsk2_ta + dTsk2_alb;
            
            trendyr_dTsnonLocalCld = dTsnonLocalCld ./ (time.date(end) - time.date(1));
            trendyr_dTsk2 = dTsk2 ./ (time.date(end) - time.date(1));
            save([dnonLocalCldPath, 'dTsnonLocalCld_', toaSfc{skyLevel}, '.mat'], 'dTsnonLocalCld_hus', 'dTsnonLocalCld_ta', 'dTsnonLocalCld_alb', 'dTsnonLocalCld')
            save([dnonLocalCldPath, 'dTsk2_', toaSfc{skyLevel}, '.mat'], 'dTsk2_hus', 'dTsk2_ta', 'dTsk2_alb', 'dTsk2')
            save([dnonLocalCldPath, 'trendyr_dTsnonLocalCld_', toaSfc{skyLevel}, '.mat'], 'trendyr_dTsnonLocalCld')
            save([dnonLocalCldPath, 'trendyr_dTsk2_', toaSfc{skyLevel}, '.mat'], 'trendyr_dTsk2')

        end

        disp([existModelName{level1}, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{p_1}, ' era is done!'])
    disp(' ')

end

t = toc; disp(t)
