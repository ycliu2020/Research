%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-08-09 14:58:33
% LastEditTime : 2021-04-02 10:59:08
% LastEditors  : Please set LastEditors
% Description  : v2 和 v1 最大的不同就在于v2采用了散点累积的积分, 而v1采用了线性关系的假设.
% FilePath     : /code/p2_processCMIP6Data/s4.nonLocalCld/s2_nonLocalCld_method3_v2.m
%
%%---------------------------------------------------------

clear; clc; tic;
latRange = 88.75; % Latitude range
lon1 = [2.5 357.5]; lat1 = [-latRange + 1 latRange - 1]; % world area
toaSfc = {'toa', 'sfc'};
lon_k = 0:2.5:357.5; nlonk = length(lon_k); % kernel lat lon
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);
load('/data1/liuyincheng/CMIP6-process/z_globalVar/ERF_rec.mat')% 'ERF_rec', 'timeERF_rec'

for exmNum = 3:4%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    nowpath = pwd;
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % exmPath
    exmPath = ['/data1/liuyincheng/CMIP6-process/', level.time1{exmNum}]; %/data1/liuyincheng/CMIP6-process/amip_2000-2014/
    %%% load ERF data
    load('/data1/liuyincheng/CMIP6-process/z_globalVar/ERF_rec.mat')% 'ERF_rec', 'timeERF_rec'

    if strcmp(Experiment{exmNum}(1:3), 'ami')
        timeERF = timeERF_rec.hist(find(timeERF_rec.hist == tLin.startYear{exmNum}):find(timeERF_rec.hist == tLin.endYear{exmNum}));
        ERForcing = ERF_rec.hist(find(timeERF_rec.hist == tLin.startYear{exmNum}):find(timeERF_rec.hist == tLin.endYear{exmNum}));
    elseif strcmp(Experiment{exmNum}(1:6), 'ssp245')
        timeERF = timeERF_rec.ssp(find(timeERF_rec.ssp == tLin.startYear{exmNum}):find(timeERF_rec.ssp == tLin.endYear{exmNum}));
        ERForcing = ERF_rec.ssp245(find(timeERF_rec.ssp == tLin.startYear{exmNum}):find(timeERF_rec.ssp == tLin.endYear{exmNum}));
    elseif strcmp(Experiment{exmNum}(1:6), 'ssp370')
        timeERF = timeERF_rec.ssp(find(timeERF_rec.ssp == tLin.startYear{exmNum}):find(timeERF_rec.ssp == tLin.endYear{exmNum}));
        ERForcing = ERF_rec.ssp370(find(timeERF_rec.ssp == tLin.startYear{exmNum}):find(timeERF_rec.ssp == tLin.endYear{exmNum}));
    else
        ERForcing = 0;
        disp('this experient doesnt need ERF or havent input ERF Data!')
    end

    % if length(ERForcing) ~= 1% if ERForcing=0 mean .eg. abrupt experint
    %     ERForcing = repmat(ERForcing, [12 1]);
    %     ERForcing = reshape(ERForcing, 1, length(ERForcing(:))); % yearx12, assume every month have the same radiative forcing
    % end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model member
    for mdlNum = 1:length(level.model2)
        % model path
        mdlName = level.model2{mdlNum};
        mdlPath = fullfile(exmPath, mdlName);
        eval(['cd ', mdlPath]);
        disp(' ')
        disp([mdlName, ' model start!'])
        % ensemble member path
        esmName = getPath_fileName(mdlPath, '.');
        eval(['cd ', nowpath]);

        %% 暂时只看esm实验
        esm = 'r1i1p1f1';

        if sum(strcmp(esmName, esm)) == 0
            disp(['the ', esm, ' ensemble of ', mdlName, ' didnt exist']);
            continue
        end

        specificNum = find(strcmp(esmName, 'r1i1p1f1') == 1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensemble member
        for esmNum = specificNum:specificNum%1:1%length(esmName)% note that r1i1p1 sometime not the first folder
            esmPath = fullfile(mdlPath, esmName{esmNum, 1});
            %% Step1: 加载数据
            varsPath = fullfile(esmPath, level.process3{1}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/rawdata_regrid
            dvarsPath = fullfile(esmPath, level.process3{2}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/anomaly
            dvarsTrendPath = fullfile(esmPath, level.process3{3}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/anomaly_trend
            kernelPath = fullfile(esmPath, level.process3{5}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/kernelsCal
            dradEffectPath = fullfile(esmPath, level.process3{6}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/radEffect/
            dnonLocalCldPath = fullfile(esmPath, level.process3{8}); %/data1/liuyincheng/CMIP6-process/2000-2014/MRI-ESM2-0/non_localCld/

            load([dradEffectPath, 'global_vars.mat'])% lat_f lon_f time plev_k readme
            ntime = length(time.date);
            % load dtsg and DeltaTsg
            load([dvarsPath, 'dtsg_nomask.mat'])% dtsg
            load([dvarsTrendPath, 'DeltaTsg_nomask.mat'])% DeltaTsg, DeltaTsgMeaning

            varUsed = zeros(nlonf, nlatf, ntime, 3); % dR
            varKerUsed = varUsed; % dR/kernels

            % load rad
            load([dradEffectPath, 'dR_cloud.mat'])% cloud radRffect

            % load surface temperature
            load([varsPath, 'ts.mat'])% ts
            load([dvarsPath, 'dts.mat'])%dts
            dts = autoRegrid3(lon_k, lat_k, time.date, dts, lon_f, lat_f, time.date);
            ts = autoRegrid3(lon_k, lat_k, time.date, ts, lon_f, lat_f, time.date);

            % load radiation
            load([dradEffectPath, 'dradEffect_toa_cld.mat'])% 'nonCloudEffect', 'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', 'taEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect', 'mainEffect'
            load([dradEffectPath, 'model_dradEffect.mat'])% dR_net_cld, dR_net_clr, lw_net, readme_realradEfect, sw_net 'dR_allsky', 'l_rad', 's_rad', 'dR_clr', 'readme_realradEfect'
            load([dradEffectPath, 'dR_residual_cld_toa.mat'])% dR_residual_cld_toa
            load([varsPath, 'rlut.mat'])% toa_outgoing_longwave_flux
            load([varsPath, 'rsut.mat'])% toa_outgoing_shortwave_flux
            load([varsPath, 'rsdt.mat'])% toa_incoming_shortwave_flux
            dR_nonCloud = nonCloudEffect;
            dR_netTOA = squeeze(dR_net_cld(:, :, :, 2));
            dR_nonCloud2 = dR_nonCloud + dR_residual_cld_toa;

            netTOA = rsdt - rlut - rsut; % net radiation in toa (lonxlatxmonth): W/m2
            netTOA = autoRegrid3(lon_k, lat_k, time.date, netTOA, lon_f, lat_f, time.date);

            % 将月时间序列转换为年时间序列(lon,lat,1020)→(lon,lat,1020/12)并进行纬向加权平均
            gblM_ts = monthlyToYearlyMean(ts, lat_f);
            gblM_dts = monthlyToYearlyMean(dts, lat_f);

            gblM_dR_netTOA = monthlyToYearlyMean(dR_netTOA, lat_f);
            gblM_R_netTOA = monthlyToYearlyMean(netTOA, lat_f);

            gblM_dR_cloud_toa = monthlyToYearlyMean(dR_cloud_toa, lat_f);
            gblM_dR_nonCloud = monthlyToYearlyMean(dR_nonCloud, lat_f);
            gblM_dR_nonCloud2 = monthlyToYearlyMean(dR_nonCloud2, lat_f);
            gblM_dR_residual = monthlyToYearlyMean(dR_residual_cld_toa, lat_f);
            % 海温基态值(control climate), 1850年为基态, 所以netTOA=0
            % 以1850为基态
            T_ref1850 = 13.9 + 273.15 - 0.373; %K, see https://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/download.html, 13.9 is mean of 1961-1990
            % 预定观测的bts0(初始表面温度)
            T_ref2015Series = [0.763 0.797 0.677 0.587 0.736 0.864];
            T_ref2015Series = T_ref2015Series + 13.9 + 273.15;
            T_ref2015Serie = mean(T_ref2015Series);
            bts0 = T_ref2015Serie - T_ref1850; % 15年观测相对于1850年的变化
            % 计算基于基态的时间序列
            ts15fix = gblM_ts(1) - T_ref2015Serie;
            dts_ctrl = gblM_ts -ts15fix - T_ref1850;
            % comput C and Cd
            [C_ref_Jim, Cd_ref_Jim, C_ref_Geo, Cd_ref_Geo] = calCCd();
            C_ref = C_ref_Jim;
            Cd_ref = Cd_ref_Jim; 
            % 计算初始深海温度(估计方法及计算模块设计见oceanVarEst文件夹)
            [b_td0] = calTd0F(Cd_ref);% Jim=0.1393 Geo=0.34895
            % b_td0 = 0.412; % 2015年相对于1850年(该预设过时且可能存在错误)

            %% Step2: 划分时间段以求拟合线
            % 二层EBM模型
            %% 多个不同时间段的变化
            startY = 0; endY = 84; interY = 1;
            % 计算gamma值以及其他的积分值等
            [integral_ts, integral_RnetTOA, dtd_ctrl, integral_td, gammaa] = EBMLayer2_gamma(C_ref, Cd_ref, gblM_R_netTOA, dts_ctrl, b_td0, startY, endY, interY);

            %% step3: 计算没有云情况下的温度变化
            %% 1) 计算非云气候反馈因子/参数
            kb_RnonCld = polyfit(gblM_dts, gblM_dR_nonCloud, 1);
            lamda_nonCld = abs(kb_RnonCld(1)); % 注意符号(根据反馈因子的定义)

            % test F-lamda*T
            % 云lamda
            kb_RCld = polyfit(gblM_dts, gblM_dR_cloud_toa, 1);
            lamda_Cld = abs(kb_RCld(1)); % 注意符号(根据反馈因子的定义)
            % residual lamda
            kb_Rres = polyfit(gblM_dts, gblM_dR_residual, 1);
            lamda_res = abs(kb_Rres(1)); % 注意符号(根据反馈因子的定义)
            dRTOA_cal = ERForcing' - lamda_nonCld .* (dts_ctrl - dts_ctrl(1)) - lamda_Cld .* (dts_ctrl - dts_ctrl(1)) - lamda_res .* (dts_ctrl - dts_ctrl(1));

            lamda_nonCld = lamda_nonCld + lamda_res;
            save([dnonLocalCldPath, 'lamda_nonCld.mat'], 'lamda_nonCld');

            %% 2) 计算delta_dtnc, dtnc_ctrl, integral_dtnc(积分)  通过迭代的方法计算
            [integral_ERF, dtnc_ctrl, integral_dtnc, delta_dtnc, Delta_dtnc] = EBMLayer2_tnc(ERForcing, lamda_nonCld, gammaa, integral_td, bts0, C_ref, startY, endY, interY);
            % [] = EBMLayer2_tnc_plotTest(startY, endY, interY, dtnc_ctrl, DeltaTsg, Delta_dtnc, gammaa, mdlName);
            DeltaTsg_cld = DeltaTsg - Delta_dtnc;
            %%%%%%%%%%%%%%%%%%%%Part2: 计算局地+全球云的作用(cal dRnonLocalCld3 and trendyr_dRnonLocalCld3(unite: per day))%%%%%%%%%%%%%%%%%%%%%%%%%
            load([dradEffectPath, 'dradEffect_sfc_cld.mat'])%albEffect, husEffect, taEffect, mainEffect, totalEffect, tsEffect, talwEffect, taswEffect
            load([dradEffectPath, 'dR_residual_cld_sfc.mat'])% dR_residual_cld_sfc
            varUsed(:, :, :, 1) = husEffect;
            varUsed(:, :, :, 2) = taEffect;
            varUsed(:, :, :, 3) = albEffect;
            varUsed(:, :, :, 4) = dR_residual_cld_sfc;
            varParitial = zeros(nlonf, nlatf, 4);
            varKerParitial = varParitial;
            % Step1: cal trendyr_dRnonLocalCld_toa/sfc(unite: per day)
            for varNum = 1:4
                tempTemp = squeeze(varUsed(:, :, :, varNum));

                for latNum = 1:nlatf

                    for lonNum = 1:nlonf
                        tempUsed = squeeze(squeeze(tempTemp(lonNum, latNum, :)));
                        k = polyfit(dtsg, tempUsed, 1); % 一元一次拟合
                        varParitial(lonNum, latNum, varNum) = k(1);
                    end

                end

            end

            dRnonLocalCld3_hus = squeeze(varParitial(:, :, 1)) * DeltaTsg_cld;
            dRnonLocalCld3_ta = squeeze(varParitial(:, :, 2)) * DeltaTsg_cld;
            dRnonLocalCld3_alb = squeeze(varParitial(:, :, 3)) * DeltaTsg_cld;
            dRnonLocalCld3_res = squeeze(varParitial(:, :, 4)) * DeltaTsg_cld;
            dRnonLocalCld3 = dRnonLocalCld3_hus + dRnonLocalCld3_ta + dRnonLocalCld3_alb + dRnonLocalCld3_res;

            trendyr_dRnonLocalCld3 = dRnonLocalCld3 ./ (time.date(end) - time.date(1));
            save([dnonLocalCldPath, 'dRnonLocalCld3_sfc.mat'], 'dRnonLocalCld3_hus', 'dRnonLocalCld3_ta', 'dRnonLocalCld3_alb', 'dRnonLocalCld3_res', 'dRnonLocalCld3')
            save([dnonLocalCldPath, 'trendyr_dRnonLocalCld3_sfc.mat'], 'trendyr_dRnonLocalCld3')

            % Step2: cal trendyr_dTsnonLocalCld3_toa/sfc(unite: per day)
            % load ts kernel
            load([kernelPath, '/kernels_sfc_cld'], 'ts_lwkernel');
            load([kernelPath, '/global_vars.mat']); % lat_k,lon_k,time
            % regrid ts_lwkernel to 144x72(unite grids)
            kernelTime = 1:12;
            [Xlon, Ylat, Ttime] = meshgrid(lat_f, lon_f, kernelTime);
            [Xlonf, Ylatf, Ttimef] = meshgrid(lat_k, lon_k, kernelTime);
            ts_lwkernel = interp3(Xlonf, Ylatf, Ttimef, ts_lwkernel, Xlon, Ylat, Ttime);
            % extend to the whole time series
            startmonth = 1;

            if startmonth ~= 1
                tempK = zeros(nlonf, nlatf, 12);
                tempK(1:12 - startmonth + 1) = ts_lwkernel(:, :, startmonth:12);
                tempK(12 - startmonth + 2:12) = ts_lwkernel(:, :, 1:startmonth - 1);
                ts_lwkernel = tempK;
            end

            nyear = ntime / 12;
            ts_lwkernelSeries = repmat(ts_lwkernel, [1 1 nyear]);

            for varNum = 1:4
                varKerUsed(:, :, :, varNum) = squeeze(varUsed(:, :, :, varNum)) ./ -ts_lwkernelSeries;
            end

            for varNum = 1:4
                tempTemp = squeeze(varKerUsed(:, :, :, varNum));

                for latNum = 1:nlatf

                    for lonNum = 1:nlonf
                        tempUsed = squeeze(squeeze(tempTemp(lonNum, latNum, :)));
                        k = polyfit(dtsg, tempUsed, 1); % 一元一次拟合
                        varKerParitial(lonNum, latNum, varNum) = k(1);
                    end

                end

            end

            dTsnonLocalCld3_hus = squeeze(varKerParitial(:, :, 1)) * DeltaTsg_cld;
            dTsnonLocalCld3_ta = squeeze(varKerParitial(:, :, 2)) * DeltaTsg_cld;
            dTsnonLocalCld3_alb = squeeze(varKerParitial(:, :, 3)) * DeltaTsg_cld;
            dTsnonLocalCld3_res = squeeze(varKerParitial(:, :, 4)) * DeltaTsg_cld;

            dTsnonLocalCld3 = dTsnonLocalCld3_hus + dTsnonLocalCld3_ta + dTsnonLocalCld3_alb + dTsnonLocalCld3_res;
            trendyr_dTsnonLocalCld3 = dTsnonLocalCld3 ./ (time.date(end) - time.date(1));

            save([dnonLocalCldPath, 'dTsnonLocalCld3_sfc.mat'], 'dTsnonLocalCld3_hus', 'dTsnonLocalCld3_ta', 'dTsnonLocalCld3_alb', 'dTsnonLocalCld3_res', 'dTsnonLocalCld3')
            save([dnonLocalCldPath, 'trendyr_dTsnonLocalCld3_sfc.mat'], 'trendyr_dTsnonLocalCld3')

            disp([esmName{esmNum, 1}, ' ensemble is done!'])
        end

        disp([mdlName, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{exmNum}, ' era is done!'])
    disp(' ')

end

t = toc; disp(t)

function [integral_ts, integral_RnetTOA, dtd_ctrl, integral_td, gammaa] = EBMLayer2_gamma(C_ref, Cd_ref, gblM_R_netTOA, dts_ctrl, b_td0, startY, endY, interY)
    % 计算gamma值
    % switch nargin
    %     case 8
    %         figName = [];
    % end

    timeExamp = startY + interY * 2 - 1:interY:endY;
    sum_timeExamp = length(timeExamp); % 总共有几组时间段

    delta_ts = zeros(sum_timeExamp + 1, 1); % Delta_ts(1)=0
    integral_delta_ts = zeros(sum_timeExamp + 1, 1);
    integral_ts = zeros(sum_timeExamp + 1, 1);

    delta_td = zeros(sum_timeExamp + 1, 1);
    integral_delta_td = zeros(sum_timeExamp + 1, 1);
    dtd_ctrl = zeros(sum_timeExamp + 1, 1);

    integral_RnetTOA = zeros(sum_timeExamp + 1, 1); % integral_RnetTOA=sum(RnetTOA)

    % 设定初始值
    delta_td(1) = 0;
    integral_delta_td(1) = 0;
    dtd_ctrl(1) = b_td0;

    delta_ts(1) = 0;
    integral_delta_ts = dts_ctrl - dts_ctrl(1);
    integral_ts(1) = dts_ctrl(1);

    integral_RnetTOA(1) = gblM_R_netTOA(1);

    % 计算每个格点的值
    for gridNum = 2:sum_timeExamp + 1% 循环每一个格点

        Seri_RnetTOA = gblM_R_netTOA(1:timeExamp(gridNum - 1) - startY + 1);
        Seri_dts = dts_ctrl(1:timeExamp(gridNum - 1) - startY + 1);

        integral_ts(gridNum) = sum(Seri_dts);
        integral_RnetTOA(gridNum) = sum(Seri_RnetTOA);
        delta_td(gridNum) = (integral_RnetTOA(gridNum) - C_ref * integral_delta_ts(gridNum)) / Cd_ref - integral_delta_td(gridNum - 1);
        integral_delta_td(gridNum) = integral_delta_td(gridNum - 1) + delta_td(gridNum);
        dtd_ctrl = dtd_ctrl(1) + integral_delta_td;
    end

    % 计算Td积分
    integral_td = zeros(sum_timeExamp + 1, 1);
    integral_td(1) = b_td0;

    for gridNum = 2:sum_timeExamp + 1%循环每一组
        Seri_dtd = dtd_ctrl(1:timeExamp(gridNum - 1) - startY + 1);
        integral_td(gridNum) = sum(Seri_dtd);
    end

    % 计算gamma (function 2, see document for details)
    xdata = integral_ts - integral_td;
    ydata = integral_delta_td .* Cd_ref;
    % gamma1=ydata./xdata;

    % 线性拟合
    xdata = xdata(40:end);
    ydata = ydata(40:end);
    kb = polyfit(xdata, ydata, 1);
    k_fit = kb(1);
    b_fit = kb(2);

    % % 线性拟合(过零点)
    % [k_fit, b_fit] = plotLinFit_point0(xdata, ydata, xLabel, yLabel, mdlName);
    gammaa = k_fit;
    resiual = b_fit;

end

function [integral_ERF, dtnc_ctrl, integral_dtnc, delta_dtnc, Delta_dtnc] = EBMLayer2_tnc(ERForcing, lamda_nonCld, gammaa, integral_td, bts0, C_ref, startY, endY, interY)
    % description.
    % 计算在没有云的变化情况下, 温度的变化

    %% 2) 计算delta_dtnc, dtnc_ctrl, integral_dtnc(积分)  通过迭代的方法计算
    timeExamp = startY + interY * 2 - 1:interY:endY;
    sum_timeExamp = length(timeExamp); % 总共有几组时间段
    % 计算ERF积分
    integral_ERF = zeros(sum_timeExamp + 1, 1); %note time=0, delta_tnc=0;
    integral_ERF(1) = ERForcing(1);

    for gridNum = 2:sum_timeExamp + 1
        Seri_ERF = ERForcing(1:timeExamp(gridNum - 1) - startY + 1);
        integral_ERF(gridNum) = sum(Seri_ERF);
    end

    % 准备迭代
    dtnc_ctrl = zeros(sum_timeExamp + 1, 1);
    integral_dtnc = zeros(sum_timeExamp + 1, 1);
    delta_dtnc = zeros(sum_timeExamp + 1, 1); %note time=0, delta_tnc=0;
    integral_delta_dtnc = zeros(sum_timeExamp + 1, 1); %note time=0, delta_tnc=0;

    dtnc_ctrl(1) = bts0;
    integral_dtnc(1) = bts0;
    integral_delta_dtnc(1) = 0;

    for gridNum = 2:sum_timeExamp + 1
        delta_dtnc(gridNum) = (integral_ERF(gridNum) + gammaa * integral_td(gridNum) - (lamda_nonCld + gammaa) * (integral_dtnc(gridNum - 1) + dtnc_ctrl(gridNum - 1))) / (C_ref + lamda_nonCld + gammaa) - integral_delta_dtnc(gridNum - 1);
        dtnc_ctrl(gridNum) = dtnc_ctrl(gridNum - 1) + delta_dtnc(gridNum);
        integral_dtnc(gridNum) = integral_dtnc(gridNum - 1) + dtnc_ctrl(gridNum);
        integral_delta_dtnc(gridNum) = integral_delta_dtnc(gridNum - 1) + delta_dtnc(gridNum);
    end

    timeSeries = startY:interY:endY;

    k_dtnc = polyfit(timeSeries(40:end)', dtnc_ctrl(40:end), 1);
    Delta_dtnc = k_dtnc(1) * (timeSeries(end) - timeSeries(1));

end

function [] = EBMLayer2_tnc_plotTest(startY, endY, interY, dtnc_ctrl, DeltaTsg, Delta_dtnc, gammaa, mdlName, figName)
    % 计算没有云情况下, 温度变化的斜率
    switch nargin
        case 8
            figName = [];
    end

    xdata = startY:interY:endY;
    ydata = dtnc_ctrl';
    yLabel = 'tnc(非云致温度变化)';
    xLabel = 'time';

    % 线性拟合
    [k_fit, ~] = plotLinFit(xdata, ydata, xLabel, yLabel, mdlName, figName);

    k_tnc = k_fit;

    DeltaTsg_cld = DeltaTsg - Delta_dtnc;

    title({['ssp370 Model: ', mdlName]; ['k\_tnc: ', num2str(k_tnc), ' gamma:', num2str(gammaa)]; ['无云情况下温度变化:', num2str(Delta_dtnc), '云致温度变化:', num2str(DeltaTsg_cld)]})
    figName = ['/home/liuyc/Research/P02.Ts_change_research/figure/02.cmip6Result/1.6/nonLocalCld3/tempChangeTest/ssp370_r1i1p1f1_', mdlName, '_Jim.png'];
    saveas(gcf, figName)

end

function [C_ref_Jim, Cd_ref_Jim, C_ref_Geo, Cd_ref_Geo] = calCCd()
    %% computation of C and Cd
    %% model parameters
    D_ref_Jim = 55; Dd_ref_Jim = 2768; % D_ref = 77; Dd_ref = 1105; in Geoffroy and D_ref = 55; Dd_ref = 2768; in Jiménez-de-la-Cuesta,
    D_ref_Geo = 77; Dd_ref_Geo = 1105;
    cp_water = 4180; % heat capacity J/kg/K
    rho_water = 1030; % desnsity of saltwater kg/m3
    f0 = 0.7; % The proportion of ocean cover the earth surface
    % cal reference C and C0
    C_ref_Jim = (D_ref_Jim * rho_water * cp_water * f0) / (86400 * 365.25); % C_ref(upper-ocean) mean the result of formula (22) in reference Geoffroy
    Cd_ref_Jim = (Dd_ref_Jim * rho_water * cp_water * f0) / (86400 * 365.25); % Cd_ref(deep-layer ocean) mean the result of formula (22) in reference Geoffroy

    C_ref_Geo = (D_ref_Geo * rho_water * cp_water * f0) / (86400 * 365.25); % C_ref(upper-ocean) mean the result of formula (22) in reference Geoffroy
    Cd_ref_Geo = (Dd_ref_Geo * rho_water * cp_water * f0) / (86400 * 365.25); % Cd_ref(deep-layer ocean) mean the result of formula (22) in reference Geoffroy

end

function [Td0] = calTd0F(Cd_ref)
    % Td0 mean deep ocean temperature in 2015
    heatCPath = '/data1/liuyincheng/Observe-rawdata/NCEI/Ocean Heat/';

    heatC_2000m = ncread([heatCPath, 'heat_content_anomaly_0-2000_yearly.nc'], 'yearl_h22_WO');
    heatC_2000m = heatC_2000m(1:end - 1);

    heatC_2000m_2005 = heatC_2000m(end - 13:end);

    % time line
    timeHC_2005 = 2005:2019 - 1;
    timeHC_1850 = 1850:2019 - 1;

    % see showProcess_calTd0.m for more details about method select
    ydata = heatC_2000m_2005;
    xdata = timeHC_2005;
    k_2000_mul = polyfit(xdata', log(ydata + 20), 1);
    xdata_predict = timeHC_1850;
    heatC_2000m_mulFit_predict = exp(polyval(k_2000_mul, xdata_predict)) - 20;

    %% computation of C and Cd
    heatC_2000_2015year = heatC_2000m(timeHC_2005 == 2015);
    heatC_2000_1850year = heatC_2000m_mulFit_predict(1);
    Td0 = (heatC_2000_2015year - heatC_2000_1850year) * 1e22 / (365.25 * 86400 * Cd_ref * 3.619e8 * 1e3^2);
    disp(['Td(2015)=', num2str(Td0)])

end