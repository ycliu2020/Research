%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-08-15 10:27:16
% LastEditTime : 2020-09-04 16:30:52
% LastEditors  : LYC
% Description  :
% FilePath     : /code/p2_processCMIP6Data/s3.nonLocalCld/s2_nonLocalCld_method3_v1.m
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

for exmNum = 4:4%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    nowpath = pwd;
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % exmPath
    exmPath = ['/data1/liuyincheng/cmip6-process/', level.time1{exmNum}]; %/data1/liuyincheng/cmip6-process/amip_2000-2014/
    %%% load ERF data
    load('/data1/liuyincheng/cmip6-process/z_globalVar/ERF_rec.mat')% 'ERF_rec', 'timeERF_rec'
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
    for mdlNum = 9:9 %1:length(level.model2)
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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% load and read 
            esmPath = fullfile(mdlPath, esmName{esmNum, 1});
            % data path
            varsPath = fullfile(esmPath, level.process3{1}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/rawdata_regrid
            dvarsPath = fullfile(esmPath, level.process3{2}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
            dvarsTrendPath = fullfile(esmPath, level.process3{3}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly_trend
            kernelPath = fullfile(esmPath, level.process3{5}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/kernelsCal
            dradEffectPath = fullfile(esmPath, level.process3{6}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
            dnonLocalCldPath = fullfile(esmPath, level.process3{8}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/non_localCld/

            load([dradEffectPath, 'global_vars.mat'])% lat_f lon_f time plev_k readme
            ntime = length(time.date);
            % load dtsg and DeltaTsg
            load([dvarsPath, 'dtsg_nomask.mat'])% dtsg
            load([dvarsTrendPath, 'DeltaTsg_nomask.mat'])% DeltaTsg, DeltaTsgMeaning

            varUsed = zeros(nlonf, nlatf, ntime, 3); % dR
            varKerUsed = varUsed; % dR/kernels

            % cal rad
            load([dradEffectPath, 'dR_cloud_toa.mat'])% cloud radRffect

            % cal surface temperature
            load([varsPath, 'ts.mat'])% ts
            load([dvarsPath, 'dts.mat'])%dts
            dts = autoRegrid3(lon_k, lat_k, time.date, dts, lon_f, lat_f, time.date);
            ts = autoRegrid3(lon_k, lat_k, time.date, ts, lon_f, lat_f, time.date);

            % cal  radiation
            load([dradEffectPath, 'dradEfect_toa_cld.mat'])% 'totalEffect', 'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', 'taEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect', 'mainEffect'
            load([dradEffectPath, 'real_dradEfect.mat'])% 'dR_allsky', 'l_rad', 's_rad', 'dR_clr', 'readme_realradEfect'
            load([dradEffectPath, 'dR_residual_cld_sfc.mat'])% dR_residual_cld_sfc
            load([dradEffectPath, 'dR_residual_cld_toa.mat'])% dR_residual_cld_toa
            load([varsPath, 'rlut.mat'])% toa_outgoing_longwave_flux
            load([varsPath, 'rsut.mat'])% toa_outgoing_shortwave_flux
            load([varsPath, 'rsdt.mat'])% toa_incoming_shortwave_flux
            dR_nonCloud = totalEffect;
            dR_netTOA = squeeze(dR_allsky(:, :, :, 2));
            
            % regrid to 144x72(unite grids)
            netTOA = rsdt - rlut - rsut; % net radiation in toa (lonxlatxmonth): W/m2
            netTOA = autoRegrid3(lon_k, lat_k, time.date, netTOA, lon_f, lat_f, time.date);
            %%%%%%%%%%%%%%%%%%%%Part1: 2层EBM模式计算非云和云致温度变化(分离云对全球的整体作用)%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Step1: 设定初始值
            % 将月时间序列转换为年时间序列(lon,lat,1020)→(lon,lat,1020/12)并进行纬向加权平均
            gblM_ts = monthlyToYearlyMean(ts, lat_f);
            gblM_dts = monthlyToYearlyMean(dts, lat_f);

            gblM_dR_netTOA = monthlyToYearlyMean(dR_netTOA, lat_f);
            gblM_R_netTOA = monthlyToYearlyMean(netTOA, lat_f);

            gblM_dR_cloud_toa = monthlyToYearlyMean(dR_cloud_toa, lat_f);
            gblM_dR_nonCloud = monthlyToYearlyMean(dR_nonCloud, lat_f);

            %% 基态值(control climate)
            % 以1850年的数值为基态
            T_ref1850 = 13.9 + 273.15 - 0.373; %K, see https://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/download.html, 13.9 is mean of 1961-1990
            % 预定观测的bts0(初始表面温度)
            T_ref2015Series = [0.763 0.797 0.677 0.587 0.736 0.864];
            T_ref2015Series = T_ref2015Series + 13.9 + 273.15;
            bts0 = mean((T_ref2015Series - T_ref1850)); % 用观测的
            % 预定初始深海温度(见oceanVarEst文件夹)
            b_td0 = 0.412; % K 相对于1850年
            % 1850年为基态, 所以netTOA=0
            
            % 计算基于基态的时间序列
            dts_ctrl = gblM_dts - T_ref1850;

            %% Step2: 二层EBM模型拟合计算gammma值
            % model parameters
            D_ref = 55; Dd_ref = 2768; %;% D_ref = 77; Dd_ref = 1105; in Geoffroy and D_ref = 55; Dd_ref = 2768; in Jiménez-de-la-Cuesta,
            cp_water = 4180; % heat capacity J/kg/K
            rho_water = 1030; % desnsity of saltwater kg/m3
            f0 = 0.7; % The proportion of ocean cover the earth surface
            % cal reference C and C0
            C_ref = (D_ref * rho_water * cp_water * f0) / (86400 * 365.25); % C_ref(upper-ocean) mean the result of formula (22) in reference Geoffroy
            Cd_ref = (Dd_ref * rho_water * cp_water * f0) / (86400 * 365.25); % Cd_ref(deep-layer ocean) mean the result of formula (22) in reference Geoffroy

            startY = 0; endY = 84; interY = 1;
            [Seri_k_dts, Seri_k_dtd, Seri_b_dRnetTOA, gammaa] = EBMLayer2_gamma2(C_ref, Cd_ref, gblM_R_netTOA, dts_ctrl, bts0, b_td0, startY, endY, interY);
            save([dnonLocalCldPath, 'gamma-EBM.mat'], 'gammaa');

            %% step3: 计算没有云情况下的温度变化
            %% 计算非云气候反馈因子/参数
            kb_RnonCld = polyfit(gblM_dts, gblM_dR_nonCloud, 1);
            lamda_nonCld=-kb_RnonCld(1); % 注意符号(根据反馈因子的定义)
            save([dnonLocalCldPath, 'lamda_nonCld.mat'], 'lamda_nonCld');
            %% 计算没有云情况下温度变化的斜率及其变化
            % [k_tnc, DeltaTsg_noncld, DeltaTsg_cld] = EBMLayer2_ktnc(gammaa, lamda_nonCld, ERForcing, C_ref, Cd_ref, dR_netTOA, dts, bts0, b_td0, Seri_k_dts, Seri_k_dtd, startY, endY, interY, DeltaTsg);
            [k_tnc, DeltaTsg_noncld, DeltaTsg_cld] = EBMLayer2_ktnc_Test(gammaa, lamda_nonCld, ERForcing, C_ref, Cd_ref, dR_netTOA, dts, bts0, b_td0, Seri_k_dts, Seri_k_dtd, startY, endY, interY, DeltaTsg, mdlName);
            % save
            save([dradEffectPath, 'DeltaTsg_cld.mat'], 'DeltaTsg_cld', 'DeltaTsg_noncld', 'DeltaTsg')

            % %%%%%%%%%%%%%%%%%%%%Part2: 计算局地+全球云的作用(cal dRnonLocalCld3 and trendyr_dRnonLocalCld3(unite: per day))%%%%%%%%%%%%%%%%%%%%%%%%%
            % load([dradEffectPath, 'dradEfect_sfc_cld.mat'])%albEffect, husEffect, taEffect, mainEffect, totalEffect, tsEffect, talwEffect, taswEffect
            % varUsed(:, :, :, 1) = husEffect;
            % varUsed(:, :, :, 2) = taEffect;
            % varUsed(:, :, :, 3) = albEffect;
            % varParitial = zeros(nlonf, nlatf, 3);
            % varKerParitial = varParitial;
            % % Step1: cal trendyr_dRnonLocalCld_toa/sfc(unite: per day)
            % for varNum = 1:3
            %     tempTemp = squeeze(varUsed(:, :, :, varNum));

            %     for latNum = 1:nlatf

            %         for lonNum = 1:nlonf
            %             tempUsed = squeeze(squeeze(tempTemp(lonNum, latNum, :)));
            %             k = polyfit(dtsg, tempUsed, 1); % 一元一次拟合
            %             varParitial(lonNum, latNum, varNum) = k(1);
            %         end

            %     end

            % end

            % dRnonLocalCld3_hus = squeeze(varParitial(:, :, 1)) * DeltaTsg_cld;
            % dRnonLocalCld3_ta = squeeze(varParitial(:, :, 2)) * DeltaTsg_cld;
            % dRnonLocalCld3_alb = squeeze(varParitial(:, :, 3)) * DeltaTsg_cld;
            % dRnonLocalCld3 = dRnonLocalCld3_hus + dRnonLocalCld3_ta + dRnonLocalCld3_alb;

            % trendyr_dRnonLocalCld3 = dRnonLocalCld3 ./ (time.date(end) - time.date(1));
            % save([dnonLocalCldPath, 'dRnonLocalCld3_sfc.mat'], 'dRnonLocalCld3_hus', 'dRnonLocalCld3_ta', 'dRnonLocalCld3_alb', 'dRnonLocalCld3')
            % save([dnonLocalCldPath, 'trendyr_dRnonLocalCld3_sfc.mat'], 'trendyr_dRnonLocalCld3')

            % % Step2: cal trendyr_dTsnonLocalCld3_toa/sfc(unite: per day)
            % % load ts kernel
            % load([kernelPath, '/kernels_sfc_cld'], 'ts_lwkernel');
            % load([kernelPath, '/global_vars.mat']); % lat_k,lon_k,time
            % % regrid ts_lwkernel to 144x72(unite grids)
            % kernelTime = 1:12;
            % [Xlon, Ylat, Ttime] = meshgrid(lat_f, lon_f, kernelTime);
            % [Xlonf, Ylatf, Ttimef] = meshgrid(lat_k, lon_k, kernelTime);
            % ts_lwkernel = interp3(Xlonf, Ylatf, Ttimef, ts_lwkernel, Xlon, Ylat, Ttime);
            % % extend to the whole time series
            % startmonth = 1;

            % if startmonth ~= 1
            %     tempK = zeros(nlonf, nlatf, 12);
            %     tempK(1:12 - startmonth + 1) = ts_lwkernel(:, :, startmonth:12);
            %     tempK(12 - startmonth + 2:12) = ts_lwkernel(:, :, 1:startmonth - 1);
            %     ts_lwkernel = tempK;
            % end

            % nyear = ntime / 12;
            % ts_lwkernelSeries = repmat(ts_lwkernel, [1 1 nyear]);

            % for varNum = 1:3
            %     varKerUsed(:, :, :, varNum) = squeeze(varUsed(:, :, :, varNum)) ./ -ts_lwkernelSeries;
            % end

            % for varNum = 1:3
            %     tempTemp = squeeze(varKerUsed(:, :, :, varNum));

            %     for latNum = 1:nlatf

            %         for lonNum = 1:nlonf
            %             tempUsed = squeeze(squeeze(tempTemp(lonNum, latNum, :)));
            %             k = polyfit(dtsg, tempUsed, 1); % 一元一次拟合
            %             varKerParitial(lonNum, latNum, varNum) = k(1);
            %         end

            %     end

            % end

            % dTsnonLocalCld3_hus = squeeze(varKerParitial(:, :, 1)) * DeltaTsg_cld;
            % dTsnonLocalCld3_ta = squeeze(varKerParitial(:, :, 2)) * DeltaTsg_cld;
            % dTsnonLocalCld3_alb = squeeze(varKerParitial(:, :, 3)) * DeltaTsg_cld;

            % dTsnonLocalCld3 = dTsnonLocalCld3_hus + dTsnonLocalCld3_ta + dTsnonLocalCld3_alb;
            % trendyr_dTsnonLocalCld3 = dTsnonLocalCld3 ./ (time.date(end) - time.date(1));

            % save([dnonLocalCldPath, 'dTsnonLocalCld3_sfc.mat'], 'dTsnonLocalCld3_hus', 'dTsnonLocalCld3_ta', 'dTsnonLocalCld3_alb', 'dTsnonLocalCld3')
            % save([dnonLocalCldPath, 'trendyr_dTsnonLocalCld3_sfc.mat'], 'trendyr_dTsnonLocalCld3')

            disp([esmName{esmNum, 1}, ' ensemble is done!'])
        end

        disp([mdlName, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{exmNum}, ' era is done!'])
    disp(' ')

end

t = toc; disp(t)

function [Seri_k_dts, Seri_k_dtd, Seri_b_dRnetTOA, gammaa] = EBMLayer2_gamma(C_ref, Cd_ref, dR_netTOA, dts, bts0, b_td0, startT, endT, interT)
    % 计算gamma值
    % switch nargin
    %     case 8
    %         figName = [];
    % end

    timeExamp = startT + interT * 4 - 1:interT:endT;
    diffTimeSum = length(timeExamp); % 总共有几组时间段

    Delta_ts = zeros(diffTimeSum, 1);
    Seri_k_dtd = zeros(diffTimeSum, 1);

    Seri_k_dts = zeros(diffTimeSum, 1);
    Seri_b_dts = zeros(diffTimeSum, 1);

    Seri_k_dRnetTOA = zeros(diffTimeSum, 1);
    Seri_b_dRnetTOA = zeros(diffTimeSum, 1);

    p_Plus = zeros(diffTimeSum, 1);
    p_Minus = zeros(diffTimeSum, 1);

    % 计算k_dtd (function 1)
    for diffTime = 1:diffTimeSum%循环每一组
        timeSeries = startT:timeExamp(diffTime);
        lenTime = length(timeSeries);

        Seri_dRnetTOA_ctrl = dR_netTOA(1:timeExamp(diffTime) - startT + 1);
        Seri_dts_ctrl = dts(1:timeExamp(diffTime) - startT + 1);

        %% cal Delta ts
        k_dts = polyfit(timeSeries', Seri_dts_ctrl, 1);
        Delta_ts(diffTime) = k_dts(1) * lenTime;
        Seri_k_dts(diffTime) = k_dts(1);
        Seri_b_dts(diffTime) = bts0;

        k_dRnetTOA = polyfit(timeSeries', Seri_dRnetTOA_ctrl, 1);
        Seri_k_dRnetTOA(diffTime) = k_dRnetTOA(1);
        Seri_b_dRnetTOA(diffTime) = k_dRnetTOA(2);
        p_Plus(diffTime) = timeExamp(diffTime) + startT; %(t1+t2)
        p_Minus(diffTime) = timeExamp(diffTime) - startT; %(t1+t2)
        Seri_k_dtd(diffTime) = (0.5 * Seri_k_dRnetTOA(diffTime) * p_Plus(diffTime) + Seri_b_dRnetTOA(diffTime) - C_ref * Seri_k_dts(diffTime)) / Cd_ref;
    end

    % 计算gamma (function 2, see document for details)
    xdata = 0.5 * (Seri_k_dts - Seri_k_dtd) .* p_Plus + Seri_b_dts - b_td0;
    ydata = Seri_k_dtd .* Cd_ref;
    % gamma1=ydata./xdata;

    % 线性拟合
    kb = polyfit(xdata, ydata, 1);
    k_fit = kb(1);
    b_fit = kb(2);

    % % 线性拟合(过零点)
    % [k_fit, b_fit] = plotLinFit_point0(xdata, ydata, xLabel, yLabel, modelName);
    gammaa = k_fit;
    resiual = b_fit;

end
function [Seri_k_dts, Seri_k_dtd, Seri_b_dRnetTOA, gammaa] = EBMLayer2_gamma2(C_ref, Cd_ref, dR_netTOA, dts, bts0, b_td0, startT, endT, interT)
    % 计算gamma值
    % switch nargin
    %     case 8
    %         figName = [];
    % end

    timeExamp = startT + interT * 4 - 1:interT:endT;
    diffTimeSum = length(timeExamp); % 总共有几组时间段

    Delta_ts = zeros(diffTimeSum, 1);
    Seri_k_dtd = zeros(diffTimeSum, 1);

    Seri_k_dts = zeros(diffTimeSum, 1);
    Seri_b_dts = zeros(diffTimeSum, 1);

    Seri_k_dRnetTOA = zeros(diffTimeSum, 1);
    Seri_b_dRnetTOA = zeros(diffTimeSum, 1);

    p_Plus = zeros(diffTimeSum, 1);
    p_Minus = zeros(diffTimeSum, 1);

    % 计算k_dtd (function 1)
    for diffTime = 1:diffTimeSum%循环每一组
        timeSeries = startT:timeExamp(diffTime);
        dtime = length(timeSeries);

        Seri_dRnetTOA_ctrl = dR_netTOA(1:timeExamp(diffTime) - startT + 1);
        Seri_dts_ctrl = dts(1:timeExamp(diffTime) - startT + 1);

        %% cal Delta ts
        k_dts = polyfit(timeSeries', Seri_dts_ctrl, 1);
        Delta_ts(diffTime) = k_dts(1) * dtime;
        Seri_k_dts(diffTime) = k_dts(1);
        Seri_b_dts(diffTime) = bts0;

        p_Plus(diffTime) = timeExamp(diffTime) + startT; %(t1+t2)
        p_Minus(diffTime) = timeExamp(diffTime) - startT; %(t1+t2)
        Seri_k_dtd(diffTime) = (sum(Seri_dRnetTOA_ctrl) / p_Minus(diffTime) - C_ref * Seri_k_dts(diffTime)) / Cd_ref;
    end

    % 计算gamma (function 2, see document for details)
    xdata = 0.5 * (Seri_k_dts - Seri_k_dtd) .* p_Plus + Seri_b_dts - b_td0;
    ydata = Seri_k_dtd .* Cd_ref;
    % gamma1=ydata./xdata;

    % 线性拟合
    kb = polyfit(xdata, ydata, 1);
    k_fit = kb(1);
    b_fit = kb(2);

    % % 线性拟合(过零点)
    % [k_fit, b_fit] = plotLinFit_point0(xdata, ydata, xLabel, yLabel, modelName);
    gammaa = k_fit;
    resiual = b_fit;

end

function [k_tnc, DeltaTsg_noncld, DeltaTsg_cld] = EBMLayer2_ktnc(gammaa, lamda_nonCld, ERForcing, C_ref, Cd_ref, dR_netTOA, dts, bts0, b_td0, Seri_k_dts, Seri_k_dtd, startT, endT, interT, DeltaTsg, figName)
    % 计算没有云情况下, 温度变化的斜率
    switch nargin
        case 16
            figName = [];
    end

    timeExamp = startT + interT * 4 - 1:interT:endT;
    diffTimeSum = length(timeExamp); % 总共有几组时间段

    sum_ERF = zeros(diffTimeSum, 1);
    
    p_Plus = zeros(diffTimeSum, 1);
    p_Minus = zeros(diffTimeSum, 1);
    
    % 计算(function 7)
    for diffTime = 1:diffTimeSum%循环每一组
        timeSeries = startT:timeExamp(diffTime);
        lenTime = length(timeSeries);

        p_Plus(diffTime) = timeExamp(diffTime) + startT; %(t1+t2)
        p_Minus(diffTime) = timeExamp(diffTime) - startT; %(t1+t2)
        
        sum_ERF(diffTime)=sum(ERForcing(1:lenTime));
    end
    % 计算Seri_k_tnc
    xdata=C_ref+0.5*p_Plus.*(lamda_nonCld+gammaa);
    ydata=(sum_ERF./p_Minus)+gammaa*(0.5*p_Plus.*Seri_k_dtd+b_td0)-(lamda_nonCld+gammaa)*bts0;
    kb = polyfit(xdata, ydata, 1);
    k_fit = kb(1);
    b_fit = kb(2);

    k_tnc=k_fit;
    resiual = b_fit;
    DeltaTsg_noncld=k_tnc*(endT - startT);
    DeltaTsg_cld=DeltaTsg-DeltaTsg_noncld;

end

function [k_tnc, DeltaTsg_noncld, DeltaTsg_cld] = EBMLayer2_ktnc_Test(gammaa, lamda_nonCld, ERForcing, C_ref, Cd_ref, dR_netTOA, dts, bts0, b_td0, Seri_k_dts, Seri_k_dtd, startT, endT, interT, DeltaTsg, modelName, figName)
    % 计算没有云情况下, 温度变化的斜率
    switch nargin
        case 16
            figName = [];
    end

    timeExamp = startT + interT * 4 - 1:interT:endT;
    diffTimeSum = length(timeExamp); % 总共有几组时间段
    
    p_Plus = zeros(diffTimeSum, 1);
    p_Minus = zeros(diffTimeSum, 1);
    
    % 计算ERF(function 7)
    for diffTime = 1:diffTimeSum%循环每一组
        timeSeries = startT:timeExamp(diffTime);
        lenTime = length(timeSeries);

        Seri_erf = ERForcing(1:timeExamp(diffTime) - startT + 1);

        %% cal Delta ERF
        k_f = polyfit(timeSeries', Seri_erf, 1);
        Delta_ts(diffTime) = k_dts(1) * dtime;
        Seri_k_f(diffTime) = k_f(1);
        Seri_b_dts(diffTime) = b;


        Seri_k_dtd(diffTime) = (sum(Seri_dRnetTOA_ctrl) / p_Minus(diffTime) - C_ref * Seri_k_dts(diffTime)) / Cd_ref;

        p_Plus(diffTime) = timeExamp(diffTime) + startT; %(t1+t2)
        p_Minus(diffTime) = timeExamp(diffTime) - startT; %(t1+t2)
        
        sum_ERF(diffTime)=sum(ERForcing(1:lenTime));
    end
    % 计算Seri_k_tnc
    xdata=C_ref+0.5*p_Plus.*(lamda_nonCld+gammaa);
    ydata=(sum_ERF./p_Minus)+gammaa*(0.5*p_Plus.*Seri_k_dtd+b_td0)-(lamda_nonCld+gammaa)*bts0;

    yLabel = 'y';
    xLabel = 'x';

    % 线性拟合
    [k_fit, b_fit] = plotLinFit(xdata, ydata, xLabel, yLabel, modelName, figName);
    
    k_tnc=k_fit;
    DeltaTsg_noncld=k_tnc*(endT - startT);
    DeltaTsg_cld=DeltaTsg-DeltaTsg_noncld;

    resiual = b_fit;
    title({['ssp370 Model: ', modelName]; ['k_tnc: ', num2str(k_tnc), '  Res(认为截距为误差): ', num2str(resiual)];['无云情况下温度变化:',num2str(DeltaTsg_noncld),'云致温度变化:',num2str(DeltaTsg_cld)]})
    figName = ['/home/liuyc/Research/P02.Ts_change_research/figure/02.cmip6Result/1.6/nonLocalCld3/tempChangeTest/ssp370_r1i1p1f1_', modelName, '_Jim.png'];
    saveas(gcf, figName)

end