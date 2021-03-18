%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-09-22 21:51:51
% LastEditTime : 2021-03-12 16:59:31
% LastEditors  : Please set LastEditors
% Description  :
% FilePath     : /code/p3_paperFigIntegrate/Fig11_nonLocalCld3_methodShow.m
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
outPutPath = '/home/liuyc/Research/P02.Ts_change_research/figure/proj3_PaperFig/v0.0/Fig11_CMIP6_showNonLocalCldCalcProcess/';

for exmNum = 4:4%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
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

            load([dnonLocalCldPath, 'lamda_nonCld.mat'])% lat_f lon_f time plev_k readme

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
            C_name = 'Jim';
            if strcmp(C_name, 'Jim') == 1
                C_ref = C_ref_Jim;
                Cd_ref = Cd_ref_Jim;
            elseif strcmp(C_name, 'Geo') == 1
                C_ref = C_ref_Geo;
                Cd_ref = Cd_ref_Geo;
            end

            % 计算初始深海温度(估计方法及计算模块设计见oceanVarEst文件夹)
            [b_td0] = calTd0F(Cd_ref); % Jim=0.1393K Geo=0.34895K
            % b_td0 = 0.412; % 2015年相对于1850年(该预设过时且可能存在错误)

            %% Step2: 划分时间段以求拟合线
            % 二层EBM模型
            %% 多个不同时间段的变化
            startY = 0; endY = 84; interY = 1;
            % 计算gamma值以及其他的积分值等
            [integral_ts, integral_RnetTOA, dtd_ctrl, integral_td, delta_td, integral_delta_td, gammaa] = EBMLayer2_gamma(C_ref, Cd_ref, gblM_R_netTOA, dts_ctrl, b_td0, startY, endY, interY);
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Fig1 :plot delta Td
            xdata = 2015:1:2099;
            ydata = integral_delta_td';

            timeSer = [1985 1990 1995 2000 2005 2015 2025 2035 2045 2055 2065 2075 2085 2095 2105];

            yLabel = '\DeltaT_d (K)';
            xLabel = 'time line';

            unitePlot(xdata, ydata, xLabel, yLabel, timeSer, 'on');

            title({['ssp370 Model: ', mdlName]; ['T\_deep(2015): ', num2str(b_td0), '  T\_deep(2099): ', num2str(ydata(end))]; ['T(2015): ', num2str(bts0), '  T(2099): ', num2str(dts_ctrl(end)), ' (均相对于1850年)']})

            figName = [outPutPath, 'Td/ssp370_r1i1p1f1_', mdlName,'_',C_name,'.png'];

            saveas(gcf, figName)
            close gcf

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Fig2 :plot delta gamma
            xdata = integral_ts - integral_td;
            ydata = integral_delta_td .* Cd_ref;
            xdata = xdata(40:end);
            ydata = ydata(40:end);
            xLabel = '\Sigma(T-T_d)';
            yLabel = 'C_d\DeltaT_d';

            unitePlot(xdata, ydata, xLabel, yLabel, timeSer, 'off');
            hold on

            kb = polyfit(xdata, ydata, 1);
            k = kb(1);
            b = kb(2);
            yfitPlot = polyval(kb, 1:1000);
            plot(1:1000, yfitPlot, '-k', 'LineWidth', 3)

            title({['ssp370 Model: ', mdlName]; ['斜率: ', num2str(k), ' 截距: ', num2str(b), '取后44年结果拟合']})

            lgd = legend('data', 'Linear fit', 'Fontsize', 14, 'Location', 'northwest');
            lgd_inf = get(lgd);

            % 计算r2并标注
            yfit = polyval(kb, xdata);
            yresid = ydata - yfit; %将残差值计算为有符号数的向量：
            SSresid = sum(yresid.^2); %计算残差的平方并相加，以获得残差平方和：
            SStotal = (length(ydata) - 1) * var(ydata); %通过将观测次数减 1(自由度) 再乘以 y 的方差，计算 y 的总平方和：
            rsq = 1 - SSresid / SStotal; %计算r2
            text_ant = {['R^2= ', num2str(rsq)], ['\gamma=slope=', num2str(k)]};

            text(lgd_inf.Position(1) - 0.12, lgd_inf.Position(2), text_ant, 'Fontsize', 14, 'units', 'normalized');
            hold on

            disp([esmName{esmNum, 1}, ' ensemble is done!'])
            figName = [outPutPath, 'Gamma/ssp370_r1i1p1f1_', mdlName,'_',C_name,'.png'];
            saveas(gcf, figName)

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Fig3 :plot Tnc
            %% 2) 计算delta_dtnc, dtnc_ctrl, integral_dtnc(积分)  通过迭代的方法计算
            [integral_ERF, dtnc_ctrl, integral_dtnc, delta_dtnc, Delta_dtnc] = EBMLayer2_tnc(ERForcing, lamda_nonCld, gammaa, integral_td, bts0, C_ref, startY, endY, interY);

            xdata = 2015:1:2099;
            ydata = dtnc_ctrl';

            yLabel = '\T_{nc} (K)';
            xLabel = 'time line';
            kb = polyfit(xdata, ydata, 1);
            k = kb(1);
            b = kb(2);
            yfit = polyval(kb, xdata);

            unitePlot(xdata, ydata, xLabel, yLabel, timeSer, 'on');
            % old measure
            DeltaTsg_cld = DeltaTsg - Delta_dtnc;
       

            title({['ssp370 Model: ', mdlName]; ['斜率: ', num2str(k), ' T_{nc}(2015): ', num2str(ydata(1)), ' T_{nc}(2099): ', num2str(ydata(end))]; ['T(2015): ', num2str(bts0), '  T(2099): ', num2str(dts_ctrl(end)), ' (均相对于1850年)']; ['无云情况下温度变化:', num2str(Delta_dtnc), '云致温度变化:', num2str(DeltaTsg_cld)]})

            figName = [outPutPath, 'Tnc/ssp370_r1i1p1f1_', mdlName,'_',C_name,'.png'];
            saveas(gcf, figName)
            close gcf

        end

        disp([mdlName, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{exmNum}, ' era is done!'])
    disp(' ')

end

function [integral_ts, integral_RnetTOA, dtd_ctrl, integral_td, delta_td, integral_delta_td, gammaa] = EBMLayer2_gamma(C_ref, Cd_ref, gblM_R_netTOA, dts_ctrl, b_td0, startY, endY, interY)
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

function [] = unitePlot(xdata, ydata, xLabel, yLabel, timeSer, axtickSwich)

    char_timeSer = cellstr(string(timeSer));

    set(0, 'DefaultFigureVisible', 'off')
    ss = get(0, 'ScreenSize');
    h1 = figure('Position', [-1021 -743 843 545]);
    set(h1, 'Color', [1 1 1]);
    plot(xdata, ydata, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', '#67C1EE', 'LineWidth', 2, 'MarkerSize', 8)

    % set x axes
    if strcmp(axtickSwich, 'on')
        datetick('x', 'yyyy', 'keepticks'); % 用日期的格式显示横坐标
        xticks(timeSer)
        xticklabels(char_timeSer)
    end

    xlim([xdata(1) - 5 xdata(end) + 5])
    ylim([min(ydata) - 0.05 * (max(ydata) - min(ydata)) max(ydata) + 0.05 * (max(ydata) - min(ydata))])
    xlabel(xLabel)
    ylabel(yLabel)

    ax = gca;
    % ax.TickLength = [0.02 0.01]; %刻度线长度      set(gca,'ticklength', [0.02 0.01]);
    ax.XColor = 'k'; ax.YColor = 'k'; % 设置刻度线颜色
    ax.FontSize = 13;

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
