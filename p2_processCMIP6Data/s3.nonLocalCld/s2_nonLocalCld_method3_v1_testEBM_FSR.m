%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-08-09 14:58:33
% LastEditTime : 2020-09-01 16:12:07
% LastEditors  : LYC
% Description  : 修正计算非云致温度变化实验可行性分析 feasibility study report
% FilePath     : /code/p2_processCMIP6Data/s3.nonLocalCld/s2_nonLocalCld_method3_testEBM_FSR.m
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
            load([varsPath, 'rlut.mat'])% toa_outgoing_longwave_flux
            load([varsPath, 'rsut.mat'])% toa_outgoing_shortwave_flux
            load([varsPath, 'rsdt.mat'])% toa_incoming_shortwave_flux
            dR_nonCloud = totalEffect;
            dR_netTOA = squeeze(dR_allsky(:, :, :, 2));

            %regrid to 144x72(unite grids)
            netTOA = rsdt - rlut - rsut; % net radiation in toa (lonxlatxmonth): W/m2
            netTOA = autoRegrid3(lon_k, lat_k, time.date, netTOA, lon_f, lat_f, time.date);

            %% Step1: 将月时间序列转换为年时间序列(lon,lat,1020)→(lon,lat,1020/12)并进行纬向加权平均
            gblM_ts = monthlyToYearlyMean(ts, lat_f);
            gblM_dts = monthlyToYearlyMean(dts, lat_f);

            gblM_dR_netTOA = monthlyToYearlyMean(dR_netTOA, lat_f);
            gblM_R_netTOA = monthlyToYearlyMean(netTOA, lat_f);

            gblM_dR_cloud_toa = monthlyToYearlyMean(dR_cloud_toa, lat_f);
            gblM_dR_nonCloud = monthlyToYearlyMean(dR_nonCloud, lat_f);

            % %% 基态值(control climate)
            % 以1850年的数值为基态
            T_ref1850 = 13.9 + 273.15 - 0.373; %K, see https://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/download.html, 13.9 is mean of 1961-1990
            % 预定观测的bts0(初始表面温度)
            T_ref2015Series = [0.763 0.797 0.677 0.587 0.736 0.864];
            T_ref2015Series = T_ref2015Series + 13.9 + 273.15;
            bts0 = mean((T_ref2015Series - T_ref1850)); % 用观测的
            % 计算基于基态的时间序列
            dts_ctrl = gblM_dts - T_ref1850;
            % 预定初始深海温度(见oceanVarEst文件夹)
            b_td0 = 0.412; % K 相对于1850年
            % 1850年为基态, 所以netTOA=0
            
            %% Step2: 划分时间段以求拟合线
            % % 一层EBM模型(拟合效果不佳)
            % startY=1;endY=80;interY=5;
            % [k1,C1, x1, y1,Integral_ts,Integral_RnetTOA,Delta_ts] = EBMLayer1( gblM_dR_netTOA, gblM_dts, startY, endY, interY, mdlName);

            % 二层EBM模型
            % model parameters
            D_ref = 55; Dd_ref = 2768; %;% D_ref = 77; Dd_ref = 1105; in Geoffroy and D_ref = 55; Dd_ref = 2768; in Jiménez-de-la-Cuesta,
            cp_water = 4180; % heat capacity J/kg/K
            rho_water = 1030; % desnsity of saltwater kg/m3
            f0 = 0.7; % The proportion of ocean cover the earth surface
            % cal reference C and C0
            C_ref = (D_ref * rho_water * cp_water * f0) / (86400 * 365.25); % C_ref(upper-ocean) mean the result of formula (22) in reference Geoffroy
            Cd_ref = (Dd_ref * rho_water * cp_water * f0) / (86400 * 365.25); % Cd_ref(deep-layer ocean) mean the result of formula (22) in reference Geoffroy

            startY = 0; endY = 84; interY = 1;
            [Seri_k_dts, Seri_b_dts, Seri_k_dtd, Seri_b_dRnetTOA, gammaa] = EBMLayer2_gammaTest(C_ref, Cd_ref, gblM_R_netTOA, dts_ctrl, startY, endY, interY, mdlName);
            % [Seri_k_dts, Seri_b_dts, Seri_k_dtd, gammaa] = EBMLayer2_gammaTest2(C_ref, Cd_ref, gblM_R_netTOA, dts_ctrl, bts0, b_td0, startY, endY, interY, mdlName);

            disp([esmName{esmNum, 1}, ' ensemble is done!'])
        end

        disp([mdlName, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{exmNum}, ' era is done!'])
    disp(' ')

end

t = toc; disp(t)

function [Seri_k_dts, Seri_b_dts, Seri_k_dtd, Seri_b_dRnetTOA, gammaa] = EBMLayer2_gammaTest(C_ref, Cd_ref, dR_netTOA, dts, bts0, b_td0, startT, endT, interT, modelName, figName)

    switch nargin
        case 10
            figName = [];
    end

    timeExamp = startT + interT * 6 - 1:interT:endT;
    diffTimeSum = length(timeExamp); % 总共有几组时间段

    Delta_ts = zeros(diffTimeSum, 1);
    Seri_k_dtd = zeros(diffTimeSum, 1);

    Seri_k_dts = zeros(diffTimeSum, 1);
    Seri_b_dts = zeros(diffTimeSum, 1);

    Seri_k_dRnetTOA = zeros(diffTimeSum, 1);
    Seri_b_dRnetTOA = zeros(diffTimeSum, 1);
    p_Plus = zeros(diffTimeSum, 1);

    % 计算k_dtd (function 1)
    for diffTime = 1:diffTimeSum%循环每一组
        timeSeries = startT:timeExamp(diffTime);
        dtime = length(timeSeries);

        Seri_dRnetTOA_ctrl = dR_netTOA(1:timeExamp(diffTime) - startT + 1);
        Seri_dts_ctrl = dts(1:timeExamp(diffTime) - startT + 1);
        Seri_dts = dts(1:timeExamp(diffTime) - startT + 1);

        %% cal Delta ts
        k_dts = polyfit(timeSeries', Seri_dts_ctrl, 1);
        Delta_ts(diffTime) = k_dts(1) * dtime;
        Seri_k_dts(diffTime) = k_dts(1);
        Seri_b_dts(diffTime) = bts0;

        k_dRnetTOA = polyfit(timeSeries', Seri_dRnetTOA_ctrl, 1);
        Seri_k_dRnetTOA(diffTime) = k_dRnetTOA(1);
        Seri_b_dRnetTOA(diffTime) = k_dRnetTOA(2);
        p_Plus(diffTime) = timeExamp(diffTime) + startT; %(t1+t2)
        Seri_k_dtd(diffTime) = (0.5 * Seri_k_dRnetTOA(diffTime) * p_Plus(diffTime) + Seri_b_dRnetTOA(diffTime) - C_ref * Seri_k_dts(diffTime)) / Cd_ref;
    end

    % 计算gamma (function 2, see document for details)
    xdata = 0.5 * (Seri_k_dts - Seri_k_dtd) .* p_Plus + Seri_b_dts - b_td0;
    ydata = Seri_k_dtd .* Cd_ref;
    % gamma1=ydata./xdata;

    yLabel = 'Cd*ktd';
    xLabel = '0.5*(k\_t-k\_td)(t1+t2)+b\_t-b\_td';

    % 线性拟合
    [k_fit, b_fit] = plotLinFit(xdata, ydata, xLabel, yLabel, modelName, figName);

    % % 线性拟合(过零点)
    % [k_fit, b_fit] = plotLinFit_point0(xdata, ydata, xLabel, yLabel, modelName);
    gammaa = k_fit;
    resiual = b_fit;
    title({['ssp370 Model: ', modelName]; ['gamma: ', num2str(gammaa), '  Res(认为截距为误差): ', num2str(resiual)]; ['T\_deep(2015): ', num2str(b_td0), '   T(2015): ', num2str(bts0), '(均相对于1850年)']; ['T-T\_deep=', num2str(bts0 - b_td0)]})
    figName = ['/home/liuyc/Research/P02.Ts_change_research/figure/02.cmip6Result/1.6/nonLocalCld2_fix/Cktest/ssp370_r1i1p1f1_', modelName, '_Jim.png'];
    saveas(gcf, figName)

end




function [Seri_k_dts, Seri_b_dts, Seri_k_dtd gammaa] = EBMLayer2_gammaTest2(C_ref, Cd_ref, dR_netTOA, dts, bts0, b_td0, startT, endT, interT, modelName, figName)

    switch nargin
        case 10
            figName = [];
    end

    timeExamp = startT + interT * 6 - 1:interT:endT;
    diffTimeSum = length(timeExamp); % 总共有几组时间段

    Delta_ts = zeros(diffTimeSum, 1);
    Seri_k_dtd = zeros(diffTimeSum, 1);

    Seri_k_dts = zeros(diffTimeSum, 1);
    Seri_b_dts = zeros(diffTimeSum, 1);

    p_Plus = zeros(diffTimeSum, 1);

    % 计算k_dtd (function 1)
    for diffTime = 1:diffTimeSum%循环每一组
        timeSeries = startT:timeExamp(diffTime);
        dtime = length(timeSeries);

        Seri_dRnetTOA_ctrl = dR_netTOA(1:timeExamp(diffTime) - startT + 1);
        Seri_dts_ctrl = dts(1:timeExamp(diffTime) - startT + 1);
        Seri_dts = dts(1:timeExamp(diffTime) - startT + 1);

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

    yLabel = 'Cd*ktd';
    xLabel = '0.5*(k\_t-k\_td)(t1+t2)+b\_t-b\_td';

    % 线性拟合
    [k_fit, b_fit] = plotLinFit(xdata, ydata, xLabel, yLabel, modelName, figName);

    % % 线性拟合(过零点)
    % [k_fit, b_fit] = plotLinFit_point0(xdata, ydata, xLabel, yLabel, modelName);
    gammaa = k_fit;
    resiual = b_fit;
    title({['ssp370 Model: ', modelName]; ['gamma: ', num2str(gammaa), '  Res(认为截距为误差): ', num2str(resiual)]; ['T\_deep(2015): ', num2str(b_td0), '   T(2015): ', num2str(bts0), '(均相对于1850年)']; ['T-T\_deep=', num2str(bts0 - b_td0)]})
    % figName = ['/home/liuyc/Research/P02.Ts_change_research/figure/02.cmip6Result/1.6/nonLocalCld2_fix/Cktest/ssp370_r1i1p1f1_', modelName, '_Jim.png'];
    % saveas(gcf, figName)

end

function [k, C, x, y, Integral_ts, Integral_RnetTOA, Delta_ts] = EBMLayer1Test(dR_netTOA, dts, startT, endT, interT, modelName, figName)
    %
    % description.
    % dts: d均为全球纬向平均
    % startT: 开始的年
    % endT: 结束的年
    % interT: 中间的间隔
    %
    switch nargin
        case 6
            figName = [];
    end

    timeExamp = startT + interT * 6 - 1:interT:endT;
    diffTimeSum = length(timeExamp); % 总共有几组时间段

    Delta_ts = zeros(diffTimeSum, 1);
    Integral_RnetTOA = zeros(diffTimeSum, 1);
    Integral_ts = zeros(diffTimeSum, 1);
    Seri_k_dts = zeros(diffTimeSum, 1);
    Seri_b_dts = zeros(diffTimeSum, 1);
    Seri_k_dRnetTOA = zeros(diffTimeSum, 1);
    Seri_b_dRnetTOA = zeros(diffTimeSum, 1);
    numerator = zeros(diffTimeSum, 1);
    p_Plus = zeros(diffTimeSum, 1);
    % 计算基于基态的变量
    ts_ctrl = mean(dts(1:10));
    dts_ctrl = dts - ts_ctrl;
    RnetTOA_ctrl = mean(dR_netTOA(1:10));
    dRnetTOA_ctrl = dR_netTOA - RnetTOA_ctrl;

    for diffTime = 1:diffTimeSum%循环每一组
        timeSeries = startT:timeExamp(diffTime);
        Seri_dRnetTOA_ctrl = dRnetTOA_ctrl(startT:timeExamp(diffTime));
        Seri_dts_ctrl = dts_ctrl(startT:timeExamp(diffTime));
        Seri_dts = dts(startT:timeExamp(diffTime));

        dtime = length(timeSeries);

        %% cal Delta ts
        k_dts = polyfit(timeSeries', Seri_dts_ctrl, 1);
        Delta_ts(diffTime) = k_dts(1) * dtime;
        Seri_k_dts(diffTime) = k_dts(1);
        Seri_b_dts(diffTime) = k_dts(2);
        k_dRnetTOA = polyfit(timeSeries', Seri_dRnetTOA_ctrl, 1);
        Seri_k_dRnetTOA(diffTime) = k_dRnetTOA(1);
        Seri_b_dRnetTOA(diffTime) = k_dRnetTOA(2);
        %% 计算积分
        % old method
        Integral_RnetTOA(diffTime) = sum(Seri_dRnetTOA_ctrl);
        Integral_ts(diffTime) = sum(Seri_dts_ctrl);
        Integral_ts1(diffTime) = sum(Seri_dts);
        % new method
        p_Plus(diffTime) = timeExamp(diffTime) + startT;
    end

    % polyfit y and x
    %  plot function Newmethod
    y = (Seri_k_dRnetTOA .* p_Plus + 2 * Seri_b_dRnetTOA) ./ (Seri_k_dts .* p_Plus + 2 * Seri_b_dts);
    x = 2 * Seri_k_dts ./ (Seri_k_dts .* p_Plus + 2 * Seri_b_dts);
    xLabel = '2k\_ts/(k\_ts(t1+t2)+2b\_ts)';
    yLabel = '(k\_netTOA(t1+t2)+2b\_netTOA)/(k\_ts(t1+t2)+2b\_ts)';
    [k_fit, b_fit] = plotLinFit(x, y, xLabel, yLabel, modelName, figName);
    C = k_fit;
    k = b_fit;

    % y = Integral_RnetTOA ./ Integral_ts;
    % x = Delta_ts ./ Integral_ts;
    % kb = polyfit(x, y, 1);
    % C = kb(1);
    % k = kb(2);
    % yfit = polyval(kb, x);
    % figure;
    % plot(x, y, 'o', x, yfit, '-')
    % xlabel('DeltaT/T(t)dt积分')
    % ylabel('RnetTOA(t)dt积分/T(t)dt积分')
    % title({['ssp370 Model: ',modelName];['C= ', num2str(C), ' k= ', num2str(k)]})
    % % 计算r2并标注
    % yresid = y - yfit; %将残差值计算为有符号数的向量：
    % SSresid = sum(yresid.^2); %计算残差的平方并相加，以获得残差平方和：
    % SStotal = (length(y) - 1) * var(y); %通过将观测次数减 1(自由度) 再乘以 y 的方差，计算 y 的总平方和：
    % rsq = 1 - SSresid / SStotal; %计算r2
    % text(0.75, 0.8, ['r2= ', num2str(rsq)], 'units', 'normalized');
    % hold on
    % legend('data', 'linear fit')
    % figName=['/home/liuyc/Research/P02.Ts_change_research/figure/02.cmip6Result/1.6/nonLocalCld2_fix/Cktest/ssp370_r1i1p1f1_',modelName,'.png'];
    % saveas(gcf, figName)
    % % close gcf

    % % teacher test
    % x = Integral_RnetTOA;
    % y = Delta_ts;

    % kb = polyfit(x, y, 1);
    % C = kb(1);
    % k = kb(2);
    % yfit = polyval(kb, x);
    % figure;
    % plot(x, y, 'o', x, yfit, '-')
    % xlabel('netTOA*dt积分')
    % ylabel('\DeltaT(相对于基态的变化)')
    % title({['ssp370 Model: ',modelName];['斜率= ', num2str(C), ' 截距= ', num2str(k)]})
    % % title({['ssp370 Model: ',modelName];['C= ', num2str(C), ' k= ', num2str(k)]})
    % % 计算r2并标注
    % yresid = y - yfit; %将残差值计算为有符号数的向量：
    % SSresid = sum(yresid.^2); %计算残差的平方并相加，以获得残差平方和：
    % SStotal = (length(y) - 1) * var(y); %通过将观测次数减 1(自由度) 再乘以 y 的方差，计算 y 的总平方和：
    % rsq = 1 - SSresid / SStotal; %计算r2
    % text(0.75, 0.8, ['r2= ', num2str(rsq)], 'units', 'normalized');
    % hold on
    % legend('data', 'linear fit')
    % figName=['/home/liuyc/Research/P02.Ts_change_research/figure/02.cmip6Result/1.6/nonLocalCld2_fix/RTtest/ssp370_r1i1p1f1_',modelName,'.png'];
    % % saveas(gcf, figName)
    % % close gcf

end
