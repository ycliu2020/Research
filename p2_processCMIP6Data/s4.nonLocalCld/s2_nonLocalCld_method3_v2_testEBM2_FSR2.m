%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-08-09 14:58:33
% LastEditTime : 2020-09-23 09:54:14
% LastEditors  : LYC
% Description  : 这里和1的区别主要采用了离散的计算方式计算gamma项  修正计算非云致温度变化实验可行性分析 feasibility study report
% FilePath     : /code/p2_processCMIP6Data/s4.nonLocalCld/s2_nonLocalCld_method3_v2_testEBM2_FSR2.m
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
    for mdlNum = 1:1%length(level.model2)
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

            % load rad
            load([dradEffectPath, 'dR_cloud_toa.mat'])% cloud radRffect

            % load surface temperature
            load([varsPath, 'ts.mat'])% ts
            load([dvarsPath, 'dts.mat'])%dts
            dts = autoRegrid3(lon_k, lat_k, time.date, dts, lon_f, lat_f, time.date);
            ts = autoRegrid3(lon_k, lat_k, time.date, ts, lon_f, lat_f, time.date);

            % load radiation
            load([dradEffectPath, 'dradEfect_toa_cld.mat'])% 'totalEffect', 'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', 'taEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect', 'mainEffect'
            load([dradEffectPath, 'real_dradEfect.mat'])% 'dR_allsky', 'l_rad', 's_rad', 'dR_clr', 'readme_realradEfect'
            load([varsPath, 'rlut.mat'])% toa_outgoing_longwave_flux
            load([varsPath, 'rsut.mat'])% toa_outgoing_shortwave_flux
            load([varsPath, 'rsdt.mat'])% toa_incoming_shortwave_flux
            dR_nonCloud = totalEffect;
            dR_netTOA = squeeze(dR_allsky(:, :, :, 2));

            netTOA = rsdt - rlut - rsut; % net radiation in toa (lonxlatxmonth): W/m2
            netTOA = autoRegrid3(lon_k, lat_k, time.date, netTOA, lon_f, lat_f, time.date);

            % 将月时间序列转换为年时间序列(lon,lat,1020)→(lon,lat,1020/12)并进行纬向加权平均
            gblM_ts = monthlyToYearlyMean(ts, lat_f);
            gblM_dts = monthlyToYearlyMean(dts, lat_f);

            gblM_dR_netTOA = monthlyToYearlyMean(dR_netTOA, lat_f);
            gblM_R_netTOA = monthlyToYearlyMean(netTOA, lat_f);

            gblM_dR_cloud_toa = monthlyToYearlyMean(dR_cloud_toa, lat_f);
            gblM_dR_nonCloud = monthlyToYearlyMean(dR_nonCloud, lat_f);

            % 海温基态值(control climate), 1850年为基态, 所以netTOA=0
            % 以1850为基态
            T_ref1850 = 13.9 + 273.15 - 0.373; %K, see https://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/download.html, 13.9 is mean of 1961-1990
            % 预定观测的bts0(初始表面温度)
            T_ref2015Series = [0.763 0.797 0.677 0.587 0.736 0.864];
            T_ref2015Series = T_ref2015Series + 13.9 + 273.15;
            T_ref2015Serie = mean(T_ref2015Series);
            bts0 = T_ref2015Serie - T_ref1850; % 用观测的
            % 计算基于基态的时间序列
            ts15fix=gblM_ts(1)-T_ref2015Serie;
            dts_ctrl = gblM_ts -ts15fix- T_ref1850;
            % 预定初始深海温度(见oceanVarEst文件夹)
            b_td0 = 0.412; % K 相对于1850年
            
            %% Step2: 划分时间段以求拟合线
            % 二层EBM模型
            % model parameters
            D_ref = 55; Dd_ref = 2768; % D_ref = 77; Dd_ref = 1105; in Geoffroy and D_ref = 55; Dd_ref = 2768; in Jiménez-de-la-Cuesta,
            cp_water = 4180; % heat capacity J/kg/K
            rho_water = 1030; % desnsity of saltwater kg/m3
            f0 = 0.7; % The proportion of ocean cover the earth surface
            % cal reference C and C0
            C_ref = (D_ref * rho_water * cp_water * f0) / (86400 * 365.25); % C_ref(upper-ocean) mean the result of formula (22) in reference Geoffroy
            Cd_ref = (Dd_ref * rho_water * cp_water * f0) / (86400 * 365.25); % Cd_ref(deep-layer ocean) mean the result of formula (22) in reference Geoffroy
            
            %% 多个不同时间段的变化
            startY = 0; endY = 84; interY = 1;
            timeExamp = startY + interY * 2 - 1:interY:endY;
            sum_timeExamp = length(timeExamp); % 总共有几组时间段
        
            delta_ts = zeros(sum_timeExamp+1, 1);% Delta_ts(1)=0
            integral_delta_ts = zeros(sum_timeExamp+1, 1);
            integral_ts = zeros(sum_timeExamp+1, 1);

            delta_td = zeros(sum_timeExamp+1, 1);
            integral_delta_td = zeros(sum_timeExamp+1, 1);
            dtd_ctrl=zeros(sum_timeExamp+1, 1);

            integral_RnetTOA = zeros(sum_timeExamp+1, 1); % integral_RnetTOA=sum(RnetTOA)
            
            % 设定初始值
            delta_td(1) = 0;
            integral_delta_td(1) = 0;
            dtd_ctrl(1)=b_td0;

            delta_ts(1) = 0;
            integral_delta_ts=dts_ctrl-dts_ctrl(1);
            integral_ts(1) = dts_ctrl(1);

            integral_RnetTOA(1) = gblM_R_netTOA(1);

            % 计算每个格点的值
            for gridNum = 2:sum_timeExamp+1 % 循环每一个格点
    
                Seri_RnetTOA = gblM_R_netTOA(1:timeExamp(gridNum-1) - startY+1);
                Seri_dts = dts_ctrl(1:timeExamp(gridNum-1) - startY+1);
        
                integral_ts(gridNum) = sum(Seri_dts);
                integral_RnetTOA(gridNum) = sum(Seri_RnetTOA);
                delta_td(gridNum)=(integral_RnetTOA(gridNum)-C_ref*integral_delta_ts(gridNum))/Cd_ref-integral_delta_td(gridNum-1);
                integral_delta_td(gridNum)=integral_delta_td(gridNum-1)+delta_td(gridNum);
                dtd_ctrl=dtd_ctrl(1)+integral_delta_td;
            end
        
        
            % 计算Td积分
            integral_td = zeros(sum_timeExamp+1, 1);
            integral_td(1) = b_td0;
            for gridNum = 2:sum_timeExamp+1%循环每一组
                Seri_dtd = dtd_ctrl(1:timeExamp(gridNum-1) - startY + 1);
                integral_td(gridNum) = sum(Seri_dtd);
            end
            
            % 计算gamma (function 2, see document for details)
            xdata = integral_ts-integral_td;
            ydata = integral_delta_td.* Cd_ref;
            % gamma1=ydata./xdata;
            xdata=xdata(1:end);
            ydata=ydata(1:end);
            yLabel = 'x';
            xLabel = 'y';
        
            % 线性拟合
            [k_fit, b_fit] = plotLinFit(xdata, ydata, xLabel, yLabel, mdlName, []);
        
            % % 线性拟合(过零点)
            % [k_fit, b_fit] = plotLinFit_point0(xdata, ydata, xLabel, yLabel, mdlName);
            gammaa = k_fit;
            resiual = b_fit;
            title({['ssp370 Model: ', mdlName]; ['gamma: ', num2str(gammaa), '  Res(认为截距为误差): ', num2str(resiual)]; ['T\_deep(2015): ', num2str(b_td0), '   T(2015): ', num2str(bts0), '(均相对于1850年)']; ['T-T\_deep=', num2str(bts0 - b_td0)]})
            figName = ['/home/liuyc/Research/P02.Ts_change_research/figure/02.cmip6Result/1.6/nonLocalCld3/Cktest/ssp370_r1i1p1f1_', mdlName, '_Jim.png'];
            saveas(gcf, figName)
        
            % [Seri_k_dts, Seri_b_dts, Seri_k_dtd, Seri_b_dRnetTOA, gammaa] = EBMLayer2_gammaTest(C_ref, Cd_ref, gblM_R_netTOA, dts_ctrl, startY, endY, interY, mdlName);
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

