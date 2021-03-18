%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-08-25 20:02:03
% LastEditTime : 2021-03-11 20:04:35
% LastEditors  : Please set LastEditors
% Description  : 根据NCEI的数据估算2015年初始深海温度较1850年变化了多少
% FilePath     : /code/p2_processCMIP6Data/s4.nonLocalCld/oceanVarEst/calTd0.m
%
%%---------------------------------------------------------
clear; clc; tic;
% test
[C_ref_Jim, Cd_ref_Jim, C_ref_Geo, Cd_ref_Geo] = calCCd();
Cd_ref = Cd_ref_Jim;
[Td0] = calTd0F(Cd_ref);
Cd_ref = Cd_ref_Geo;
[Td0] = calTd0F(Cd_ref);

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