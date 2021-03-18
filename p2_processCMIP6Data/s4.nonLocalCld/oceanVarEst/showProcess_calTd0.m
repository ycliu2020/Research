%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-08-25 20:02:03
% LastEditTime : 2021-03-11 19:51:24
% LastEditors  : Please set LastEditors
% Description  : 根据NCEI的数据估算2015年初始深海温度较1850年变化了多少
% FilePath     : /code/p2_processCMIP6Data/s4.nonLocalCld/oceanVarEst/showProcess_calTd0.m
%
%%---------------------------------------------------------

clear; clc; tic;

% read Data
heatCPath = '/data1/liuyincheng/Observe-rawdata/NCEI/Ocean Heat/';
heatC_700m = ncread([heatCPath, 'heat_content_anomaly_0-700_yearly.nc'], 'yearl_h22_WO');
heatC_700m_pentad = ncread([heatCPath, 'heat_content_anomaly_0-700_pentad.nc'], 'pent_h22_WO');

heatC_2000m = ncread([heatCPath, 'heat_content_anomaly_0-2000_yearly.nc'], 'yearl_h22_WO');
heatC_2000m_pentad = ncread([heatCPath, 'heat_content_anomaly_0-2000_pentad.nc'], 'pent_h22_WO');
heatC_700m = heatC_700m(1:end - 1);
heatC_2000m = heatC_2000m(1:end - 1);

heatC_700m_1955 = heatC_700m;
heatC_700m_1955_2004 = heatC_700m(1:end - 14);
heatC_700m_2005 = heatC_700m(end - 13:end);
heatC_2000m_2005 = heatC_2000m(end - 13:end);

% time line
timeHC_2005 = 2005:2019 - 1;
timeHC_1955 = 1955:2019 - 1;
timeHC_1955_2004 = 1955:2019 - 15;
timeHC_1850 = 1850:2019 - 1;

%% test The correlation of 700m and 2000m
figure
plot(timeHC_2005, heatC_2000m_2005, 'ro')
hold on
plot(timeHC_1955, heatC_700m, 'bo')
title('Correlation test of HC 2000m and 700m')

lgd_size = 12;
lgd = legend('HC 2000m', 'HC 700m', 'FontName', 'Microsoft YaHei', 'Fontsize', lgd_size, 'Location', 'northwest');
lgd_inf = get(lgd);
cc0 = corrcoef(heatC_2000m_2005, heatC_700m_2005, 'Rows', 'complete');
text_ant = {['cc= ', num2str(cc0(1, 2))]};
text(lgd_inf.Position(1) - 0.13, lgd_inf.Position(2), text_ant, 'FontName', 'Microsoft YaHei', 'Fontsize', lgd_size, 'units', 'normalized');

%% test 那种拟合线最适合 700m, 由相关性可知, 同一类型的拟合线也适用于2000m
% method 1 liner fit
ydata = heatC_700m_1955;
xdata = timeHC_1955;
k_700_liner = polyfit(xdata', ydata, 1);
heatC_700m_linFit = polyval(k_700_liner, xdata);
figure
plot(timeHC_2005, heatC_2000m_2005, 'ro')
hold on;
plot(xdata, ydata, 'bo')
hold on;
plot(xdata, heatC_700m_linFit, 'b-')
hold on;

% method 2 Polynomial curve fitting
ydata = heatC_700m_1955;
xdata = timeHC_1955;
k_700_mul = polyfit(xdata', log(ydata + 10), 1);
heatC_700m_mulFit = exp(polyval(k_700_mul, xdata)) - 10;
plot(xdata, heatC_700m_mulFit, 'b*')

% legend set
lgd_size = 12;
lgd = legend('HC 2000m', 'HC 700m', 'Liner fit', 'multiple fit', 'FontName', 'Microsoft YaHei', 'Fontsize', lgd_size, 'Location', 'northwest');
lgd_inf = get(lgd);

% 计算r2并标注
ydata = heatC_700m;
yfit = heatC_700m_linFit';
rsqLiner = calRSquare(ydata, yfit);

ydata = heatC_700m;
yfit = heatC_700m_mulFit';
rsqMul = calRSquare(ydata, yfit);
text_ant = {['R^2(LinerFit)= ', num2str(rsqLiner)], ['R^2(MulFit)= ', num2str(rsqMul)]};
text(lgd_inf.Position(1) - 0.12, lgd_inf.Position(2) - 0.05, text_ant, 'FontName', 'Microsoft YaHei', 'Fontsize', lgd_size, 'units', 'normalized');
hold on
title('Selection of fitting method')

%% 上述实验表明指数估计最为准确, 采用指数拟合模拟HC 2000m 的情况
ydata = heatC_2000m_2005;
xdata = timeHC_2005;
k_2000_mul = polyfit(xdata', log(ydata + 20), 1);

xdata_predict = timeHC_1850;
heatC_2000m_mulFit_predict = exp(polyval(k_2000_mul, xdata_predict)) - 20;
heatC_2000m_mulFit = heatC_2000m_mulFit_predict(end - 13:end);
figure
plot(xdata, ydata, 'ro')
hold on
plot(xdata_predict, heatC_2000m_mulFit_predict, 'r-.')

% legend set
lgd_size = 12;
lgd = legend('HC 2000m', 'multiple fit', 'FontName', 'Microsoft YaHei', 'Fontsize', lgd_size, 'Location', 'northwest');
lgd_inf = get(lgd);

% 计算r2并标注
ydata = heatC_2000m;
yfit = heatC_2000m_mulFit';
rsqLiner = calRSquare(ydata, yfit);

text_ant = {['R^2(MulFit)= ', num2str(rsqLiner)]}; %,['R^2(MulFit)= ', num2str(rsqMul)]};
text(lgd_inf.Position(1) - 0.12, lgd_inf.Position(2) - 0.05, text_ant, 'FontName', 'Microsoft YaHei', 'Fontsize', lgd_size, 'units', 'normalized');
hold on
title('Fitting HC2000m')

% choose exact guy's parameter
[C_ref_Jim, Cd_ref_Jim, C_ref_Geo, Cd_ref_Geo] = calCCd();
C_ref = C_ref_Geo;
Cd_ref = Cd_ref_Jim;
%% computation of C and Cd
heatC_2000_2015year = heatC_2000m(timeHC_2005 == 2015);
heatC_2000_1850year = heatC_2000m_mulFit_predict(1);
Td0 = (heatC_2000_2015year - heatC_2000_1850year) * 1e22 / (365.25 * 86400 * Cd_ref * 3.619e8 * 1e3^2);
disp(['Td(2015)=',num2str(Td0)])

function rsq = calRSquare(ydata, yfit)
    yresid = ydata - yfit; %将残差值计算为有符号数的向量：
    SSresid = sum(yresid.^2); %计算残差的平方并相加，以获得残差平方和：
    SStotal = (length(ydata) - 1) * var(ydata); %通过将观测次数减 1(自由度) 再乘以 y 的方差，计算 y 的总平方和：
    rsq = 1 - SSresid / SStotal; %计算r2
    rsq = roundn(rsq, -3);
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
