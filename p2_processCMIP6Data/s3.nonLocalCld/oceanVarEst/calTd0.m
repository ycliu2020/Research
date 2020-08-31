%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-08-25 20:02:03
% LastEditTime : 2020-08-26 10:41:30
% LastEditors  : LYC
% Description  : 根据NCEI的数据估算2015年初始深海温度较1850年变化了多少
% FilePath     : /code/p2_processCMIP6Data/s3.nonLocalCld/oceanVarEst/calTd0.m
%  
%%---------------------------------------------------------
clear; clc; tic;
heatCPath='/data1/liuyincheng/Observe-rawdata/NCEI/Ocean Heat/';
heatC_700m=ncread([heatCPath,'heat_content_anomaly_0-700_yearly.nc'],'yearl_h22_WO'); 
heatC_2000m=ncread([heatCPath,'heat_content_anomaly_0-2000_yearly.nc'],'yearl_h22_WO');
heatC_700m_2005=heatC_700m(end-14:end-1);
heatC_2000m_2005=heatC_2000m(end-14:end-1);

% 假设0-2000m的热容在2005年之前按照0-500m的线性斜率变化, 2005年后按观测数据变化
heatC_700m_fix=heatC_700m(1:end-15)+(heatC_2000m_2005(1)-heatC_700m(end-14));
heatC_2000m_1955Pre=heatC_700m_fix;
heatC_2000m_1955Pre(end+1:end+14)=heatC_2000m_2005;

timeHC=1955:2019-1;
timeHC_obs=1955:2019-15;
timeHC_pre=1850:2015;

k_2000=polyfit(timeHC_obs',heatC_700m_fix,1);
heatC_2000m_1850Pre = polyval(k_2000,timeHC_pre);
heatC_2000m_1850Pre(end-13:end)=heatC_2000m_2005;

figure
plot(timeHC_pre,heatC_2000m_1850Pre,'o') 
figure
plot(timeHC,heatC_2000m_1955Pre,'o') 
% timeHC_pre=1850:2015;
% 
% k_700(1)
% heatC_2000fit = polyval(k,timeHC_pre); 
% % cal test

mean2015_hc2000m=mean(heatC_2000m_1850Pre(find(timeHC_pre==2014):end));
Td0=(mean2015_hc2000m-heatC_2000m_1850Pre(timeHC_pre==1850))*1e22/(365.25*86400*105.5285*3.619e8*1e3^2);