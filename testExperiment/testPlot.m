%%---------------------------------------------------------
% Author       : LYC
% Date         : 2021-03-21 15:11:13
% LastEditors  : Please set LastEditors
% Description  : 
% FilePath     : /code/testExperiment/testPlot.m
% symbol_custom_string_obkoro1:  
%%---------------------------------------------------------
clear; clc;
inputPath='/data1/liuyincheng/CMIP6-rawdata/ssp126/zos_Omon_AWI-CM-1-1-MR_ssp126_r1i1p1f1_gn_201501-202012.nc';
% 首先添加m_map工具箱的路径

% read data
time=ncread(inputPath,'time');
lat=ncread(inputPath,'lat');
lon=ncread(inputPath,'lon');
zos=ncread(inputPath,'zos');

lat=reshape(lat,3389,245);
lon=reshape(lon,3389,245);
zos=reshape(zos,3389,245,72);
% plot
% select 2015-01 as exmaple
% rangeL=1000;
% rangeR=10000;
% plot_zos=squeeze(zos(rangeL:rangeR,1));
% plot_lat=squeeze(lat(rangeL:rangeR,1));
% plot_lon=squeeze(lon(rangeL:rangeR,1));
figure
plot_zos=squeeze(zos(:,:,1));

lon1=[0 360];lat1=[-75 75];%画布经纬度范围lat1=[-75.5 77.5]
m_proj('Mercator','lon',lon1, 'lat',lat1);%投影坐标等预设参数Mercator,Equidistant cylindrical,Miller Cylindrical
m_grid('linestyle','none','tickdir','out','yaxislocation','left','fontsize',10,'color','k');%边框坐标轴选项


m_pcolor(lon,lat,plot_zos);
m_coast;       
% m_pcolor(plot_lat,plot_lon,plot_zos');       
colormap(m_colmap('jet','step',10));
h=colorbar('northoutside');
