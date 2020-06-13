%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-13 20:22:31
% LastEditTime : 2020-06-13 21:03:18
% LastEditors  : LYC
% Description  : lamda_cloud mean cloud climate feedback parameter, Data souce: Zelinka MD, et al. (2020) Causes of higher climate sensitivity in CMIP6 models. Geophysical Research Letters 47(1):e2019GL085782.
% FilePath     : /Research/p2_processCMIP6Data/s1.modelDataProcess/Info_lamda_cloud.m
%
%%---------------------------------------------------------
% model name
lamda_cloud.modelName = {'BCC-CSM2-MR', 'BCC-ESM1', 'CAMS-CSM1-0', 'CESM2', 'CESM2-WACCM', 'CNRM-CM6-1', ...
    'CNRM-CM6-1-HR', 'CNRM-ESM2-1', 'CanESM5', 'E3SM-1-0', 'EC-Earth3', 'EC-Earth3-Veg', 'FGOALS-f3-L', ...
    'GFDL-CM4', 'GISS-E2-1-G', 'GISS-E2-1-H', 'HadGEM3-GC31-LL', 'INM-CM4-8', 'IPSL-CM6A-LR', 'MIROC-ES2L', ...
    'MIROC6', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'NESM3', 'NorESM2-LM', 'SAM0-UNICON', 'UKESM1-0-LL'};
% model value
lamda_cloud.modelVar = [0.51, 0.52, -0.36, 0.96, 1.17, 0.55, 0.54, 0.56, 0.80, 0.94, 0.29, 0.29, -0.01, 0.56, 0.00, -0.03, 0.79, -0.13, 0.38, -0.02, 0.12, 0.20, 0.38, 0.40, 0.36, 0.69, 0.81];
lamda_cloud.readme = 'lamda_cloud mean cloud climate feedback parameter, follow zelinka cal. for more info read /Research/p2_processCMIP6Data/s1.modelDataProcess/Info_lamda_cloud.m';

outPath='/data1/liuyincheng/cmip6-process/z_globalVar/';
save([outPath, 'lamda_cloud.mat'], 'lamda_cloud');
% usage example 
% lamda_cloudSpeci=lamda_cloud.modelVar(strcmp(lamda_cloud.modelName,'E3SM-1-0')==1);
