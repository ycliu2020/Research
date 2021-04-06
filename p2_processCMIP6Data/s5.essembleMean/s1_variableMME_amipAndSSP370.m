%%---------------------------------------------------------
% Author       : LYC
% Date         : 2021-04-02 15:48:14
% LastEditors  : Please set LastEditors
% Description  : Multi-model ensemble (MME) mean of CMIP6
% FilePath     : /code/p2_processCMIP6Data/s5.essembleMean/s1_variableMME_amipAndSSP370.m
%
%%---------------------------------------------------------

clc; clear; tic;
nowpath = pwd;
dbstop if error

% load mask map
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat') % load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat') % load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')

latRange = 90; % Latitude range
lon1 = [2.5 357.5]; lat1 = [-latRange + 1 latRange - 1]; % world area
toaSfc = {'toa', 'sfc'};
lon_k = 0:2.5:357.5; nlonk = length(lon_k); % kernel lat lon
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);

%% cal intersection of AMIP and SSP370
exmMdlAMIP2000 = onlySee_ripf1(1); % amip
exmMdlSSP370 = onlySee_ripf1(4); % ssp370

[MME_Models.name, exmMdlAMIP2000.index, exmMdlSSP370.index] = intersect(exmMdlAMIP2000.name', exmMdlSSP370.name');
MME_Models.readme = 'This essmble set contains the intersect of AMIP and SSP370, i.e., the total of 9 models';

varsOrig.outputPath = '/data1/liuyincheng/CMIP6-process/z_ensembleMean/MME1/';
% experiment
for exmNum = [1 2 4]
    % %% dradEffect_sfc_cld
    % varsOrig.Name = {'tsEffect', 'taEffect', 'albEffect', 'husEffect'};
    % varsOrig.fileFolder = 'radEffect/';
    % varsOrig.fileName = 'dradEffect_sfc_cld.mat';
    % % main auto calculation function
    % MME_autoCal3(varsOrig, MME_Models, exmNum);

    % %% dR_cloud
    % varsOrig.Name = {'dR_cloud_sfc', 'dR_cloud_toa'};
    % varsOrig.fileFolder = 'radEffect/';
    % varsOrig.fileName = 'dR_cloud.mat';
    % % main auto calculation function
    % MME_autoCal3(varsOrig, MME_Models, exmNum);

    % %% drlds
    % varsOrig.Name = {'drlds'};
    % varsOrig.fileFolder = 'anomaly/';
    % varsOrig.fileName = 'drlds.mat';
    % % main auto calculation function
    % MME_autoCal3(varsOrig, MME_Models, exmNum);

    % %% drsds
    % varsOrig.Name = {'drsds'};
    % varsOrig.fileFolder = 'anomaly/';
    % varsOrig.fileName = 'drsds.mat';
    % % main auto calculation function
    % MME_autoCal3(varsOrig, MME_Models, exmNum);

    % %% drsus
    % varsOrig.Name = {'drsus'};
    % varsOrig.fileFolder = 'anomaly/';
    % varsOrig.fileName = 'drsus.mat';
    % % main auto calculation function
    % MME_autoCal3(varsOrig, MME_Models, exmNum);

    %% trendRadEff_toa
    varsOrig.Name = {'trendyr_dRtoa_CRF', ...
                    'trendyr_dRtoa_ta', ...
                    'trendyr_dRtoa_tas1', ...
                    'trendyr_dRtoa_tas2', ...
                    'trendyr_dRtoa_nonTs', ...
                    'trendyr_dRtoa_residual', ...
                    'trendyr_dRtoa_cloud', ...
                    'trendyr_dRtoa_q', ...
                    'trendyr_dRtoa_alb', ...
                    'trendyr_dRtoa_ts'};
    varsOrig.fileFolder = 'radEffect_trend/';
    varsOrig.fileName = 'trend_dradEffect_toa_cld.mat';
    % main auto calculation function
    MME_autoCal2(varsOrig, MME_Models, exmNum);

    %% trendRadEff_sfc
    varsOrig.Name = {'trendyr_dRsfc_CRF', ...
                    'trendyr_dRsfc_ta', ...
                    'trendyr_dRsfc_taOnly', ...
                    'trendyr_dRsfc_tas1', ...
                    'trendyr_dRsfc_tas', ...
                    'trendyr_dRsfc_tas2', ...
                    'trendyr_dRatm_tsAtom', ...
                    'trendyr_dRsfc_nonTs', ...
                    'trendyr_dRsfc_residual', ...
                    'trendyr_dRsfc_cloud', ...
                    'trendyr_dRsfc_q', ...
                    'trendyr_dRsfc_alb', ...
                    'trendyr_dRsfc_ts'};
    varsOrig.fileFolder = 'radEffect_trend/';
    varsOrig.fileName = 'trend_dradEffect_sfc_cld.mat';
    % main auto calculation function
    MME_autoCal2(varsOrig, MME_Models, exmNum);

    %% trendanomaly_dts
    varsOrig.Name = {'trendyr_dts'};
    varsOrig.fileFolder = 'anomaly_trend/';
    varsOrig.fileName = 'trend_dts.mat';
    % main auto calculation function
    MME_autoCal2(varsOrig, MME_Models, exmNum);

    %% trendanomaly_drhs
    varsOrig.Name = {'trendyr_drhs'};
    varsOrig.fileFolder = 'anomaly_trend/';
    varsOrig.fileName = 'trend_drhs.mat';
    % main auto calculation function
    MME_autoCal2(varsOrig, MME_Models, exmNum);

    %% trendanomaly_drhsPlus
    varsOrig.Name = {'trendyr_drhsPlus'};
    varsOrig.fileFolder = 'anomaly_trend/';
    varsOrig.fileName = 'trend_drhsPlus.mat';
    % main auto calculation function
    MME_autoCal2(varsOrig, MME_Models, exmNum);

    %% trendanomaly_dhFlux
    varsOrig.Name = {'trendyr_dhFlux'};
    varsOrig.fileFolder = 'anomaly_trend/';
    varsOrig.fileName = 'trend_dhFlux.mat';
    % main auto calculation function
    MME_autoCal2(varsOrig, MME_Models, exmNum);

    %% trendanomaly_dnetTOA
    varsOrig.Name = {'trendyr_dnetTOA'};
    varsOrig.fileFolder = 'anomaly_trend/';
    varsOrig.fileName = 'trend_dnetTOA.mat';
    % main auto calculation function
    MME_autoCal2(varsOrig, MME_Models, exmNum);
end

eval(['cd ', nowpath]);
t = toc; disp(t)
