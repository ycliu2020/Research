%%---------------------------------------------------------
% Author       : LYC
% Date         : 2021-04-02 15:48:14
% LastEditors  : Please set LastEditors
% Description  : Multi-model ensemble (MME) mean of CMIP6
% FilePath     : /code/p2_processCMIP6Data/s5.essembleMean/s1_varMME_joint.m
%
%%---------------------------------------------------------

clc; clear; tic;
nowpath = pwd;
dbstop if error
run '/home/liuyc/lib/tools/matlab/myTools/autoScript/preLoadVar.m'
MMEType='MME3';
MMENum=str2num(MMEType(end));
% experiment
for exmNum = [1 2 4]
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    level.readme_MME
    MME_Models.name = level.modlist_MME{MMENum, exmNum};
    varsOrig.outputPath = fullfile(level.path_MME, MMEType);%/data1/liuyincheng/CMIP6-process/2000-2014/
    MME_Models.readme=['This essmble set contains the selection of ',Experiment{exmNum}];
    
    %% dradEffect_sfc_cld
    varsOrig.Name = {'tsEffect', 'taEffect', 'albEffect', 'husEffect'};
    varsOrig.fileFolder = 'radEffect/';
    varsOrig.fileName = 'dradEffect_sfc_cld.mat';
    % main auto calculation function
    MME_autoCal3(varsOrig, MME_Models, exmNum);

    %% dR_cloud
    varsOrig.Name = {'dR_cloud_sfc', 'dR_cloud_toa'};
    varsOrig.fileFolder = 'radEffect/';
    varsOrig.fileName = 'dR_cloud.mat';
    % main auto calculation function
    MME_autoCal3(varsOrig, MME_Models, exmNum);

    %% drlds
    varsOrig.Name = {'drlds'};
    varsOrig.fileFolder = 'anomaly/';
    varsOrig.fileName = 'drlds.mat';
    % main auto calculation function
    MME_autoCal3(varsOrig, MME_Models, exmNum);

    %% drsds
    varsOrig.Name = {'drsds'};
    varsOrig.fileFolder = 'anomaly/';
    varsOrig.fileName = 'drsds.mat';
    % main auto calculation function
    MME_autoCal3(varsOrig, MME_Models, exmNum);

    %% drsus
    varsOrig.Name = {'drsus'};
    varsOrig.fileFolder = 'anomaly/';
    varsOrig.fileName = 'drsus.mat';
    % main auto calculation function
    MME_autoCal3(varsOrig, MME_Models, exmNum);

    %% drlus
    varsOrig.Name = {'drlus'};
    varsOrig.fileFolder = 'anomaly/';
    varsOrig.fileName = 'drlus.mat';
    % main auto calculation function
    MME_autoCal3(varsOrig, MME_Models, exmNum);

    %% dhfss
    varsOrig.Name = {'dhfss'};
    varsOrig.fileFolder = 'anomaly/';
    varsOrig.fileName = 'dhfss.mat';
    % main auto calculation function
    MME_autoCal3(varsOrig, MME_Models, exmNum);

    %% dhfls
    varsOrig.Name = {'dhfls'};
    varsOrig.fileFolder = 'anomaly/';
    varsOrig.fileName = 'dhfls.mat';
    % main auto calculation function
    MME_autoCal3(varsOrig, MME_Models, exmNum);
    
    % trendRadEff_toa
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

    %% dR_residual_cld_sfc
    varsOrig.Name = {'dR_residual_cld_sfc', 'lw_residual_cld_sfc', 'sw_residual_cld_sfc'};
    varsOrig.fileFolder = 'radEffect/';
    varsOrig.fileName = 'dR_residual_cld_sfc.mat';
    % main auto calculation function
    MME_autoCal3(varsOrig, MME_Models, exmNum);

    %% dR_residual_cld_toa
    varsOrig.Name = {'dR_residual_cld_toa', 'lw_residual_cld_toa', 'sw_residual_cld_toa'};
    varsOrig.fileFolder = 'radEffect/';
    varsOrig.fileName = 'dR_residual_cld_toa.mat';
    % main auto calculation function
    MME_autoCal3(varsOrig, MME_Models, exmNum);
end

eval(['cd ', nowpath]);
t = toc; disp(t)
