%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-10-14 15:3nValue:49
% LastEditTime : 2021-06-06 23:41:57
% LastEditors  : Please set LastEditors
% Description  :
% FilePath     : /code/p2_processCMIP6Data/s2.radEffTrend/contribAnalysis/rhs_Rts/MME_cntAls_joint_directMean.m
%
%%---------------------------------------------------------
clc; clear;
nowpath = pwd;

%% loop and process
[mlabels, areaStr] = cmipPlotParameters('sfc', 'land', 'radEffect'); % plot parameters

MMEType = 'MME1';
timeType = 'monthly'; % monthly or yearly
MMENum = str2num(MMEType(end));

exmStart = 1; exmEnd = 1;

for exmNum = exmStart:exmEnd

    if strcmp(MMEType, 'MME1')
        exmMdlAMIP2000 = onlySee_ripf1(1); % amip
        exmMdlSSP370 = onlySee_ripf1(4); % ssp370
        [MME_Models.name, ~, ~] = intersect(exmMdlAMIP2000.name', exmMdlSSP370.name');
    else
        [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
        level.readme_MME
        MME_Models.name = level.modlist_MME{MMENum, exmNum};
    end

    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    exmPath = fullfile('/data1/liuyincheng/CMIP6-process/', level.time1{exmNum});
    level.model2 = {};
    level.model2 = MME_Models.name;

    % model loop
    mdlStart = 1; mdlEnd = 1;%length(level.model2); % differnt models %length(level.model2)
    lineStart = 1; % 写入excel的起始位置
    mdlCount = 1; % 初始模式计数
    cov_MdlSet_dTs_x_sfc = zeros(3, 3, 4, length(level.model2)); % 方差矩阵, areaNum, modlNum

    for mdlNum = mdlStart:mdlEnd
        % model path
        mdlName = level.model2{mdlNum};
        mdlPath = fullfile(exmPath, level.model2{mdlNum});
        eval(['cd ', mdlPath]);
        disp(' ')
        disp([level.model2{mdlNum}, ' model start!'])

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

        for esmNum = specificNum:specificNum %1:1 %length(esmName) % note that r1i1p1 sometime not the first folder
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% load and read
            esmPath = fullfile(mdlPath, esmName{esmNum, 1});
            dradTrendPath = fullfile(esmPath, level.process3{7}); %/data1/liuyincheng/cmip6-process/amip_1980-2014/CESM2/Effect_trend
            dvarsPath = fullfile(esmPath, level.process3{2}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
            vsTsEffectPath = fullfile(esmPath, level.process3{9}); %/data1/liuyincheng/cmip6-process/amip_1980-2014/CESM2/vsTsEffect
            rawVarPath =  fullfile(esmPath, level.process3{1}); %/data1/liuyincheng/cmip6-process/amip_1980-2014/CESM2/vsTsEffect
            load([dradTrendPath, 'global_vars.mat']) % lat_f lon_f time plevf readme
            load([dvarsPath, 'global_vars.mat']) % lat_k lon_k time plevk readme
            % load([dvarsPath, 'dts.mat']) % dts clim_ts
            load([vsTsEffectPath, 'dTs_x_sfc.mat']) % dTs_alb, dTs_cloud, dTs_hus, dTs_residual, dTs_ta
            nlonf = length(lon_f); nlatf = length(lat_f); ntime = length(time.date);
            % dts = autoRegrid3(lon_k, lat_k, time.date, dts, lon_f, lat_f, time.date);

            % cal rhs(RHeating)
            load([dvarsPath, 'drlds.mat']) % surface_downwelling_longwave_flux_in_air
            load([dvarsPath, 'drlus.mat']) % surface_downwelling_longwave_flux_in_air
            load([dvarsPath, 'drsds.mat']) % surface_downwelling_shortwave_flux_in_air
            load([dvarsPath, 'drsus.mat']) % surface_upwelling_shortwave_flux_in_air
            drlds = autoRegrid3(lon_k, lat_k, time.date, drlds, lon_f, lat_f, time.date);
            drsds = autoRegrid3(lon_k, lat_k, time.date, drsds, lon_f, lat_f, time.date);
            drsus = autoRegrid3(lon_k, lat_k, time.date, drsus, lon_f, lat_f, time.date);
            drlus = autoRegrid3(lon_k, lat_k, time.date, drlus, lon_f, lat_f, time.date);

            dR_swnet = drsds - drsus; % sfc net shortwave flux
            drhting = drlds + dR_swnet; % equilibrium equation's RHS, nearly equal to sfc upward rad
            % load LH and SH
            load([dvarsPath, 'dhfls.mat']) % Surface Upward Latent Heat Flux
            load([dvarsPath, 'dhfss.mat']) % Surface Upward Sensible Heat Flux
            % unite define a vector component which is positive when directed downward
            dhfls = autoRegrid3(lon_k, lat_k, time.date, dhfls, lon_f, lat_f, time.date);
            dhfss = autoRegrid3(lon_k, lat_k, time.date, dhfss, lon_f, lat_f, time.date);
            dhFlux = -dhfls - dhfss; % LH+SH

            % read raw hfls hfss
            load([rawVarPath, 'hfls.mat']) % Surface Upward Latent Heat Flux
            load([rawVarPath, 'hfss.mat']) % Surface Upward Sensible Heat Flux

            load([rawVarPath, 'rlds.mat']) % surface_downwelling_longwave_flux_in_air
            load([rawVarPath, 'rlus.mat']) % surface_downwelling_longwave_flux_in_air
            load([rawVarPath, 'rsds.mat']) % surface_downwelling_shortwave_flux_in_air
            load([rawVarPath, 'rsus.mat']) % surface_upwelling_shortwave_flux_in_air


            %% test
            % read sfc values
            dradEffectPath = fullfile(esmPath, level.process3{6}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
            load([dradEffectPath, 'dR_cloud.mat']) % dR_cloud_sfc
            load([dradEffectPath, 'dradEffect_sfc_cld.mat']) % albEffect, husEffect, nonCloudAndTsEffect, taEffect, taOnlyEffect, tasEffect, tasEffect1, tasEffect2, totalEffect, tsEffect, wvlwEffect, wvswEffect
            load([dradEffectPath, 'dR_residual_cld_sfc.mat']) % dR_resiual_cld_sfc
            % load([dradEffectPath, 'model_dradEffect.mat']) % 'l_rad', 's_rad', 'dR_allsky', 'dR_clr', 'readme_realradEfect'
            drhs=drlus;%-tsEffect;
            %% Test1: Rheating variance of a sum test
            dTs_x_sfc = zeros(nlonf, nlatf, ntime, 2);
            dTs_x_sfc(:, :, :, 1) = drhs;
            dTs_x_sfc(:, :, :, 2) = drhs-drhting+dhFlux;
            dTs_x_sfc(:, :, :, 3) = drhting;
            dTs_x_sfc(:, :, :, 4) = -dhFlux;

            a1=dTs_x_sfc(:, :, 1, 1);
            a2=dTs_x_sfc(:, :, 1, 2);
            a3=dTs_x_sfc(:, :, 1, 3);
            a4=dTs_x_sfc(:, :, 1, 4);

            ac_hfls=clim_hfls(:,:,1);
            ac_hfss=clim_hfss(:,:,1);
            
            ah_hfls=dhfls(:,:,1);
            ah_hfss=dhfss(:,:,1);

            aRaw_hfls=hfls(:,:,1);
            aRaw_hfss=hfss(:,:,1);

            RawBalance= rsds - rsus+rlds-rlus+hfls+hfss;
            dvarBalance= drsds - drsus+drlds-drlus+dhfls+dhfss;
            aD1=dvarBalance(:,:,1);
            aR1=RawBalance(:,:,1);
            % dTs_x_sfc(:, :, :, 1) = drhs;
            % dTs_x_sfc(:, :, :, 2) = drhting - (nonCloudAndTsEffect + dR_cloud_sfc)-dhFlux; % or -squeeze(dR_allsky(:,:,:,1))+dR_residual_cld_sfc;
            % dTs_x_sfc(:, :, :, 3) = dhFlux; % or -squeeze(dR_allsky(:,:,:,1))+dR_residual_cld_sfc;
            % dTs_x_sfc(:, :, :, 4) = albEffect;
            % dTs_x_sfc(:, :, :, 5) = dR_cloud_sfc;
            % dTs_x_sfc(:, :, :, 6) = husEffect;
            % dTs_x_sfc(:, :, :, 7) = taEffect;
            % dTs_x_sfc(:,:,:,3)=dR_residual_cld_sfc;
            %
            % cal Function
            latRange = 90;
            areaStr = {'world', 'china east', 'USA east', 'EUR west'};
            
            for areaNum = 1:length(areaStr)
                [cov_glbMoth_dTs_x_sfc] = calCovContribution(latRange, lat_f, timeType, ntime, areaStr, areaNum, dTs_x_sfc);
                cov_MdlSet_dTs_x_sfc(:, :, areaNum, mdlCount) = cov_glbMoth_dTs_x_sfc;
            end

            disp([esmName{esmNum, 1}, ' ensemble is done!'])
            mdlCount = mdlCount + 1;

        end

        disp([level.model2{mdlNum}, ' model is done!'])
        disp(' ')
    end
    % cal average
    cov_MdlSetMean_dTs_x_sfc=mean(cov_MdlSet_dTs_x_sfc, 4);
    % save to excel
    outPutFile = ['/home/liuyc/Research/P02.Ts_change_research/table/Radiation_Tscontribution/',MMEType,'_DirectMean_radContrib_',  level.time1{exmNum}(6:end - 1), '_', timeType, '_Rts.xlsx'];
    for areaNum = 1:length(areaStr)
        [lineStart] = saveToExcl(MMEType, areaNum, areaStr, dTs_x_sfc, outPutFile, lineStart, cov_MdlSetMean_dTs_x_sfc(:,:,areaNum));
    end

    disp([level.time1{exmNum}, ' era is done!'])
    disp(' ')
end

eval(['cd ', nowpath]);

function [cov_glbMoth_dTs_x_sfc] = calCovContribution(latRange, lat_f, timeType, ntime, areaStr, areaNum, dTs_x_sfc)

    size_dTs_x_sfc = size(dTs_x_sfc);
    nValue = size_dTs_x_sfc(4);
    % Mask Part
    for x_ind = 1:nValue
        [dTs_x_sfc(:, :, :, x_ind), ~, ~] = maskArea(dTs_x_sfc(:, :, :, x_ind), lat_f, latRange, -latRange, areaStr{areaNum});
    end

    % cal the weight average value

    glb_dTs_x_sfc = zeros(ntime, nValue);

    for timeNum = 1:ntime

        for x_ind = 1:nValue
            glb_dTs_x_sfc(timeNum, x_ind) = areaMeanLatWeight(dTs_x_sfc(:, :, timeNum, x_ind), lat_f);
        end

    end

    glbMoth_dTs_x_sfc = glb_dTs_x_sfc;

    if strcmp(timeType, 'yearly')
        % cal intel annual mean
        ntime_year = ntime / 12; % year num
        glbMoth_dTs_x_sfc = zeros(ntime_year, nValue);
        countNum = 1;

        for timeNum = 1:ntime_year
            glbMoth_dTs_x_sfc(timeNum, :) = sum(glb_dTs_x_sfc(countNum:countNum + 12 - 1, :), 1);
            countNum = countNum + 12;
        end

    end

    % test sum = 0
    % for timeNum = 1:ntime
    %     glb_dTs_x_sfc_test(timeNum)=sum(glb_dTs_x_sfc(timeNum,:));

    % end

    %% cal variance and covariance
    cov_glbMoth_dTs_x_sfc = cov(glbMoth_dTs_x_sfc(:, 2:end));
    var_glbMoth_dTs_sfc = var(glbMoth_dTs_x_sfc(:, 1));
    cov_glbMoth_dTs_x_sfc = cov_glbMoth_dTs_x_sfc ./ var_glbMoth_dTs_sfc;
    % transfor into paper used(lower triangular matrix)
    cov_triu1 = triu(cov_glbMoth_dTs_x_sfc, 1) .* 2;
    cov_diag = diag(diag(cov_glbMoth_dTs_x_sfc));
    cov_glbMoth_dTs_x_sfc = cov_triu1 + cov_diag;
    cov_glbMoth_dTs_x_sfc = fliplr(cov_glbMoth_dTs_x_sfc);
    disp([areaStr{areaNum}, ' Rheating decompose: '])
    disp(cov_glbMoth_dTs_x_sfc)
    % check sum equal to 1
    check_sum = sum(sum(cov_glbMoth_dTs_x_sfc));
    disp('check sum: ')
    disp(check_sum)
end

function [lineStart] = saveToExcl(MMEType, areaNum, areaStr, dTs_x_sfc, outPutFile, lineStart, cov_glbMoth_dTs_x_sfc)
    disp([areaStr{areaNum}, ' Rheating decompose: '])
    disp(cov_glbMoth_dTs_x_sfc)
    % check sum equal to 1
    check_sum = sum(sum(cov_glbMoth_dTs_x_sfc));
    disp('check sum: ')
    disp(check_sum)
    
    size_dTs_x_sfc = size(dTs_x_sfc);
    nValue = size_dTs_x_sfc(4);

    % save as to excel

    % Model and regional name

    % value name
    if areaNum == 1
        outPutTxt = {MMEType; 'E1'; 'SH+LH'; 'Alb'; 'Cld'; 'WV'; 'Ta'};
    else
        outPutTxt = {''; 'E1'; 'SH+LH'; 'Alb'; 'Cld'; 'WV'; 'Ta'};
    end

    outPutTxt_rev = fliplr(outPutTxt');
    nValue = nValue - 1; % 减去方程左式的一项
    writecell(outPutTxt, outPutFile, 'Sheet', 1, 'Range', ['A', num2str(lineStart), ':A', num2str(lineStart + nValue)])
    writecell(outPutTxt_rev(1:end - 1), outPutFile, 'Sheet', 1, 'Range', ['B', num2str(lineStart + nValue + 1), ':G', num2str(lineStart + nValue + 1)])
    writematrix(cov_glbMoth_dTs_x_sfc, outPutFile, 'Sheet', 1, 'Range', ['B', num2str(lineStart + 1), ':G', num2str(lineStart + nValue)])

    regionTxt = areaStr{areaNum};
    writematrix(regionTxt, outPutFile, 'Sheet', 1, 'Range', ['B', num2str(lineStart), ':B', num2str(lineStart)])

    lineStart = lineStart + 8;

end

% if exmNum == 1
%     % cut to 2000.03-2014.02
%     timeStr = string(datestr(datenum(time.date), 'yyyy-mm'));
%     cutStart = find(timeStr == '2000-03');
%     cutEnd = find(timeStr == '2014-02');
%     time = time.date(cutStart:cutEnd);
%     ntime = length(time);
%     dTs_x_sfc = dTs_x_sfc(:, :, cutStart:cutEnd, :);
% end
