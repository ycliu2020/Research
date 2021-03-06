%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-10-14 15:3nValue:49
% LastEditTime : 2020-11-02 16:24:58
% LastEditors  : LYC
% Description  :
% FilePath     : /code/p2_processCMIP6Data/s2.radEffTrend/contribAnalysis/contribAnalysis_eastUSA.m
%
%%---------------------------------------------------------
clc; clear;
nowpath = pwd;

%% calculate mask file
% load mask map
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/中国干湿以及青藏高原分区地图/main_use/mask720.mat')
% 国境线
bou_china = shaperead('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/map_data/中国及各省shp/国界线/bou1_4p.shp'); % load china  boundary
bou_chinaX = [bou_china(:).X]; bou_chinaY = [bou_china(:).Y];
% 国境线(带南海线)
bou_china_line = shaperead('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/map_data/中国及各省shp/国界线/bou1_4l.shp'); % load china  boundary
bou_china_lineX = [bou_china_line(:).X]; bou_china_lineY = [bou_china_line(:).Y];
% 省界线
bou_chinaProvince = shaperead('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/map_data/中国及各省shp/省界线/bou2_4p.shp'); % load china  boundary
bou_chinaProvinceX = [bou_chinaProvince(:).X]; bou_chinaProvinceY = [bou_chinaProvince(:).Y];

%save maskFile to area of china east
maskUSA_east = maskworld_cp;
maskUSA_east(lonw < -90 | lonw>-72, :) = 0;
maskUSA_east(:, latw < 30 | latw>45) = 0;
readme_mask_USA_east = 'area of east USA, Range: 288>lonw>270E,30<latw<45N; Based on maskworld_cp(lon, lat same as it)';
save('/home/liuyc/lib/tools/matlab/plot/myMap/03.other_area/mask14472_eastUSA.mat', 'maskUSA_east', 'lat14472', 'lon14472', 'readme_mask_USA_east')

%% loop and process
[mlabels, areaNum] = cmipPlotParameters('atm', 'land', 'radEffect'); % plot parameters

esm = 'r1i1p1f1';
exmStart = 1; exmEnd = 1;

for exmNum = exmStart:exmEnd
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % mPath.input:E:/data/cmip6-process/2000-2014/
    exmPath = fullfile('/data1/liuyincheng/cmip6-process/', level.time1{exmNum});
    % mPath.output:a_research/P02.Ts_change_research/figure/04.cmip6Result/2000-2014/ 根据输出的格式进行修改
    % mPath.uniOutput = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/02.cmip6Result/1.6/', lower(mlabels.level), level.time1{exmNum});
    % mPath.Output = fullfile(mPath.uniOutput);
    % auto_mkdir(mPath.Output)

    % model loop
    mdlStart = 1; mdlEnd = length(level.model2); % differnt models%length(level.model2)

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
        for esmNum = specificNum:specificNum%1:1%length(esmName)% note that r1i1p1 sometime not the first folder
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% load and read
            esmPath = fullfile(mdlPath, esmName{esmNum, 1});
            dradTrendPath = fullfile(esmPath, level.process3{7}); %/data1/liuyincheng/cmip6-process/amip_1980-2014/CESM2/Effect_trend
            dvarsPath = fullfile(esmPath, level.process3{2}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
            vsTsEffectPath = fullfile(esmPath, level.process3{9}); %/data1/liuyincheng/cmip6-process/amip_1980-2014/CESM2/vsTsEffect
            load([dradTrendPath, 'global_vars.mat'])% lat_f lon_f time plevf readme
            load([dvarsPath, 'global_vars.mat'])% lat_k lon_k time plevk readme
            load([dvarsPath, 'dts.mat'])% dts clim_ts
            load([vsTsEffectPath, 'dTs_x_sfc.mat'])% dTs_alb, dTs_cloud, dTs_hus, dTs_residual, dTs_ta
            nlonf = length(lon_f); nlatf = length(lat_f); ntime = length(time.date);
            dts = autoRegrid3(lon_k, lat_k, time.date, dts, lon_f, lat_f, time.date);

            % cal rhs(RHeating)
            load([dvarsPath, 'drlds.mat'])% surface_downwelling_longwave_flux_in_air
            load([dvarsPath, 'drsds.mat'])% surface_downwelling_shortwave_flux_in_air
            load([dvarsPath, 'drsus.mat'])% surface_upwelling_shortwave_flux_in_air
            dR_swnet = drsds - drsus; % sfc net shortwave flux
            drhs = drlds + dR_swnet; % equilibrium equation's RHS, nearly equal to sfc upward rad
            drhs = autoRegrid3(lon_k, lat_k, time.date, drhs, lon_f, lat_f, time.date);
            
            %%test
            dradEffectPath = fullfile(esmPath, level.process3{6}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
            load([dradEffectPath, 'dR_cloud_sfc.mat'])% dR_cloud_sfc
            load([dradEffectPath, 'dradEfect_sfc_cld.mat'])% albEffect, husEffect, mainEffect, taEffect, taOnlyEffect, taOnlyEffect2, tasEffect, tasEffect2, totalEffect, tsEffect, wvlwEffect, wvswEffect
            load([dradEffectPath, 'dR_residual_cld_sfc.mat'])% dR_resiual_cld_sfc
            load([dradEffectPath, 'real_dradEfect.mat'])% 'l_rad', 's_rad', 'dR_allsky', 'dR_clr', 'readme_realradEfect'

            % add all vars into one var
            % dTs_x_sfc=zeros(nlonf, nlatf, ntime,6);
            % dTs_x_sfc(:,:,:,1)=dts;
            % dTs_x_sfc(:,:,:,2)=dTs_cloud;
            % dTs_x_sfc(:,:,:,3)=dTs_ta;
            % dTs_x_sfc(:,:,:,4)=dTs_hus;
            % dTs_x_sfc(:,:,:,5)=dTs_alb;
            % dTs_x_sfc(:,:,:,6)=dTs_residual;

            dTs_x_sfc = zeros(nlonf, nlatf, ntime, 2);
            dTs_x_sfc(:, :, :, 1) = -tsEffect;
            dTs_x_sfc(:, :, :, 2) = -tsEffect - mainEffect - dR_cloud_sfc; %-squeeze(dR_allsky(:,:,:,1))+dR_residual_cld_sfc;
            dTs_x_sfc(:, :, :, 3) = albEffect;
            dTs_x_sfc(:, :, :, 4) = dR_cloud_sfc;
            dTs_x_sfc(:, :, :, 5) = husEffect;
            dTs_x_sfc(:, :, :, 6) = taEffect;
            % dTs_x_sfc(:,:,:,3)=dR_residual_cld_sfc;
            size_dTs_x_sfc = size(dTs_x_sfc);
            nValue = size_dTs_x_sfc(4);

            % mask east china
            latRange = 80;
            areaNum = 'USA east';

            for x_ind = 1:nValue
                [dTs_x_sfc(:, :, :, x_ind), ~, ~] = maskArea(dTs_x_sfc(:, :, :, x_ind), lat_f, latRange, -latRange, areaNum);
            end

            % cal the weight average value of east china

            glb_dTs_x_sfc = zeros(ntime, nValue);

            for timeNum = 1:ntime

                for x_ind = 1:nValue
                    glb_dTs_x_sfc(timeNum, x_ind) = areaMeanLatWeight(dTs_x_sfc(:, :, timeNum, x_ind), lat_f);
                end

            end

            % cal intel annual mean
            ntime_year = ntime / 12; % year num
            glbMoth_dTs_x_sfc = zeros(ntime_year, nValue);
            countNum = 1;

            for timeNum = 1:ntime_year
                glbMoth_dTs_x_sfc(timeNum, :) = sum(glb_dTs_x_sfc(countNum:countNum + 12 - 1, :), 1);
                countNum = countNum + 12;
            end

            % glbMoth_dTs_x_sfc=glb_dTs_x_sfc;
            % % test sum = 0
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
            disp('Rheating decompose: ')
            disp(cov_glbMoth_dTs_x_sfc)
            % check sum equal to 1
            check_sum=sum(sum(cov_glbMoth_dTs_x_sfc));
            disp('check sum: ')
            disp(check_sum)

            % save as to excel
            % add var annotation
            % outPutTxt={level.model2{mdlNum};'Res';'Alb';'Cld';'WV';'Ta'};
            % outPutTxt_rev=fliplr(outPutTxt');
            % outPutFile=['/home/liuyc/Research/P02.Ts_change_research/figure/proj2_cmip6Result/Radiation_Tscontribution/radContrib_',level.time1{exmNum}(1:end-1),'.xlsx'];
            % lineStart=8*(mdlNum-1)+1;
            % writecell(outPutTxt,outPutFile,'Sheet',1,'Range',['A',num2str(lineStart),':A',num2str(lineStart+6)])
            % writecell(outPutTxt_rev(1:end-1),outPutFile,'Sheet',1,'Range',['B',num2str(lineStart+nValue),':G',num2str(lineStart+nValue)])
            % writematrix(cov_glbMoth_dTs_x_sfc,outPutFile,'Sheet',1,'Range',['B',num2str(lineStart+1),':G',num2str(lineStart+6)])

            disp([esmName{esmNum, 1}, ' ensemble is done!'])

        end

        disp([level.model2{mdlNum}, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{exmNum}, ' era is done!'])
    disp(' ')
end

eval(['cd ', nowpath]);
