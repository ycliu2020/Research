%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-08-31 17:00:15
% LastEditTime : 2020-09-14 17:15:01
% LastEditors  : LYC
% Description  : 计算不同变量间的时间序列相关性并画图
% FilePath     : /code/p2_processCMIP6Data/s2.radEffTrend/timeSeriAnalysis/timeSeriAns_globalDistribution.m
%
%%---------------------------------------------------------
clc; clear; tic;
nowpath = pwd;
% load mask map
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')

latRange = 60; % Latitude range
lon1 = [2.5 357.5]; lat1 = [-latRange + 1 latRange - 1]; % world area
toaSfc = {'toa', 'sfc'};
lon_k = 0:2.5:357.5; nlonk = length(lon_k); % kernel lat lon
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment
for exmNum = 1:4%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    %CAMS-CSM1-0 didn't have sfc clear sky radiation, delete it
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % exmPath
    exmPath = ['/data1/liuyincheng/cmip6-process/', level.time1{exmNum}]; %/data1/liuyincheng/cmip6-process/2000-2014/
    mPath.uniOutput = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/02.cmip6Result/timeSeriesAns/timeCC_globalDistribution', level.time1{exmNum});%['dRTs_',lower(mlabels.level)],
    mPath.Output = fullfile(mPath.uniOutput);
    auto_mkdir(mPath.Output)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model
    for mdlNum = 1:length(level.model2)
        % model path
        mdlPath = fullfile(exmPath, level.model2{mdlNum});
        eval(['cd ', mdlPath]);
        disp(' ')
        disp([level.model2{mdlNum}, ' model start!'])

        % ensemble member path
        esmName = getPath_fileName(mdlPath, '.');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensemble member
        for esmNum = 1:length(esmName)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% load and read
            esmPath = fullfile(mdlPath, esmName{esmNum, 1});
            % data path
            varsPath = fullfile(esmPath, level.process3{1}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/rawdata_regrid
            dvarsPath = fullfile(esmPath, level.process3{2}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
            dvarsTrendPath = fullfile(esmPath, level.process3{3}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly_trend
            kernelPath = fullfile(esmPath, level.process3{5}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/kernelsCal
            dradEffectPath = fullfile(esmPath, level.process3{6}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
            dnonLocalCldPath = fullfile(esmPath, level.process3{8}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/non_localCld/
            load([dradEffectPath, 'global_vars.mat'])% lat_f lon_f time plev_k readme
            ntime = length(time.date);

            % cal surface temperature
            load([varsPath, 'ts.mat'])% ts
            load([dvarsPath, 'dts.mat'])%dts
            dts = autoRegrid3(lon_k, lat_k, time.date, dts, lon_f, lat_f, time.date);
            ts = autoRegrid3(lon_k, lat_k, time.date, ts, lon_f, lat_f, time.date);

            % cal rhs(RHeating)
            load([dvarsPath, 'drlds.mat'])% surface_downwelling_longwave_flux_in_air
            load([dvarsPath, 'drsds.mat'])% surface_downwelling_shortwave_flux_in_air
            load([dvarsPath, 'drsus.mat'])% surface_upwelling_shortwave_flux_in_air
            dR_swnet = drsds - drsus; % sfc net shortwave flux
            drhs = drlds + dR_swnet; % equilibrium equation's RHS, nearly equal to sfc upward rad
            drhs = autoRegrid3(lon_k, lat_k, time.date, drhs, lon_f, lat_f, time.date);

            % cloud radRffect
            load([dradEffectPath, 'dR_cloud_toa.mat'])
            load([dradEffectPath, 'dR_cloud_sfc.mat'])

            % TotalEffect
            load([dradEffectPath, 'dradEfect_toa_cld.mat'])% 'totalEffect', 'wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect', 'taEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect', 'mainEffect'
            load([dradEffectPath, 'real_dradEfect.mat'])% 'dR_allsky', 'l_rad', 's_rad', 'dR_clr', 'readme_realradEfect'
            load([dradEffectPath, 'dR_residual_cld_sfc.mat'])% dR_residual_cld_sfc
            load([dradEffectPath, 'dR_residual_cld_toa.mat'])% dR_residual_cld_toa
            dR_netSfc = squeeze(dR_allsky(:, :, :, 1));
            dR_netTOA = squeeze(dR_allsky(:, :, :, 2));
            varUsed = zeros(nlonf, nlatf, ntime, 2);
            varUsed(:, :, :, 1) = dts;
            varUsed(:, :, :, 2) = drhs;
            varNames = {'dTs', 'dRHeating', 'ta RadEfect', 'ts RadEfect', 'q RadEfect', 'alb RadEfect'};
            yLabel = {'K', 'W/m2', 'W/m2', 'W/m2', 'W/m2', 'W/m2'};

            % varUsed(:, :, :, 1) = dR_netSfc;
            % varUsed(:, :, :, 2) = dR_cloud_sfc;
            % varUsed(:, :, :, 3) = taEffect;
            % varUsed(:, :, :, 4) = tsEffect;
            % varUsed(:, :, :, 5) = wvlwEffect+wvswEffect;
            % varUsed(:, :, :, 6) = albEffect;
            % varNames = {'Total radiation', 'cloud RadEfect', 'ta RadEfect', 'ts RadEfect', 'q RadEfect', 'alb RadEfect'};
            % transfor to yearly Data
            sizeVarUsed = size(varUsed);
            sizeLon=sizeVarUsed(1);sizeLat=sizeVarUsed(2);sizeTime=sizeVarUsed(3);sizeVar=sizeVarUsed(4);
            varUsedYearly = zeros(sizeLon, sizeLat, sizeTime / 12, sizeVar);

            for varNum = 1:sizeVar
                varUsedYearly(:, :, :, varNum) = monthlyToYearly(varUsed(:, :, :, varNum));
            end
            % mask only -60-60 land.
            [varUsedYearly, ~, ~] = maskArea(varUsedYearly,lat_f,latRange,-latRange, 'world');
            
            % cal areaMeanLatWeight
            sizeVarUsedYearly=size(varUsedYearly);
            varUsedYearly_weightMean=zeros(sizeVarUsedYearly(3),sizeVar);
            for varNum = 1:sizeVar
                for timeNum = 1 : sizeVarUsedYearly(3)
                    varUsedYearly_weightMean(timeNum,varNum)=areaMeanLatWeight(varUsedYearly(:,:,timeNum,varNum), lat_f);
                end
            end
            time.vec = datevec(time.date);
            timeYearly = unique(time.vec(:, 1));
            % cal cc of areaMeanLatWeight
            sizeWeightMean = size(varUsedYearly_weightMean);
            cc_weightMean=cell(1,sizeWeightMean(2));
            pp_weightMean=cell(1,sizeWeightMean(2));
            for varNum = 1 : sizeWeightMean(2)
                [cc0,pp0] = corrcoef(varUsedYearly_weightMean(:,1),varUsedYearly_weightMean(:,varNum),'Rows','complete');
                cc_weightMean{varNum} = roundn(cc0(1,2),-2);%保留俩位小数
                pp_weightMean{varNum} = pp0(1,2);% confidence interval
            end
        
            % cal cc of every point of map
            cc_global=zeros(nlonf,nlatf,sizeVar);
            pp_global=zeros(nlonf,nlatf,sizeVar);
            for varNum = 1 : sizeVar
                for latNum = 1:nlatf
                    for lonNum =  1:nlonf
                        varUsedTemp=squeeze(squeeze(varUsedYearly(lonNum, latNum, :, :)));
                        [cc0,pp0] = corrcoef(varUsedYearly(lonNum, latNum, :, 1), varUsedYearly(lonNum, latNum, :, varNum), 'Rows', 'complete');
                        cc_global(lonNum,latNum,varNum)= roundn(cc0(1,2),-2);%保留俩位小数
                        pp_global(lonNum,latNum,varNum)= pp0(1,2);% confidence interval
                    end
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % plot cc global
            f_row = 1; f_col = 1; % 设置画图的行列
            set(0, 'DefaultFigureVisible', 'on')
            ss = get(0, 'ScreenSize'); 
            coef_amplify = 3;
            h = figure('Position', [ss(4)/2-100 ss(3) / 35 ss(3)/5*f_col*coef_amplify (ss(4)-80)/5*f_row*coef_amplify]);
            % clf reset;
            set(h, 'Color', [1 1 1]);
            f_matrix = reshape(1:f_row*f_col, [f_col, f_row])';

            % figure
            for varNum = 1:f_row * f_col
                trendz = squeeze(cc_global(:, :, 2));
                [plotRow, plotCol] = find(f_matrix == varNum);
                subplot_yc(f_row, f_col, plotRow, plotCol); 
                hold on
                m_proj('Mercator', 'lon_f', lon1, 'lat_f', lat1); %Mercator,Equidistant cylindrical,lambert,Miller Cylindrical
                m_pcolor(lon_f, lat_f, trendz');
                colormap(mycolor(18)); %mycolor(100)is soden color????????colormap(flipud(mycolor(13)));%colormap(jet(4))
                col_SeriesNum=10;
                % [colorbar_Series] = findSuit_colorInt(trendz, col_SeriesNum);
                max_color=1;
                min_color=-max_color;
                caxis([min_color max_color]);
                hold on
                m_line(world_mapx(:), world_mapy(:), 'color', [0 0 0], 'LineWidth', 0.5);
                m_grid('linestyle', 'none', 'tickdir', 'out', 'yaxislocation', 'left', 'fontsize', 8, 'color', 'k');

                headLineTxt = {'The correlation coefficient of dTs and dRHeating',['Level:', toaSfc{2}, ', Era: ', level.time1{exmNum}(1:end - 1),', Model:', level.model2{mdlNum}, ', Ensemble: ', esmName{esmNum}], ['global mean CC= ',num2str(cc_weightMean{2})]};
                title(headLineTxt,'Interpreter','none','fontsize', 14); % cc=',num2str(corr))
                hold on
                c = colorbar;
                % c.TickLength = 0.0245;
                % c.Limits = [min_color min_color];
            end
    

            % save figures
            figName = [level.time1{exmNum}(1:end - 1), '_', level.model2{mdlNum},'_',esmName{esmNum}];
            figurePath = [mPath.Output, '/', figName, '.png'];
            saveas(gcf, figurePath)
            % save_png(figurePath)%high resolution
            close gcf
    
        end

    end

end
eval(['cd ', nowpath]);
t = toc; disp(t)


