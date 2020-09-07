%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-17 15:12:08
% LastEditTime : 2020-09-06 15:39:32
% LastEditors  : LYC
% Description  : only cal ensember r1i1p1f1
% FilePath     : /code/p2_processCMIP6Data/s3.nonLocalCld/plot/s_plot3sub_nonLocalCld3.m
%  
%%---------------------------------------------------------

clc; clear;
nowpath = pwd;
% load mask map
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')

[mlabels, areaNum] = cmipPlotParameters('sfc', 'land', 'nonLocalCld3_speci'); % plot parameters mlabels.level(SFC/TOA)
esm = 'r1i1p1f1';
% Latitude range
latRange = 60;
lon1 = [2.5 357.5]; lat1 = [-latRange + 1 latRange - 1]; % world area
set(0, 'defaultfigurecolor', 'w')

exm_left = 4; exm_right = 4;
for exmName = exm_left:exm_right
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmName);
    % mPath.input:E:/data/cmip6-process/2000-2014/

    mPath.input = fullfile('/data1/liuyincheng/cmip6-process/', level.time1{exmName});
    % mPath.output:a_research/P02.Ts_change_research/figure/04.cmip6Result/2000-2014/
    mPath.uniOutput = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/02.cmip6Result/nonLocalCld/nonLocalCld3/3subPlot_Res', level.time1{exmName});%['dRTs_',lower(mlabels.level)],
    mPath.Output = fullfile(mPath.uniOutput);
    auto_mkdir(mPath.Output)

    % model loop
    mdl_left = 1; mdl_right = length(level.model2); % differnt models%length(level.model2)
    for mdlName = mdl_left:mdl_right
        mdlPath=fullfile(mPath.input, level.model2{mdlName});
        % load data
        varsPath = fullfile(mdlPath, esm, level.process3{1}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/rawdata
        dvarsPath = fullfile(mdlPath, esm, level.process3{2}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
        dvarsTrendPath = fullfile(mdlPath, esm, level.process3{3}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly_trend
        kernelPath = fullfile(mdlPath, esm, level.process3{5}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/kernelsCal
        dradEffectPath = fullfile(mdlPath, esm, level.process3{6}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
        dradTrendPath = fullfile(mdlPath, esm, level.process3{7}); %/data1/liuyincheng/cmip6-process/amip_1980-2014/CESM2/radEffect_trend
        dnonLocalCldPath = fullfile(mdlPath, esm, level.process3{8}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/non_localCld/
        vsTsEffectTrendPath = fullfile(mdlPath, esm, level.process3{10}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/vsTsEffect_trend/
        if ~exist(dradTrendPath,'dir')
            disp(['the ', esm, ' ensemble of ',level.model2{mdlName}, ' didnt exist']);
            continue
        end

        load([dvarsTrendPath, 'global_vars.mat'])% lat_f lon_f time plevk readme
        load([dradTrendPath, 'trend_dradEfect_toa_cld.mat'])% 10 vars:'trendyr_dRtoa_ta','trendyr_dRtoa_taOnly2', 'trendyr_dRtoa_tas2., 'trendyr_dRtoa_tsAtom', 'trendyr_dRtoa_mainEffect', 'trendyr_dRtoa_residual', 'trendyr_dRtoa_cloud', 'trendyr_dRtoa_q', 'trendyr_dRtoa_alb', 'trendyr_dRtoa_ts'
        load([dradTrendPath, 'trend_dradEfect_sfc_cld.mat'])% 10 vars:'trendyr_dRsfc_ta','trendyr_dRsfc_taOnly2', 'trendyr_dRsfc_tas2., 'trendyr_dRsfc_tsAtom', 'trendyr_dRsfc_mainEffect', 'trendyr_dRsfc_residual', 'trendyr_dRsfc_cloud', 'trendyr_dRsfc_q', 'trendyr_dRsfc_alb', 'trendyr_dRsfc_ts'
        load([dvarsTrendPath, 'trend_dnetTOA.mat'])% trendyr_dnetTOA
        load([dvarsTrendPath, 'trend_drhs.mat'])% trendyr_drhs
        load([dvarsTrendPath, 'trend_dhFlux.mat'])% trendyr_dhFlux
        load([dvarsTrendPath, 'trend_drhsPlus.mat'])% trendyr_dhFlux
        load([dvarsTrendPath, 'trend_dts.mat'])% trendyr_dts
        load([dnonLocalCldPath, ['trendyr_dRnonLocalCld3_',lower(mlabels.level),'.mat']])% trendyr_dRnonLocalCld3
        load([dnonLocalCldPath, ['trendyr_dTsnonLocalCld3_',lower(mlabels.level),'.mat']])% trendyr_dTsnonLocalCld3
        load([vsTsEffectTrendPath, ['trend_dTs_x_',lower(mlabels.level),'.mat']])% trendm_dTs_alb, trendm_dTs_cld, trendm_dTs_hus, trendm_dTs_ta, trends_dTs_alb, trends_dTs_cld, trends_dTs_hus, trends_dTs_ta, trendyr_dTs_alb, trendyr_dTs_cld, trendyr_dTs_hus, trendyr_dTs_ta
        nlon = length(lon_f); nlat = length(lat_f); 
        
        % cal dRsumNonlocalCldK2 and dRsumNonlocalCld, rad effect
        if strcmp(mlabels.level,'SFC')==1
            trendyr_dRsumNonlocalCld3=trendyr_dRsfc_cloud+trendyr_dRnonLocalCld3;
        elseif strcmp(mlabels.level,'TOA')==1
            trendyr_dRsumNonlocalCld3=trendyr_dRtoa_cloud+trendyr_dRnonLocalCld3;
        end
        
        % cal dTssumNonlocalCld and dTssumNonlocalCldK2, Ts effect
        trendyr_dTssumNonlocalCld3=trendyr_dTs_cld+trendyr_dTsnonLocalCld3;

        % cal dTs_mainEffect
        trendyr_dTs_mainEffect=trendyr_dTs_cld+trendyr_dTs_hus+trendyr_dTs_ta+trendyr_dTs_alb;

        % use one var to plot
        trendyr = zeros(nlon, nlat, 12);
        for jj = 1:3
            eval(['trendyr(:,:,jj)=', mlabels.vars{jj}, ';'])
        end

        trendyr = trendyr * 365 * 10;
        % mask and cal the cc
        [trendyr, yr_cc, yr_pp] = maskArea(trendyr, lat_f, latRange, -latRange, areaNum);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %plot
        f_row = 3; f_col = 1; % 设置画图的行列
        set(0, 'DefaultFigureVisible', 'on')
        ss = get(0, 'ScreenSize'); 
        coef_amplify = 1.5;
        h = figure('Position', [ss(4)/2-100 ss(3) / 35 ss(3)/5*f_col*coef_amplify (ss(4)-80)/5*f_row*coef_amplify]);
        % clf reset;
        set(h, 'Color', [1 1 1]);
        f_matrix = reshape(1:f_row*f_col, [f_col, f_row])';
        % figure
        for varNum = 1:f_row * f_col
            trendz = squeeze(trendyr(:, :, varNum));
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

            title({[mlabels.component{varNum}, mlabels.unite{varNum}]},'Interpreter','latex','fontsize', 10); % cc=',num2str(corr))
            % sequence
            figSeq=char(97+varNum-1);
            yLoc=-0.15;
            if plotRow==f_row
                yLoc=-0.15;
            end
            txt=text(0.5, yLoc, ['(',figSeq,')'],'Units','normalized','FontSize',10);
            hold on
            c = colorbar;
            % c.TickLength = 0.0245;
            % c.Limits = [min_color min_color];
        end

        headLineTxt = {['Level:', mlabels.level, ', Era: ', level.time1{exmName}(1:end - 1),', Trend(year mean)'], ['Model:', level.model2{mdlName}, ', Ensemble: ', esm]};
        sgtt = sgtitle(headLineTxt, 'Fontsize', 14, 'Interpreter', 'none');
        figTitle = [level.time1{exmName}(1:end - 1), '_', level.model2{mdlName}, '_', mlabels.fileN1,'_',esm];
        figurename = [mPath.Output, '/', figTitle, '.png'];
        saveas(gcf, figurename)
        % save_png(figurename)%high resolution
        close gcf
    end

end
