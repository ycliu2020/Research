%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-15 13:43:44
% LastEditTime : 2020-06-16 13:52:53
% LastEditors  : LYC
% Description  :
% FilePath     : /code/p2_processCMIP6Data/s99.plot/s2_plotEachBars_nonLocalCld_dTsTrend.m
%
%%---------------------------------------------------------
%        title(['$k2\cdot \frac{\partial R}{\partial Tsg}$'],'Interpreter','latex','FontSize',20)

clc; clear;
nowpath = pwd;
% load mask map
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')

[mlabels, areaNum] = plotParameters('sfc', 'land', 'nonLocalCld_Ts'); % plot parameters mlabels.level(SFC/TOA)

%% global set
% Experent ID
p1_left = 3; p1_right = 4;
% colorRange ={[-5 -5 -4 -4 -2 -2];[-5 -5 -4 -4 -1 -1];[-3 -4 -1 -2 -.5 -1];[-5 -7 -1.5 -3 -1 -1]};

% Latitude range
latRange = 60;
lon1 = [2.5 357.5]; lat1 = [-latRange + 1 latRange - 1]; % world area
set(0, 'defaultfigurecolor', 'w')

for p_1 = p1_left:p1_right
    [readme, Experiment, level, tLin, mPlev, vars] = modelParameters(p_1);
    % mPath.input:E:/data/cmip6-process/2000-2014/

    mPath.input = fullfile('/data1/liuyincheng/cmip6-process/', level.time1{p_1});
    % mPath.output:a_research/P02.Ts_change_research/figure/04.cmip6Result/2000-2014/
    mPath.uniOutput = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/02.cmip6Result/nonLocalCld/',['dTs_',lower(mlabels.level)], level.time1{p_1});
    mPath.Output = fullfile(mPath.uniOutput);
    auto_mkdir(mPath.Output)

    % load k1_cld, k2_cld, existModelName
    load([mPath.input, 'z_assembleData/non_localCld/k_cld.mat'])
    
    % model loop
    p2_left = 1; p2_right = length(existModelName); % differnt models%length(existModelName)
    for level1 = p2_left:p2_right
        % load data
        varsPath = [mPath.input, existModelName{level1}, '/', level.process3{1}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/rawdata
        dvarsPath = [mPath.input, existModelName{level1}, '/', level.process3{2}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
        dvarsTrendPath = [mPath.input, existModelName{level1}, '/', level.process3{3}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly_trend
        kernelPath = [mPath.input, existModelName{level1}, '/', level.process3{5}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/kernelsCal
        dradEffectPath = [mPath.input, existModelName{level1}, '/', level.process3{6}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
        dradTrendPath = [mPath.input, existModelName{level1}, '/', level.process3{7}]; %/data1/liuyincheng/cmip6-process/amip_1980-2014/CESM2/radEffect_trend
        dnonLocalCldPath = [mPath.input, existModelName{level1}, '/', level.process3{8}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/non_localCld/
        vsTsEffectTrendPath = [mPath.input, existModelName{level1}, '/', level.process3{10}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/vsTsEffect_trend/

        load([dvarsTrendPath, 'global_vars.mat'])% lat lon time plevf readme
        load([dradTrendPath, 'trend_dradEfect_toa_cld.mat'])% 10 vars:'trendyr_dRtoa_ta','trendyr_dRtoa_taOnly2', 'trendyr_dRtoa_tas2., 'trendyr_dRtoa_tsAtom', 'trendyr_dRtoa_mainEffect', 'trendyr_dRtoa_residual', 'trendyr_dRtoa_cloud', 'trendyr_dRtoa_q', 'trendyr_dRtoa_alb', 'trendyr_dRtoa_ts'
        load([dradTrendPath, 'trend_dradEfect_sfc_cld.mat'])% 10 vars:'trendyr_dRsfc_ta','trendyr_dRsfc_taOnly2', 'trendyr_dRsfc_tas2., 'trendyr_dRsfc_tsAtom', 'trendyr_dRsfc_mainEffect', 'trendyr_dRsfc_residual', 'trendyr_dRsfc_cloud', 'trendyr_dRsfc_q', 'trendyr_dRsfc_alb', 'trendyr_dRsfc_ts'
        load([dvarsTrendPath, 'trend_dnetTOA.mat'])% trendyr_dnetTOA
        load([dvarsTrendPath, 'trend_drhs.mat'])% trendyr_drhs
        load([dvarsTrendPath, 'trend_dhFlux.mat'])% trendyr_dhFlux
        load([dvarsTrendPath, 'trend_drhsPlus.mat'])% trendyr_dhFlux
        load([dvarsTrendPath, 'trend_dts.mat'])% trendyr_dts
        load([dnonLocalCldPath, ['trendyr_dRnonLocalCld_',lower(mlabels.level),'.mat']])% trendyr_dRnonLocalCld
        load([dnonLocalCldPath, ['trendyr_dRk2_',lower(mlabels.level),'.mat']])% trendyr_dRk2
        load([dnonLocalCldPath, ['trendyr_dTsnonLocalCld_',lower(mlabels.level),'.mat']])% trendyr_dTsnonLocalCld
        load([dnonLocalCldPath, ['trendyr_dTsk2_',lower(mlabels.level),'.mat']])% trendyr_dTsk2
        load([vsTsEffectTrendPath, ['trend_dTs_x_',lower(mlabels.level),'.mat']])% trendm_dTs_alb, trendm_dTs_cld, trendm_dTs_hus, trendm_dTs_ta, trends_dTs_alb, trends_dTs_cld, trends_dTs_hus, trends_dTs_ta, trendyr_dTs_alb, trendyr_dTs_cld, trendyr_dTs_hus, trendyr_dTs_ta
        nlon = length(lon); nlat = length(lat); 
        % trendyr_dTsnonLocalCld=-trendyr_dTsnonLocalCld;
        % trendyr_dTs_cld=-trendyr_dTs_cld;
        % trendyr_dTs_alb=-trendyr_dTs_alb;
        % trendyr_dTs_hus=-trendyr_dTs_hus;
        % trendyr_dTs_ta=-trendyr_dTs_ta;
        % trendyr_dTs_ta=-trendyr_dTs_ta;
        % cal dRsumNonlocalCldK2 and dRsumNonlocalCld, rad effect
        if strcmp(mlabels.level,'SFC')==1
            trendyr_dRsumNonlocalCld=trendyr_dRsfc_cloud+trendyr_dRnonLocalCld;
        elseif strcmp(mlabels.level,'TOA')==1
            trendyr_dRsumNonlocalCld=trendyr_dRtoa_cloud+trendyr_dRnonLocalCld;
        end
        trendyr_dRsumNonlocalCldK2=trendyr_dRsumNonlocalCld+trendyr_dRk2;
        
        % cal dTssumNonlocalCld and dTssumNonlocalCldK2, Ts effect

        trendyr_dTssumNonlocalCld=trendyr_dTs_cld+trendyr_dTsnonLocalCld;
        trendyr_dTssumNonlocalCldK2=trendyr_dTssumNonlocalCld+trendyr_dTsk2;
        
        % cal dTs_mainEffect
        trendyr_dTs_mainEffect=trendyr_dTs_cld+trendyr_dTs_hus+trendyr_dTs_ta+trendyr_dTs_alb;

        % use one var to plot
        trendyr = zeros(nlon, nlat, 12);
        for jj = 1:12
            eval(['trendyr(:,:,jj)=', mlabels.vars{jj}, ';'])
        end

        trendyr = trendyr * 365 * 10;
        trendyr(:, :, 1) = trendyr(:, :, 1);
        trendyr(:, :, 10) = trendyr(:, :, 10);
        % mask and cal the cc
        [trendyr, yr_cc, yr_pp] = maskArea(trendyr, lat, latRange, -latRange, areaNum);

        set(0, 'DefaultFigureVisible', 'off')
        ss = get(0, 'ScreenSize'); 
        h = figure('Position', [ss(4) / 2 ss(3) / 35 ss(3) * 3/9.5 ss(4) * 4/5]); 
        % set(h,'visible','off');
        % clf reset;
        set(h, 'Color', [1 1 1]);
        f_matrix = reshape(1:12, [3, 4])';
        % figure
        for ii = 1:12
            trendz = squeeze(trendyr(:, :, ii));
            [i, j] = find(f_matrix == ii);
            subplot_yc(4, 3, i, j); 
            hold on
            m_proj('Mercator', 'lon', lon1, 'lat', lat1); %Mercator,Equidistant cylindrical,lambert,Miller Cylindrical
            m_pcolor(lon, lat, trendz');
            % [x, y]=find(trendz(:,:)<=0);
            % m_plot(lon(x),lat(y),'g.','markersize',5,'color','k');
            % [Xlon,Ylat] = meshgrid(lon,lat);
            % m_contourf(Xlon,Ylat,trendz',30,'LineStyle','none')
            colormap(mycolor(18)); %mycolor(100)is soden color????????colormap(flipud(mycolor(13)));%colormap(jet(4))
            col_SeriesNum=10;
            [colorbar_Series] = findSuit_colorInt(trendz, col_SeriesNum);
            % caxis([mmin(ii) mmax(ii)]);
            caxis([min(colorbar_Series) max(colorbar_Series)]);
            hold on
            m_line(world_mapx(:), world_mapy(:), 'color', [0 0 0], 'LineWidth', 0.5);
            if j==1&&i==4
                m_grid('linestyle', 'none', 'tickdir', 'out', 'yaxislocation', 'left', 'fontsize', 8, 'color', 'k'); 
            elseif j==1&&i~=4
                m_grid('linestyle', 'none', 'tickdir', 'out', 'xticklabels',[], 'yaxislocation', 'left', 'fontsize', 8, 'color', 'k');
            elseif j~=1&&i==4 
                m_grid('linestyle', 'none', 'tickdir', 'out', 'yticklabels',[], 'fontsize', 8, 'color', 'k');
            elseif j~=1&&i~=4 
                m_grid('linestyle', 'none', 'tickdir', 'out', 'xticklabels',[], 'yticklabels',[], 'fontsize', 8, 'color', 'k');
            end
            
            title({[mlabels.component{ii}, mlabels.unite{ii}]; ['year mean (cc = ', num2str(yr_cc{ii}), ')']},'Interpreter','latex','fontsize', 10); % cc=',num2str(corr))
            % c=colorbar;
            % % c.Limits=[mmin(ii) mmax(ii)];% 
            % c.Box='off';
            hold on
            c = colorbar;
            % c.TickLength = 0.0245;
            c.Ticks=colorbar_Series(2:end-1);
            c.Limits = [min(colorbar_Series) max(colorbar_Series)];
        end

        tt = ['Level:', mlabels.level, ', Era: ', level.time1{p_1}(1:end - 1), ', Model:', existModelName{level1}];
        sgtt = sgtitle(headLineTxt, 'Fontsize', 14, 'Interpreter', 'none');
        figTitle = [level.time1{p_1}(1:end - 1), '_', existModelName{level1}, '_', mlabels.fileN1];
        figurename = [mPath.Output, '/', figTitle, '.png'];
        saveas(gcf, figurename)
        % save_png(figurename)%high resolution
        close gcf
    end

end
