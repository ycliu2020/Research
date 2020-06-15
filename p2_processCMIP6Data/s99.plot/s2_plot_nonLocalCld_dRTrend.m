%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-15 13:43:44
% LastEditTime : 2020-06-15 19:26:26
% LastEditors  : LYC
% Description  :
% FilePath     : /Research/p2_processCMIP6Data/s99.plot/s2_plot_nonLocalCld_dRTrend.m
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

[mlabels, areaNum] = plotParameters('sfc', 'land', 'non_localCld_radEffect'); % plot parameters mlabels.level(SFC/TOA)

%% global set
% Experent ID
p1_left = 3; p1_right = 3;
% colorRange ={[-5 -5 -4 -4 -2 -2];[-5 -5 -4 -4 -1 -1];[-3 -4 -1 -2 -.5 -1];[-5 -7 -1.5 -3 -1 -1]};
mmin = -2.5; %colorRange{p_1};
mmax = -mmin;
% Latitude range
p_3 = 60;
lon1 = [2.5 357.5]; lat1 = [-p_3 + 1 p_3 - 1]; % world area
set(0, 'defaultfigurecolor', 'w')

for p_1 = p1_left:p1_right
    [readme, Experiment, level, tLin, mPlev, vars] = modelParameters(p_1);
    % mPath.input:E:/data/cmip6-process/2000-2014/

    mPath.input = fullfile('/data1/liuyincheng/cmip6-process/', level.time1{p_1});
    % mPath.output:a_research/P02.Ts_change_research/figure/04.cmip6Result/2000-2014/
    mPath.uniOutput = fullfile('/home/liuyc/research/P02.Ts_change_research/figure/02.cmip6Result/1.4/', lower(mlabels.level), level.time1{p_1});
    mPath.Output = fullfile(mPath.uniOutput);
    auto_mkdir(mPath.Output)

    % load k1_cld, k2_cld, existModelName
    load([inputPath, 'z_assembleData/non_localCld/k_cld.mat'])
    % model loop
    p2_left = 1; p2_right = length(existModelName); % differnt models%length(existModelName)

    for level1 = p2_left:p2_right
        % load data
        varsPath = [inputPath, existModelName{level1}, '/', level.process3{1}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/rawdata
        dvarsPath = [inputPath, existModelName{level1}, '/', level.process3{2}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
        dvarsTrendPath = [inputPath, existModelName{level1}, '/', level.process3{3}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly_trend
        kernelPath = [inputPath, existModelName{level1}, '/', level.process3{5}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/kernelsCal
        dradEffectPath = [inputPath, existModelName{level1}, '/', level.process3{6}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
        dradTrendPath = [mPath.input, existModelName{level1}, '/', level.process3{7}]; %/data1/liuyincheng/cmip6-process/amip_1980-2014/CESM2/radEffect_trend
        dnonLocalCldPath = [inputPath, existModelName{level1}, '/', level.process3{8}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/non_localCld/
        vsTsEffectTrendPath = [inputPath, existModelName{level1}, '/', level.process3{10}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/vsTsEffect_trend/

        load([dradTrendPath, 'global_vars.mat'])% lat lon time plevf readme
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
        load([vsTsEffectTrendPath, ['trendyr_dTs_x_',lower(mlabels.level),'.mat']])% trendm_dTs_alb, trendm_dTs_cld, trendm_dTs_hus, trendm_dTs_ta, trends_dTs_alb, trends_dTs_cld, trends_dTs_hus, trends_dTs_ta, trendyr_dTs_alb, trendyr_dTs_cld, trendyr_dTs_hus, trendyr_dTs_ta
        nlon = length(lon); nlat = length(lat); 
        
        % cal dRsumNonlocalCldK2 and dRsumNonlocalCld, rad effect
        if strcmp(mlabels.level,'SFC')==1
            trendyr_dRsumNonlocalCld=trendyr_dRsfc_cloud+trendyr_dRnonLocalCld;
        elseif strcmp(mlabels.level,'TOA')==1
            trendyr_dRsumNonlocalCld=trendyr_dRtoa_cloud+trendyr_dRnonLocalCld;
        end
        trendyr_dRsumNonlocalCldK2=trendyr_dRsumNonlocalCld+trendyr_dRk2;
        
        % cal dTssumNonlocalCld and dTssumNonlocalCldK2, Ts effect
        if strcmp(mlabels.level,'SFC')==1
            trendyr_dTssumNonlocalCld=trendyr_dTssfc_cloud+trendyr_dTsnonLocalCld;
        elseif strcmp(mlabels.level,'TOA')==1
            trendyr_dTssumNonlocalCld=trendyr_dTstoa_cloud+trendyr_dTsnonLocalCld;
        end
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
        [trendyr, yr_cc, yr_pp] = maskArea(trendyr, lat, p_3, -p_3, areaNum);

        set(0, 'DefaultFigureVisible', 'off')
        ss = get(0, 'ScreenSize'); % ???????????
        h = figure('Position', [ss(4) / 2 ss(3) / 35 ss(3) * 3/9.5 ss(4) * 4/5]); % ???????????????
        % set(h,'visible','off');
        % figure('Position',[ss(4)*2 ss(3)/35 ss(3)*3/9.5 ss(4)*4/5]);   % ???????????????
        % clf reset;
        set(h, 'Color', [1 1 1]);
        f_matrix = reshape(1:12, [3, 4])';
        % figure
        for ii = 1:12
            trendz = squeeze(trendyr(:, :, ii));
            [i, j] = find(f_matrix == ii);
            subplot_yc(4, 3, i, j); % ?????
            hold on
            m_proj('Mercator', 'lon', lon1, 'lat', lat1); %?????????????Mercator,Equidistant cylindrical,lambert,Miller Cylindrical
            m_pcolor(lon, lat, trendz');
            % [x, y]=find(trendz(:,:)<=0);
            % m_plot(lon(x),lat(y),'g.','markersize',5,'color','k');
            % [Xlon,Ylat] = meshgrid(lon,lat);
            % m_contourf(Xlon,Ylat,trendz',30,'LineStyle','none')
            colormap(mycolor(18)); %mycolor(100)is soden color????????colormap(flipud(mycolor(13)));%colormap(jet(4))%????????????????????
            % caxis([mmin(ii) mmax(ii)]);
            caxis([mmin mmax]);
            hold on
            m_line(world_mapx(:), world_mapy(:), 'color', [0 0 0], 'LineWidth', 0.5);
            m_grid('linestyle', 'none', 'tickdir', 'out', 'yaxislocation', 'left', 'fontsize', 10, 'color', 'k'); %????????????
            title({[mlabels.component{ii}, mlabels.unite{ii}]; ['year mean (cc = ', num2str(yr_cc{ii}), ')']}); % cc=',num2str(corr))
            % c=colorbar;
            % % c.Limits=[mmin(ii) mmax(ii)];
            % c.Box='off';
            hold on
            pos = get(gca, 'Position');

        end

        c = colorbar('southoutside', 'Position', [pos(1) - 0.467 pos(2) - 0.05 0.6 pos(4) / 8]);
        c.TickLength = 0.0245;

        if mmin == -5
            c.Ticks = [-4, -3, -2, -1, 0, 1, 2, 3, 4];
        elseif mmin == -4
            c.Ticks = [-3.2, -2.4, -1.6, -.8, 0, .8, 1.6, 2.4, 3.2];
        elseif mmin == -3
            c.Ticks = [-2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5];
        elseif mmin == -2.5
            c.Ticks = [-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2];
        end

        % c.LineWidth = 2;
        c.Limits = [mmin mmax];
        % c.LineWidth = 'white';
        % c.Box='off';

        tt = ['Level:', mlabels.level, ', Era: ', level.time1{p_1}(1:end - 1), ', Model:', existModelName{level1}];
        sgtt = sgtitle(tt, 'Fontsize', 14, 'Interpreter', 'none');
        f_tt = [level.time1{p_1}(1:end - 1), '_', existModelName{level1}, '_', mlabels.fileN1];
        % figurename = [mPath.Output, '/', f_tt, '.png'];
        % saveas(gcf, figurename)
        % % save_png(figurename)%high resolution
        % close gcf
    end

end
