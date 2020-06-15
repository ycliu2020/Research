%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% mothly data
% plot the Ts trend and rhs,Ts, Ta, wv, albedo trend(TOA)
% note: the trendyr already be masked and caled cc
%
% experiment information:
% time:2000.01-2014.12(interval:15*12);1980.01-2014.12(interval:35*12); 2015.01-2099.12(interval:85*12)
% initial time in hist(1740 total): 1,561 of 1740(2000.03);1,321 of 1740(1980.01)
% initial time in futrue(1032 total): 1 of 1032(2015.01);
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clc; clear;
nowpath = pwd;
% load mask map
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')
[mlabels, areaNum] = plotParameters('toa', 'land', 'cld-clr_nodp'); % plot parameters 

%% global set
% Experent ID
p1_left = 5; p1_right = 5; 
% colorRange ={[-5 -5 -4 -4 -2 -2];[-5 -5 -4 -4 -1 -1];[-3 -4 -1 -2 -.5 -1];[-5 -7 -1.5 -3 -1 -1]};
mmin = -3; %colorRange{p_1};
mmax = -mmin;
% Latitude range
p_3 = 60;  
lon1 = [2.5 357.5]; lat1 = [-p_3+1 p_3-1]; % world area
set(0, 'defaultfigurecolor', 'w')

for p_1 = p1_left:p1_right
    [readme, Experiment, level, tLin, mPlev, vars] = modelParameters(p_1);
    % mPath.input:E:/data/cmip6-process/2000-2014/
    mPath.input = fullfile('/data1/liuyincheng/cmip6-process/', mlabels.algmPath, level.time1{p_1});
    % mPath.output:a_research/P02.Ts_change_research/figure/04.cmip6Result/2000-2014/
    mPath.uniOutput = ['/home/liuyc/research/P02.Ts_change_research/figure/04.cmip6Result/1.3/toa/', level.time1{p_1}];
    mPath.Output = fullfile(mPath.uniOutput, mlabels.fileN1);
    auto_mkdir(mPath.Output)

    % model loop
    p2_left = 1; p2_right = 1;%length(level.model2);% differnt models%length(level.model2)
    for level1 = p2_left:p2_right
        % load data
        dradTrendPath = [mPath.input, level.model2{level1}, '/', level.process3{7}]; %/data1/liuyincheng/cmip6-process/amip_1980-2014/CESM2/Effect_trend
        danomTrendPath = [mPath.input, level.model2{level1}, '/', level.process3{3}]; %/data1/liuyincheng/cmip6-process/amip_1980-2014/CESM2/anomaly_trend
        load([dradTrendPath, 'global_vars.mat'])% lat lon time plevf readme
        load([dradTrendPath, 'trend_dradEfect_toa_cld.mat'])% 10 vars:'trendyr_dRtoa_ta','trendyr_dRtoa_taOnly2', 'trendyr_dRtoa_tas2., 'trendyr_dRtoa_tsAtom', 'trendyr_dRtoa_mainEffect', 'trendyr_dRtoa_residual', 'trendyr_dRtoa_cloud', 'trendyr_dRtoa_q', 'trendyr_dRtoa_alb', 'trendyr_dRtoa_ts'
        load([dradTrendPath, 'trend_dradEfect_sfc_cld.mat'])% 10 vars:'trendyr_dRsfc_ta','trendyr_dRsfc_taOnly2', 'trendyr_dRsfc_tas2., 'trendyr_dRsfc_tsAtom', 'trendyr_dRsfc_mainEffect', 'trendyr_dRsfc_residual', 'trendyr_dRsfc_cloud', 'trendyr_dRsfc_q', 'trendyr_dRsfc_alb', 'trendyr_dRsfc_ts'
        load([danomTrendPath, 'trend_dnetTOA.mat'])% trendyr_dnetTOA
        load([danomTrendPath, 'trend_drhs.mat'])% trendyr_drhs
        load([danomTrendPath, 'trend_dhFlux.mat'])% trendyr_dhFlux
        load([danomTrendPath, 'trend_dts.mat'])% trendyr_dts
        nlon = length(lon); nlat = length(lat);
        % use one var to plot
        trendyr = zeros(nlon, nlat, 12);
        for jj = 1:12
            eval(['trendyr(:,:,jj)=', mlabels.vars{jj}, ';'])
        end
        trendyr = trendyr * 365 * 10;
        trendyr(:, :, 1) = trendyr(:, :, 1) .* 5;
        trendyr(:, :, 10) = trendyr(:, :, 10);
        % mask and cal the cc
        [trendyr, yr_cc, yr_pp] = maskArea(trendyr, lat, p_3, -p_3, areaNum);

        set(0, 'DefaultFigureVisible', 'off')
        ss = get(0, 'ScreenSize'); % ???????????
        h = figure('Position', [ss(4)/2 ss(3)/35 ss(3)*3/9.5 ss(4)*4/5]); % ???????????????
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
            colormap(mycolor(100)); %????????colormap(flipud(mycolor(13)));%colormap(jet(4))%????????????????????
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

        c = colorbar('southoutside', 'Position', [pos(1)-0.467 pos(2)-0.05 0.6 pos(4)/8]);
        c.TickLength = 0.02;
        c.Ticks = [mmin, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, mmax];
        % c.LineWidth = 2;
        c.Limits = [mmin mmax];
        % c.LineWidth = 'white';
        % c.Box='off';

        tt = ['Level:', mlabels.level, ', Era: ', level.time1{p_1}(1:end - 1), ', Model:', level.model2{level1}];
        sgtt = sgtitle(tt, 'Fontsize', 14, 'Interpreter', 'none');
        f_tt = [level.time1{p_1}(1:end - 1), '_', level.model2{level1},'_',mlabels.fileN1];
        figurename=[mPath.Output,'/',f_tt,'.png'];
        saveas(gcf,figurename)
        % save_png(figurename)%high resolution
        close gcf
    end

end
