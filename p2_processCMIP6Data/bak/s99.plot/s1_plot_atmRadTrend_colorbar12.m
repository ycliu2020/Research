%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% mothly data
% plot the Ts trend and rhs,Ts, Ta, wv, albedo trend(ATM)
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

for p_1 = 1:2
    [readme, Experiment, level, tLin, mPlev, vars] = modelParameters(p_1);
    % inputPath:E:/data/cmip6-process/2000-2014/
    inputPath = ['/data1/liuyincheng/cmip6-process/', level.time1{p_1}];
    % outputPath:a_research/P02.Ts_change_research/figure/04.cmip6Result/2000-2014/
    outputPath = ['/home/lyc/research/P02.Ts_change_research/figure/04.cmip6Result/1.4/atm/', level.time1{p_1}];
    auto_mkdir(outputPath)
    load('/home/lyc/lib/tools/matlab/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
    load('/home/lyc/lib/tools/matlab/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
    load('/home/lyc/lib/tools/matlab/myMap/02.world_map/mat_file/correct_worldmap.mat')%
    load('/home/lyc/lib/tools/matlab/myMap/01.china_map/mat_file/mask14472.mat')
    
    % labels
    month = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
    season = {'MAM', 'JJA', 'SON', 'DJF'};
    unite = {' trend(0.2K/10a)', ' trend(W\cdotm-2/10a)', ' trend(W\cdotm-2/10a)', ...
        ' trend(W\cdotm-2/10a)', ' trend(W\cdotm-2/10a)', ' trend(W\cdotm-2/10a)', ...
        ' trend(W\cdotm-2/10a)', ' trend(W\cdotm-2/10a)', ' trend(W\cdotm-2/10a)', ...
        ' trend(W\cdotm-2/10a)', ' trend(W\cdotm-2/10a)', ' trend(W\cdotm-2/10a)'};
    component = {'Ts', 'RHeating+LH+SH(sfc)', 'RHeating(sfc)', ...
        'Ta Effect', 'Ta Effect(near sfc)', 'Ta Effect(without near sfc) ', ...
        'Ts Effect on atmos', 'Ta+q+Alb+Cloud Effect', 'Kernel Residual', ...
        'Cloud Effect', 'q Effect', 'Albedo Effect'}; %4*3
    component_labelToa = {'dts', 'dnetTOA', 'drhs', ...
        'dRtoa_ta', 'dRtoa_tas2', 'dRtoa_taOnly2', ...
        'dRatm_tsAtom', 'dRtoa_mainEffect', 'dRtoa_residual', ...
        'dRtoa_cloud', 'dRtoa_q', 'dRtoa_alb', 'dRtoa_ts'};
    component_labelSfc = {'dts', 'drhsPlus', 'drhs', ...
        'dRsfc_ta', 'dRsfc_tas2', 'dRsfc_taOnly2', ...
        'dRatm_tsAtom', 'dRsfc_mainEffect', 'dRsfc_residual', ...
        'dRsfc_cloud', 'dRsfc_q', 'dRsfc_alb', 'dRsfc_ts'};
    vars_labelToa = strcat('trendyr_', component_labelToa);
    vars_labelSfc = strcat('trendyr_', component_labelSfc);

    p_2 = 3;
    level_label = {'SFC', 'TOA', 'ATM(TOA-SFC)'};
    % figure set
    set(0, 'defaultfigurecolor', 'w')
    areaNum = 1; % world land
    latRange = 60; % Latitude range
    lon1 = [2.5 357.5]; lat1 = [-latRange+1 latRange-1]; % world area
    % colorRange={[-5 -5 -4 -4 -2 -2];[-5 -5 -4 -4 -1 -1];[-3 -4 -1 -2 -.5 -1];[-5 -7 -1.5 -3 -1 -1]};
    mmin = -4; %colorRange{p_1};
    mmax = -mmin;

    % model loop
    for level1 = 1:length(level.model2)
        % load data
        dradTrendPath = [inputPath, level.model2{level1}, '/', level.process3{7}]; %/data1/liuyincheng/cmip6-process/amip_1980-2014/CESM2/Effect_trend
        danomTrendPath = [inputPath, level.model2{level1}, '/', level.process3{3}]; %/data1/liuyincheng/cmip6-process/amip_1980-2014/CESM2/anomaly_trend
        load([dradTrendPath, 'global_vars.mat'])% lat lon time plevf readme
        load([dradTrendPath, 'trend_dradEfect_toa_cld.mat'])% 10 vars:'trendyr_dRtoa_ta','trendyr_dRtoa_taOnly2', 'trendyr_dRtoa_tas2., 'trendyr_dRtoa_tsAtom', 'trendyr_dRtoa_mainEffect', 'trendyr_dRtoa_residual', 'trendyr_dRtoa_cloud', 'trendyr_dRtoa_q', 'trendyr_dRtoa_alb', 'trendyr_dRtoa_ts'
        load([dradTrendPath, 'trend_dradEfect_sfc_cld.mat'])% 10 vars:'trendyr_dRsfc_ta','trendyr_dRsfc_taOnly2', 'trendyr_dRsfc_tas2., 'trendyr_dRsfc_tsAtom', 'trendyr_dRsfc_mainEffect', 'trendyr_dRsfc_residual', 'trendyr_dRsfc_cloud', 'trendyr_dRsfc_q', 'trendyr_dRsfc_alb', 'trendyr_dRsfc_ts'
        load([danomTrendPath, 'trend_dnetTOA.mat'])% trendyr_dnetTOA
        load([danomTrendPath, 'trend_drhs.mat'])% trendyr_drhs
        load([danomTrendPath, 'trend_drhsPlus.mat'])% trendyr_drhs
        load([danomTrendPath, 'trend_dhFlux.mat'])% trendyr_dhFlux
        load([danomTrendPath, 'trend_dts.mat'])% trendyr_dts
        nlon = length(lon); nlat = length(lat);
        % use one var to plot
        trendyr = zeros(nlon, nlat, 12);

        for jj = 1:12
            eval(['trendyrToa(:,:,jj)=', vars_labelToa{jj}, ';'])
            eval(['trendyrSfc(:,:,jj)=', vars_labelSfc{jj}, ';'])
        end

        trendyrToa = trendyrToa * 365 * 10;
        trendyrToa(:, :, 1) = trendyrToa(:, :, 1) .* 5;

        trendyrSfc = trendyrSfc * 365 * 10;
        trendyrSfc(:, :, 1) = trendyrSfc(:, :, 1) .* 5;
        trendyrAtm = trendyrToa - trendyrSfc;
        trendyrAtm(:, :, [1 2 3 7]) = trendyrSfc(:, :, [1 2 3 7]);
        % mask and cal the cc
        [trendyrAtm, yr_cc, yr_pp] = maskArea(trendyrAtm, lat, latRange, -latRange, areaNum);

        set(0, 'DefaultFigureVisible', 'off')
        ss = get(0, 'ScreenSize'); % ???????????
        h = figure('Position', [ss(4)/2 ss(3)/35 ss(3)*3/9.5 ss(4)*4/5]); % ???????????????
        % figure('Position',[ss(4)*2 ss(3)/35 ss(3)*3/9.5 ss(4)*4/5]);   % ???????????????
        clf reset;
        set(h, 'Color', [1 1 1]);
        f_matrix = reshape(1:12, [3, 4])';
        % figure
        for ii = 1:12
            trendz = squeeze(trendyrAtm(:, :, ii));
            [i, j] = find(f_matrix == ii);

            subplot_yc(4, 3, i, j); % ?????
            hold on
            m_proj('Mercator', 'lon', lon1, 'lat', lat1); %?????????????Mercator,Equidistant cylindrical,lambert,Miller Cylindrical
            m_pcolor(lon, lat, trendz');
            %contourf
            % [Xlon,Ylat] = meshgrid(lon,lat);
            % m_contourf(Xlon,Ylat,trendz',30,'LineStyle','none')
            colormap(mycolor(18)); %????????colormap(flipud(mycolor(13)));%colormap(jet(4))%????????????????????
            % caxis([mmin(ii) mmax(ii)]);
            caxis([mmin mmax]);
            hold on
            m_line(world_mapx(:), world_mapy(:), 'color', [0 0 0], 'LineWidth', 0.5);
            m_grid('linestyle', 'none', 'tickdir', 'out', 'yaxislocation', 'left', 'fontsize', 10, 'color', 'k'); %????????????
            title({[component{ii}, unite{ii}]; ['year mean (cc = ', num2str(yr_cc{ii}), ')']}); % cc=',num2str(corr))
            % c=colorbar;
            % % c.Limits=[mmin(ii) mmax(ii)];
            % c.Box='off';
            hold on
            pos = get(gca, 'Position');

        end

        c = colorbar('southoutside', 'Position', [pos(1)-0.467 pos(2)-0.05 0.6 pos(4)/8]);
        c.TickLength = 0.0245;
        if mmin==-5
            c.Ticks = [ -4, -3, -2, -1, 0, 1, 2, 3, 4];
        elseif mmin==-4
            c.Ticks = [ -3.2, -2.4, -1.6, -.8, 0, .8, 1.6, 2.4, 3.2];
        end
        % c.LineWidth = 2;
        c.Limits = [mmin mmax];
        % c.LineWidth = 'white';
        % c.Box='off';

        tt = ['Level:', level_label{p_2}, ', Era: ', level.time1{p_1}(1:end - 1), ', Model:', level.model2{level1}];
        sgtt = sgtitle(headLineTxt, 'Fontsize', 14, 'Interpreter', 'none');
        figTitle = [level.time1{p_1}(1:end - 1), '_', level.model2{level1}];
        figurename = [outputPath, figTitle, '.png'];
        saveas(gcf, figurename)
        % save_png(figurename)%high resolution
        close gcf
    end

end
