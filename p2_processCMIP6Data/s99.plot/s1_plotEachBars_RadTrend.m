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
[mlabels, areaNum] = cmipPlotParameters('toa', 'landsea', 'radEffect'); % plot parameters 
esm = 'r1i1p1f1';
%% global set
% Latitude range
latRange = 60;  
lon1 = [2.5 357.5]; lat1 = [-latRange+1 latRange-1]; % world area
set(0, 'defaultfigurecolor', 'w')
% Experent ID
exm_left = 4; exm_right = 4; 
for exmName = exm_left:exm_right
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmName);
    % mPath.input:E:/data/cmip6-process/2000-2014/
    mPath.input = fullfile('/data1/liuyincheng/cmip6-process/', level.time1{exmName});
    % mPath.output:a_research/P02.Ts_change_research/figure/04.cmip6Result/2000-2014/
    mPath.uniOutput = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/02.cmip6Result/1.5/',lower(mlabels.level),level.time1{exmName});
    mPath.Output = fullfile(mPath.uniOutput);
    auto_mkdir(mPath.Output)

    % model loop
    mdl_left = 1; mdl_right = 1;%length(level.model2);% differnt models%length(level.model2)
    for mdlName = mdl_left:mdl_right
        mdlPath=fullfile(mPath.input, level.model2{mdlName});

        % load data
        dradTrendPath = fullfile(mdlPath, esm, level.process3{7}); %/data1/liuyincheng/cmip6-process/amip_1980-2014/CESM2/Effect_trend
        danomTrendPath = fullfile(mdlPath, esm, level.process3{3}); %/data1/liuyincheng/cmip6-process/amip_1980-2014/CESM2/anomaly_trend
        if ~exist(dradTrendPath,'dir')
            disp(['the ', esm, ' ensemble of ',level.model2{mdlName}, ' didnt exist']);
            continue
        end
        load([dradTrendPath, 'global_vars.mat'])% lat lon time plevf readme
        load([dradTrendPath, 'trend_dradEfect_toa_cld.mat'])% 10 vars:'trendyr_dRtoa_ta','trendyr_dRtoa_taOnly2', 'trendyr_dRtoa_tas2., 'trendyr_dRtoa_tsAtom', 'trendyr_dRtoa_mainEffect', 'trendyr_dRtoa_residual', 'trendyr_dRtoa_cloud', 'trendyr_dRtoa_q', 'trendyr_dRtoa_alb', 'trendyr_dRtoa_ts'
        load([dradTrendPath, 'trend_dradEfect_sfc_cld.mat'])% 10 vars:'trendyr_dRsfc_ta','trendyr_dRsfc_taOnly2', 'trendyr_dRsfc_tas2., 'trendyr_dRsfc_tsAtom', 'trendyr_dRsfc_mainEffect', 'trendyr_dRsfc_residual', 'trendyr_dRsfc_cloud', 'trendyr_dRsfc_q', 'trendyr_dRsfc_alb', 'trendyr_dRsfc_ts'
        load([danomTrendPath, 'trend_dnetTOA.mat'])% trendyr_dnetTOA
        load([danomTrendPath, 'trend_drhs.mat'])% trendyr_drhs
        load([danomTrendPath, 'trend_dhFlux.mat'])% trendyr_dhFlux
        load([danomTrendPath, 'trend_drhsPlus.mat'])% trendyr_dhFlux
        load([danomTrendPath, 'trend_dts.mat'])% trendyr_dts
        nlon = length(lon); nlat = length(lat);
        % use one var to plot
        trendyr = zeros(nlon, nlat, 12);
        for varNum = 1:12
            eval(['trendyr(:,:,varNum)=', mlabels.vars{varNum}, ';'])
        end
        trendyr = trendyr * 365 * 10;
        trendyr(:, :, 1) = trendyr(:, :, 1);
        trendyr(:, :, 10) = trendyr(:, :, 10);
        % mask and cal the cc
        [trendyr, yr_cc, yr_pp] = maskArea(trendyr, lat, latRange, -latRange, areaNum);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %plot
        f_row = 4; f_col = 3; % 设置画图的行列
        set(0, 'DefaultFigureVisible', 'on')
        ss = get(0, 'ScreenSize'); 
        h = figure('Position', [ss(4)/2-100 ss(3) / 35 ss(3)/5*f_col (ss(4)-80)/5*f_row]);
        % set(h,'visible','off');
        % figure('Position',[ss(4)*2 ss(3)/35 ss(3)*3/9.5 ss(4)*4/5]);   
        % clf reset;
        set(h, 'Color', [1 1 1]);
        f_matrix = reshape(1:12, [3, 4])';
        % figure
        for varNum = 1:f_row * f_col
            trendz = squeeze(trendyr(:, :, varNum));
            [plotRow, plotCol] = find(f_matrix == varNum);
            subplot_yc(4, 3, plotRow, plotCol); % ?????
            hold on
            m_proj('Mercator', 'lon', lon1, 'lat', lat1); 
            m_pcolor(lon, lat, trendz');
            % [x, y]=find(trendz(:,:)<=0);
            % m_plot(lon(x),lat(y),'g.','markersize',5,'color','k');
            % [Xlon,Ylat] = meshgrid(lon,lat);
            % m_contourf(Xlon,Ylat,trendz',30,'LineStyle','none')
            colormap(mycolor(18)); %mycolor(100)is soden color????????colormap(flipud(mycolor(13)));%colormap(jet(4))
            col_SeriesNum = 10;
            [colorbar_Series] = findSuit_colorInt(trendz, col_SeriesNum);
            % caxis([mmin(varNum) mmax(varNum)]);
            caxis([min(colorbar_Series) max(colorbar_Series)]);
            hold on
            % plot ticklabels only on the margin
            m_line(world_mapx(:), world_mapy(:), 'color', [0 0 0], 'LineWidth', 0.5);
            if plotCol == 1 && plotRow == f_row
                m_grid('linestyle', 'none', 'tickdir', 'out', 'yaxislocation', 'left', 'fontsize', 8, 'color', 'k');
            elseif plotCol == 1 && plotRow ~= f_row
                m_grid('linestyle', 'none', 'tickdir', 'out', 'xticklabels', [], 'yaxislocation', 'left', 'fontsize', 8, 'color', 'k');
            elseif plotCol ~= 1 && plotRow == f_row
                m_grid('linestyle', 'none', 'tickdir', 'out', 'yticklabels', [], 'fontsize', 8, 'color', 'k');
            elseif plotCol ~= 1 && plotRow ~= f_row
                m_grid('linestyle', 'none', 'tickdir', 'out', 'xticklabels', [], 'yticklabels', [], 'fontsize', 8, 'color', 'k');
            end            
            title({[mlabels.component{varNum}, mlabels.unite{varNum}]; ['spatial cc = ', num2str(yr_cc{varNum})]}, 'Interpreter', 'latex', 'fontsize', 10); % cc=',num2str(corr))
            % c=colorbar;
            % % c.Limits=[mmin(varNum) mmax(varNum)];
            % c.Box='off';
            hold on
            c = colorbar;
            % c.TickLength = 0.0245;
            c.Ticks = colorbar_Series(2:end - 1);
            c.Limits = [min(colorbar_Series) max(colorbar_Series)];

        end
        headLine = {['Level:', mlabels.level, ', Era: ', level.time1{exmName}(1:end - 1), ', Trend(year mean)'],['Model:', level.model2{mdlName}, ', Ensemble: ', esm]};
        sgtt = sgtitle(headLine, 'Fontsize', 14, 'Interpreter', 'none');
        figTitle = [level.time1{exmName}(1:end - 1), '_', level.model2{mdlName},'_',mlabels.fileN1,'_',esm];
        figurename=[mPath.Output,'/',figTitle,'.png'];
        % saveas(gcf,figurename)
        % % save_png(figurename)%high resolution
        % close gcf
    end

end


% old setting of the color range
% pos = get(gca, 'Position');
% c = colorbar('southoutside', 'Position', [pos(1)-0.467 pos(2)-0.05 0.6 pos(4)/8]);
% if mmin == -5
%     c.Ticks = [-4, -3, -2, -1, 0, 1, 2, 3, 4];
% elseif mmin == -4
%     c.Ticks = [-3.2, -2.4, -1.6, -.8, 0, .8, 1.6, 2.4, 3.2];
% elseif mmin == -3
%     c.Ticks = [-2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5];
% elseif mmin == -2.5
%     c.Ticks = [-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2];
% end
% % c.LineWidth = 2;
% c.Limits = [mmin mmax];
% % c.LineWidth = 'white';
% % c.Box='off';

