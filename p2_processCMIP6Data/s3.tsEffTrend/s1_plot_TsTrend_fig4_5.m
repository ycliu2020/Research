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
[mlabels, areaNum] = cmipPlotParameters('sfc', 'land', 'radEffect'); % plot parameters

%% global set
% Latitude range
latRange = 60;
lon1 = [2.5 357.5]; lat1 = [-latRange + 1 latRange - 1]; % world area

% experiment ID
exmStart = 4; exmEnd = 4;

for exmNum = exmStart:exmEnd

    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % mPath.input:E:/data/cmip6-process/2000-2014/
    exmPath = fullfile('/data1/liuyincheng/cmip6-process/', level.time1{exmNum});
    % mPath.output:a_research/P02.Ts_change_research/figure/04.cmip6Result/2000-2014/
    mPath.uniOutput = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/proj2_cmip6Result/Ts_effect/', level.time1{exmNum});
    mPath.Output = fullfile(mPath.uniOutput);
    auto_mkdir(mPath.Output)

    % model loop
    mdlStart = 9; mdlEnd = 9;%length(level.model2); % differnt models%length(level.model2)

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
            %load data
            varsPath = fullfile(esmPath, level.process3{1}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/rawdata
            dvarsPath = fullfile(esmPath, level.process3{2}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
            dvarsTrendPath = fullfile(esmPath, level.process3{3}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly_trend
            kernelPath = fullfile(esmPath, level.process3{5}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/kernelsCal
            dradEffectPath = fullfile(esmPath, level.process3{6}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
            dradTrendPath = fullfile(esmPath, level.process3{7}); %/data1/liuyincheng/cmip6-process/amip_1980-2014/CESM2/radEffect_trend
            dnonLocalCldPath = fullfile(esmPath, level.process3{8}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/non_localCld/
            vsTsEffectTrendPath = fullfile(esmPath, level.process3{10}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/vsTsEffect_trend/

            load([dradTrendPath, 'global_vars.mat'])% lat lon time plevf readme
            load([dradTrendPath, 'trend_dradEfect_toa_cld.mat'])% 10 vars:'trendyr_dRtoa_ta','trendyr_dRtoa_taOnly2', 'trendyr_dRtoa_tas2., 'trendyr_dRtoa_tsAtom', 'trendyr_dRtoa_mainEffect', 'trendyr_dRtoa_residual', 'trendyr_dRtoa_cloud', 'trendyr_dRtoa_q', 'trendyr_dRtoa_alb', 'trendyr_dRtoa_ts'
            load([dradTrendPath, 'trend_dradEfect_sfc_cld.mat'])% 10 vars:'trendyr_dRsfc_ta','trendyr_dRsfc_taOnly2', 'trendyr_dRsfc_tas2., 'trendyr_dRsfc_tsAtom', 'trendyr_dRsfc_mainEffect', 'trendyr_dRsfc_residual', 'trendyr_dRsfc_cloud', 'trendyr_dRsfc_q', 'trendyr_dRsfc_alb', 'trendyr_dRsfc_ts'
            load([dvarsTrendPath, 'trend_dnetTOA.mat'])% trendyr_dnetTOA
            load([dvarsTrendPath, 'trend_drhs.mat'])% trendyr_drhs
            load([dvarsTrendPath, 'trend_dhFlux.mat'])% trendyr_dhFlux
            load([dvarsTrendPath, 'trend_drhsPlus.mat'])% trendyr_dhFlux
            load([dvarsTrendPath, 'trend_dts.mat'])% trendyr_dts
            load([vsTsEffectTrendPath, ['trend_dTs_x_', lower(mlabels.level), '.mat']])% trendm_dTs_alb, trendm_dTs_cld, trendm_dTs_hus, trendm_dTs_ta, trends_dTs_alb, trends_dTs_cld, trends_dTs_hus, trends_dTs_ta, trendyr_dTs_alb, trendyr_dTs_cld, trendyr_dTs_hus, trendyr_dTs_ta
            nlon = length(lon_f); nlat = length(lat_f);

            % use one var to plot
            trendyr = zeros(nlon, nlat, 6);
            trendyr(:,:,1)=trendyr_dts;
            trendyr(:,:,2)=trendyr_dTs_cld;
            trendyr(:,:,3)=trendyr_dTs_ta;
            trendyr(:,:,4)=trendyr_dTs_hus;
            trendyr(:,:,5)=trendyr_dTs_alb;
            trendyr(:,:,6)=trendyr_dTs_residual;
            mlabels.component={'Total surface temperature change', 'Cloud contribution', 'Air temperature contribution', 'Hus contribution', 'Albedo contribution', 'Residual contribution'};
            trendyr = trendyr * 365 * 10;
            % mask and cal the cc
            [trendyr, yr_cc, yr_pp] = maskArea(trendyr, lat_f, latRange, -latRange, areaNum);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %plot
            f_row = 3; f_col = 2; % 设置画图的行列
            set(0, 'DefaultFigureVisible', 'on')
            ss = get(0, 'ScreenSize'); %
            h = figure('Position', [ss(4) / 2 - 100 ss(3) / 35 ss(3) / 5 * f_col (ss(4) - 80) / 5 * f_row]);
            % set(h,'visible','off');
            % figure('Position',[ss(4)*2 ss(3)/35 ss(3)*3/9.5 ss(4)*4/5]);   %
            % clf reset;
            set(h, 'Color', [1 1 1]);
            f_matrix = reshape(1:f_row*f_col, [f_col, f_row])';
            if exmNum == 1 || exmNum == 2
                mmin = -1; 
            elseif exmNum == 3 || exmNum == 4
                mmin = -0.7; 
            end
            mmax = -mmin;

            % figure
            for varNum = 1:f_row * f_col
                trendz = squeeze(trendyr(:, :, varNum));
                [plotRow, plotCol] = find(f_matrix == varNum);
                subplot_yc(f_row, f_col, plotRow, plotCol);
                hold on
                m_proj('Mercator', 'lon', lon1, 'lat', lat1); %Mercator,Equidistant cylindrical,lambert,Miller Cylindrical
                m_pcolor(lon_f, lat_f, trendz');
                colormap(mycolor(18)); %mycolor(100)is soden color????????colormap(flipud(mycolor(13)));%colormap(jet(4))
                caxis([mmin mmax]);
                hold on
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
                
                title(mlabels.component{varNum}, 'Interpreter', 'latex', 'fontsize', 12); % cc=',num2str(corr))
                % c=colorbar;
                % % c.Limits=[mmin(varNum) mmax(varNum)];
                % c.Box='off';
                hold on
                pos = get(gca, 'Position');

            end

            c = colorbar('southoutside', 'Position', [pos(1)-0.35  pos(2) - 0.05 0.6 pos(4) / 8]);
            c.TickLength = 0.032;
            c.Ticks = mmax*[-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8];

            % c.LineWidth = 2;
            c.Limits = [mmin mmax];
            % c.LineWidth = 'white';
            % c.Box='off';

            headLineTxt = {['Level:', mlabels.level, ', Era: ', level.time1{exmNum}(1:end - 1), ', Trend(year mean)'], ['Model:', level.model2{mdlNum}, ', Ensemble: ', esm, ', Unit: K/10a']};
            sgtt = sgtitle(headLineTxt, 'Fontsize', 14, 'Interpreter', 'none');
            figTitle = [level.time1{exmNum}(1:end - 1), '_', level.model2{mdlNum}, '_', mlabels.fileN1, '_', esm];
            figurename = [mPath.Output, '/', figTitle, '.png'];
            saveas(gcf, figurename)
            % save_png(figurename)%high resolution
            % close gcf

            disp([esmName{esmNum, 1}, ' ensemble is done!'])

        end

        disp([level.model2{mdlNum}, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{exmNum}, ' era is done!'])
    disp(' ')
end

eval(['cd ', nowpath]);
