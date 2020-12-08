%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-08-31 17:00:15
% LastEditTime : 2020-11-25 17:04:57
% LastEditors  : Please set LastEditors
% Description  : 同时画时间序列和相关性分布图
% FilePath     : /code/p3_paperFigIntegrate/Fig1_2_radClosure/CMIP6_plot_RadClos.m
%
%%---------------------------------------------------------
clc; clear; tic;
nowpath = pwd;
dbstop if error

% load mask map
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')

latRange = 90; % Latitude range
lon1 = [2.5 357.5]; lat1 = [-latRange + 1 latRange - 1]; % world area
toaSfc = {'toa', 'sfc'};
lon_k = 0:2.5:357.5; nlonk = length(lon_k); % kernel lat lon
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment
for exmNum = 1:1%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    %CAMS-CSM1-0 didn't have sfc clear sky radiation, delete it
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % exmPath
    exmPath = ['/data1/liuyincheng/CMIP6-process/', level.time1{exmNum}]; %/data1/liuyincheng/cmip6-process/2000-2014/
    mPath.uniOutput = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/proj3_PaperFig/v0.0/Fig2_CMIP6_200003-201402_world_radClos'); %['dRTs_', lower(mlabels.level)],
    mPath.Output = fullfile(mPath.uniOutput);
    auto_mkdir(mPath.Output)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model
    for mdlNum = 1:length(level.model2)
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
            % data path
            varsPath = fullfile(esmPath, level.process3{1}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/rawdata_regrid
            dvarsPath = fullfile(esmPath, level.process3{2}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
            dvarsTrendPath = fullfile(esmPath, level.process3{3}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly_trend
            kernelPath = fullfile(esmPath, level.process3{5}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/kernelsCal
            dradEffectPath = fullfile(esmPath, level.process3{6}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
            dnonLocalCldPath = fullfile(esmPath, level.process3{8}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/non_localCld/
            dvsTsEffectPath = fullfile(esmPath, level.process3{9}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/vsTsEffect/
            load([dradEffectPath, 'global_vars.mat'])% lat_f lon_f time plev_k readme
            ntime = length(time.date);

            %% model dradEffect
            % (lat, lon, time, net/long/short)
            dRclr_sfc = zeros(144, 72, ntime, 3);
            dRclr_toa = zeros(144, 72, ntime, 3);
            load([dradEffectPath, 'model_dradEffect.mat'])%dR_net_cld, dR_net_clr, lw_net, readme_realradEfect, sw_net
            dRclr_sfc(:, :, :, 1) = dR_net_clr(:, :, :, 1);
            dRclr_toa(:, :, :, 1) = dR_net_clr(:, :, :, 2);
            clear dR_net_cld dR_net_clr
            dRclr_sfc(:, :, :, 2) = lw_net(:, :, :, 2);
            dRclr_toa(:, :, :, 2) = lw_net(:, :, :, 4);
            clear lw_net
            dRclr_sfc(:, :, :, 3) = sw_net(:, :, :, 2);
            dRclr_toa(:, :, :, 3) = sw_net(:, :, :, 4);
            clear sw_net

            %% kernel calc dradEffect
            dRclr_sfcKern = zeros(144, 72, ntime, 3);
            dRclr_toaKern = zeros(144, 72, ntime, 3);
            load([dradEffectPath, 'dradEffect_sfc_clr.mat'], 'lw_nonCloudEffect', 'sw_nonCloudEffect', 'nonCloudEffect')
            dRclr_sfcKern(:, :, :, 1) = nonCloudEffect;
            dRclr_sfcKern(:, :, :, 2) = lw_nonCloudEffect;
            dRclr_sfcKern(:, :, :, 3) = sw_nonCloudEffect;
            load([dradEffectPath, 'dradEffect_toa_clr.mat'], 'lw_nonCloudEffect', 'sw_nonCloudEffect', 'nonCloudEffect')
            dRclr_toaKern(:, :, :, 1) = nonCloudEffect;
            dRclr_toaKern(:, :, :, 2) = lw_nonCloudEffect;
            dRclr_toaKern(:, :, :, 3) = sw_nonCloudEffect;

            %% residual
            dRclr_sfcRes = zeros(144, 72, ntime, 3);
            dRclr_toaRes = zeros(144, 72, ntime, 3);
            load([dradEffectPath, 'dR_residual_clr_sfc.mat'])
            dRclr_sfcRes(:, :, :, 1) = dR_residual_clr_sfc;
            dRclr_sfcRes(:, :, :, 2) = lw_residual_clr_sfc;
            dRclr_sfcRes(:, :, :, 3) = sw_residual_clr_sfc;
            load([dradEffectPath, 'dR_residual_clr_toa.mat'])
            dRclr_toaRes(:, :, :, 1) = dR_residual_clr_toa;
            dRclr_toaRes(:, :, :, 2) = lw_residual_clr_toa;
            dRclr_toaRes(:, :, :, 3) = sw_residual_clr_toa;

            % mask
            areaStr = 'world';

            for varNum = 1:3
                [dRclr_sfc(:, :, :, varNum), ~, ~] = maskArea(squeeze(dRclr_sfc(:, :, :, varNum)), lat_f, latRange, -latRange, areaStr);
                [dRclr_sfcKern(:, :, :, varNum), ~, ~] = maskArea(squeeze(dRclr_sfcKern(:, :, :, varNum)), lat_f, latRange, -latRange, areaStr);
                [dRclr_sfcRes(:, :, :, varNum), ~, ~] = maskArea(squeeze(dRclr_sfcRes(:, :, :, varNum)), lat_f, latRange, -latRange, areaStr);
                [dRclr_toa(:, :, :, varNum), ~, ~] = maskArea(squeeze(dRclr_toa(:, :, :, varNum)), lat_f, latRange, -latRange, areaStr);
                [dRclr_toaKern(:, :, :, varNum), ~, ~] = maskArea(squeeze(dRclr_toaKern(:, :, :, varNum)), lat_f, latRange, -latRange, areaStr);
                [dRclr_toaRes(:, :, :, varNum), ~, ~] = maskArea(squeeze(dRclr_toaRes(:, :, :, varNum)), lat_f, latRange, -latRange, areaStr);
            end

            % lat weightMean
            dRclrMean_sfc = zeros(ntime, 3);
            dRclrMean_sfcKern = zeros(ntime, 3);
            dRclrMean_sfcRes = zeros(ntime, 3);
            dRclrMean_toa = zeros(ntime, 3);
            dRclrMean_toaKern = zeros(ntime, 3);
            dRclrMean_toaRes = zeros(ntime, 3);

            for varNum = 1:3

                for timeNum = 1:ntime
                    dRclrMean_sfc(timeNum, varNum) = areaMeanLatWeight(squeeze(squeeze(dRclr_sfc(:, :, timeNum, varNum))), lat_f);
                    dRclrMean_sfcKern(timeNum, varNum) = areaMeanLatWeight(squeeze(squeeze(dRclr_sfcKern(:, :, timeNum, varNum))), lat_f);
                    dRclrMean_sfcRes(timeNum, varNum) = areaMeanLatWeight(squeeze(squeeze(dRclr_sfcRes(:, :, timeNum, varNum))), lat_f);
                    dRclrMean_toa(timeNum, varNum) = areaMeanLatWeight(squeeze(squeeze(dRclr_toa(:, :, timeNum, varNum))), lat_f);
                    dRclrMean_toaKern(timeNum, varNum) = areaMeanLatWeight(squeeze(squeeze(dRclr_toaKern(:, :, timeNum, varNum))), lat_f);
                    dRclrMean_toaRes(timeNum, varNum) = areaMeanLatWeight(squeeze(squeeze(dRclr_toaRes(:, :, timeNum, varNum))), lat_f);
                end

            end

            clear dRclr_*

            % put used vars into one var
            meanarea = zeros(ntime, 2, 3, 3); % (time, sfc/toa, net/long/short, real/kerlCal/Res)
            meanarea(:, 1, :, 1) = dRclrMean_sfc;
            meanarea(:, 1, :, 2) = dRclrMean_sfcKern;
            meanarea(:, 1, :, 3) = dRclrMean_sfcRes;
            meanarea(:, 2, :, 1) = dRclrMean_toa;
            meanarea(:, 2, :, 2) = dRclrMean_toaKern;
            meanarea(:, 2, :, 3) = dRclrMean_toaRes;
            clear dRclr*
            
            
            % cut to 2000.03-2014.02
            timeStr=string(datestr(datenum(time.date),'yyyy-mm'));
            cutStart=find(timeStr=='2000-03');
            cutEnd=find(timeStr=='2014-02');
            time=time.date(cutStart:cutEnd);
            ntime = length(time);
            meanarea=meanarea(cutStart:cutEnd,:,:,:);

            % cal the coef
            cc = zeros(2, 3);

            for sfcToa = 1:2% sfc and toa

                for waveProp = 1:3% total, long,short
                    temp1 = meanarea(:, sfcToa, waveProp, 1);
                    temp2 = meanarea(:, sfcToa, waveProp, 2);
                    cc0 = corrcoef(temp1, temp2, 'Rows', 'complete');
                    cc(sfcToa, waveProp) = cc0(1, 2);
                    cc(sfcToa, waveProp) = roundn(cc(sfcToa, waveProp), -2);
                end

            end

            sfcToalevel = {'sfc', 'toa'}; upperLevel = upper(sfcToalevel);
            Radname = {'dR_{NET} ', 'dR_{LW} ', 'dR_{SW} '};
            set(0, 'defaultfigurecolor', 'w'); %设置画布底色为白色
            set(0, 'DefaultFigureVisible', 'off')

            ss = get(0, 'ScreenSize');
            h = figure('Position', [-1075 188 934 764]); %[离左边缘 离下边缘 自身宽 自身高][ss(4) / 2 - 100 ss(3) / 35 ss(3) / 5 * 2.5, (ss(4) - 80) / 5 * 4]

            for sfcToa = 1:2% sfc and toa

                for waveProp = 1:3% total, short,long
                    subplot_yc(3, 2, waveProp, sfcToa);
                    hold on
                    lineWdth=1.5;
                    plot(time, squeeze(meanarea(:, sfcToa, waveProp, 1)), 'color', 'b','LineWidth',lineWdth) % 深蓝
                    plot(time, squeeze(meanarea(:, sfcToa, waveProp, 3)), 'k','LineWidth',lineWdth)
                    plot(time, squeeze(meanarea(:, sfcToa, waveProp, 2)), 'color', '#C2162E','LineWidth',lineWdth) % 红色
    
    

                    timeSer = [1985 1990 1995 2002 2006 2010 2014 2025 2035 2045 2055 2065 2075 2085 2095 2105];
                    char_timeSer = cellstr(string(timeSer));
                    timeNum = datenum(timeSer, 02, 01);

                    timeSerFull=1985:1:2105;
                    char_timeSerFull = cellstr(string(timeSerFull));
                    timeNumFull = datenum(timeSerFull, 02, 01);
    
    
                    % x 轴坐标
                    datetick('x', 'yyyy', 'keepticks'); % 用日期的格式显示横坐标
                    xlim([time(1) time(ntime)]);

                    xticks(timeNum)
                    xticklabels([])
                    if waveProp==3
                        xticklabels(char_timeSer)
                    end
                    % y 轴坐标
                    ylim([-5 5])

                    if sfcToa==2
                        yticklabels([])
                    end
    
                    ax = gca;
                    ax.XMinorTick = 'on'; ax.YMinorTick = 'on'; % 开启次刻度线
                    ax.XAxis.MinorTickValues=timeNumFull;
                    ax.TickLength = [0.03 0.02]; %刻度线长度      set(gca,'ticklength', [0.02 0.01]);
                    ax.XColor = 'k'; ax.YColor = 'k'; % 设置刻度线颜色
                    ax.FontSize=13;
                    if waveProp == 1
                        text(0.45, 1.1, upperLevel{sfcToa}, 'Fontsize', 15, 'units', 'normalized');
                    end
    
                    if sfcToa == 1
                        ylabel({'Rad Anomaly (Wm^{-2})'}, 'Fontsize', 12.5)
                    end
    
                    legend({Radname{waveProp}, 'Residual', ['Kernel calculated, cc =', num2str(cc(sfcToa, waveProp))]}, 'Fontsize', 8,'FontWeight', 'bold', 'Location', 'northwest','NumColumns', 2); %'Location', 'best'
                    legend('boxoff')%删除图例背景和轮廓
                end
    
            end
            figureTitle = 'Radiation closure experiment (90N-90S mean)';
            sgtitle({figureTitle, ['Model:', level.model2{mdlNum} ', Ensemble: ', esmName{esmNum},', Era: ', level.time1{exmNum}(1:end - 10),' 2000.03-2014.02']})
            figureName = ['Fig2_', level.model2{mdlNum}, '_', esmName{esmNum}, '_200003-201402_world_radClos'];
            saveFileName = [mPath.Output, '/', figureName, '.png'];
            % saveas(gcf, saveFileName)
            save_png(saveFileName)%high resolution
            close gcf
                
            
        end

    end

end

eval(['cd ', nowpath]);
t = toc; disp(t)
