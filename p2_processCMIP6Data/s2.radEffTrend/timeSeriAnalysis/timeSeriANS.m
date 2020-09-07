%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-08-31 17:00:15
% LastEditTime : 2020-09-01 14:53:30
% LastEditors  : LYC
% Description  : 计算不同变量间的时间序列相关性并画图
% FilePath     : /code/p2_processCMIP6Data/s2.radEffTrend/timeSeriAnalysis/timeSeriANS.m
%
%%---------------------------------------------------------
clc; clear; tic;
nowpath = pwd;
% load mask map
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')

latRange = 88.75; % Latitude range
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
    exmPath = ['/data1/liuyincheng/cmip6-process/', level.time1{exmNum}]; %/data1/liuyincheng/cmip6-process/2000-2014/
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model
    for mdlNum = 1:1%length(level.model2)
        % model path
        mdlPath = fullfile(exmPath, level.model2{mdlNum});
        eval(['cd ', mdlPath]);
        disp(' ')
        disp([level.model2{mdlNum}, ' model start!'])

        % ensemble member path
        esmName = getPath_fileName(mdlPath, '.');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensemble member
        for esmNum = 1:1%length(esmName)
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
            varNames = {'dTs', 'RHeating', 'ta RadEfect', 'ts RadEfect', 'q RadEfect', 'alb RadEfect'};

            % varUsed(:, :, :, 1) = dR_netSfc;
            % varUsed(:, :, :, 2) = dR_cloud_sfc;
            % varUsed(:, :, :, 3) = taEffect;
            % varUsed(:, :, :, 4) = tsEffect;
            % varUsed(:, :, :, 5) = wvlwEffect+wvswEffect;
            % varUsed(:, :, :, 6) = albEffect;
            % varNames = {'Total radiation', 'cloud RadEfect', 'ta RadEfect', 'ts RadEfect', 'q RadEfect', 'alb RadEfect'};
            % transfor to yearly Data
            sizeVarUsed = size(varUsed);
            varUsedYearly = zeros(sizeVarUsed(1), sizeVarUsed(2), sizeVarUsed(3) / 12, sizeVarUsed(4));

            for varNum = 1:sizeVarUsed(4)
                varUsedYearly(:, :, :, varNum) = monthlyToYearly(varUsed(:, :, :, varNum));
            end

            % plot time series and calculate the correlation
            time.vec = datevec(time.date);
            timeYearly = unique(time.vec(:, 1));
            for latNum = 45:45 %nlatf
                for lonNum =  36: 36%nlonNum
                    varUsedTemp=squeeze(squeeze(varUsedYearly(lonNum, latNum, :, :)));
                    plotTimeSeries(varUsedTemp, varNames, timeYearly,toaSfc{2})

                end
                
            end
        end

    end

end
eval(['cd ', nowpath]);
t = toc; disp(t)

function [] = plotTimeSeries(varUsed, varNames, timeYearly, sfcToa)
    %
    % description.
    % varUsed dim is (time,x), x is variable member
    %
    sizeVarUsed = size(varUsed);
    
    % cal cc
    var_cc=cell(1,sizeVarUsed(2));
    var_pp=cell(1,sizeVarUsed(2));
    for varNum = 1 : sizeVarUsed(2)
        [cc0,pp0] = corrcoef(varUsed(:,1),varUsed(:,varNum),'Rows','complete');
        var_cc{varNum}=roundn(cc0(1,2),-2);%保留俩位小数
        var_pp{varNum} = pp0(1,2);% confidence interval
    
    end
    % fig set
    f_row = sizeVarUsed(2); f_col = 1; % 设置画图的行列
    set(0, 'DefaultFigureVisible', 'on')
    ss = get(0, 'ScreenSize');
    h = figure('Position', [ss(4) / 2 - 100 ss(3) / 35 ss(3) / 5 * f_col (ss(4) - 80) / 5 * f_row]);
    % clf reset;
    set(h, 'Color', [1 1 1]);

    for varNum = 1:f_row * f_col
        subplot_yc(f_row, f_col, varNum, f_col);
        hold on
        varUsedTemp=squeeze(varUsed(:, varNum));
        
        plot(timeYearly, varUsedTemp', 'b')
        % 设置坐标轴
        datetick('x', 'yyyy', 'keepticks'); % 用日期的格式显示横坐标
        xlim([timeYearly(1) timeYearly(end)])
        ylim([-10 10])
        %设置显示的x坐标轴标签
        xticks([2005 2010 2015])
        xticklabels({'2005', '2010', '2015'})
        ax = gca;
        ax.XMinorTick = 'on'; ax.YMinorTick = 'on'; % 开启次刻度线
        ax.TickLength = [0.02 0.01]; %刻度线长度      set(gca,'ticklength', [0.02 0.01]);
        ax.XColor = 'k'; ax.YColor = 'k'; % 设置刻度线颜色
        title([varNames{varNum},', &Total cc=',num2str(var_cc{varNum})])

        % ylabel({Level{i}, ' Rad Anomaly (Wm^{-2})'}, 'Fontsize', 10)
        % legend({Radname{j}, ['Diagnosed cc =', num2str(cc(i, j))], 'Residual'}, 'Fontsize', 8, 'Location', 'northwest')
        % legend('boxoff')%删除图例背景和轮廓
    end
    headLineTxt=['Level:' num2str(sfcToa),' Unite:  $W\cdot m^{-2} $'];
    sgtt = sgtitle(headLineTxt, 'Fontsize', 12, 'Interpreter', 'latex');




end
