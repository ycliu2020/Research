%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-07-08 19:48:12
% LastEditTime : 2020-07-09 20:04:47
% LastEditors  : LYC
% Description  : plot radClos
% FilePath     : /code/p1_processObserveData/ERAi/checkRadClos/ERAi_plot_RadClos.m
%
%%---------------------------------------------------------
clc; clear;
lon_k = 0:2.5:357.5; nlonk = length(lon_k); % kernel lat lon
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);
areaName = {'china', 'world'};
areaLabel = {'cn', 'un'};
Label.wave = {'net', 'sw', 'lw'};

figureTitle = {'Radiation closure experiment (60N-60S mean)', 'Radiation closure experiment (China mean)'};

[readme, level, tLin, vars] = obsParameters('ERAi');

for p_1 = 1:2
    radEfectPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi', level.standVarPath{5}); %radEffect
    outPutPath = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/01.observe/1.3/', tLin.time{p_1}, 'ERAi', 'fig_radClos');
    auto_mkdir(outPutPath)

    % load
    load([radEfectPath, 'RadMean_china.mat'])%'lon_f', 'lat_f', 'time', 'dRclrMean_sfc_cn', 'dRclrMean_toa_cn', 'dRclrMean_sfcKern_cn', 'dRclrMean_toaKern_cn', 'dRclrMean_sfcRes_cn', 'dRclrMean_toaRes_cn'
    load([radEfectPath, 'RadMean_world.mat'])% 'lon_f', 'lat_f', 'time', 'dRclrMean_sfc_un', 'dRclrMean_toa_un', 'dRclrMean_sfcKern_un', 'dRclrMean_toaKern_un', 'dRclrMean_sfcRes_un', 'dRclrMean_toaRes_un'
    ntime = length(time);

    %% plot
    for areaNum = 1:2% 1 Mean world 2 mean china
        meanarea = zeros(ntime, 2, 3, 3); % (time, sfc/toa, net/long/short, real/kerlCal/Res)

        if areaNum == 1% world
            meanarea(:, 1, :, 1) = dRclrMean_sfc_un;
            meanarea(:, 1, :, 2) = dRclrMean_sfcKern_un;
            meanarea(:, 1, :, 3) = dRclrMean_sfcRes_un;
            meanarea(:, 2, :, 1) = dRclrMean_toa_un;
            meanarea(:, 2, :, 2) = dRclrMean_toaKern_un;
            meanarea(:, 2, :, 3) = dRclrMean_toaRes_un;
        elseif areaNum == 2% china
            meanarea(:, 1, :, 1) = dRclrMean_sfc_cn;
            meanarea(:, 1, :, 2) = dRclrMean_sfcKern_cn;
            meanarea(:, 1, :, 3) = dRclrMean_sfcRes_cn;
            meanarea(:, 2, :, 1) = dRclrMean_toa_cn;
            meanarea(:, 2, :, 2) = dRclrMean_toaKern_cn;
            meanarea(:, 2, :, 3) = dRclrMean_toaRes_cn;

        end

        % cal the coef
        cc = zeros(2, 3);

        for sfcToa = 1:2% sfc and toa

            for waveProp = 1:3% total, long,short
                temp1 = meanarea(:, sfcToa, waveProp, 1);
                temp2 = meanarea(:, sfcToa, waveProp, 2);
                cc0 = corrcoef(temp1, temp2, 'Rows', 'complete');
                cc(sfcToa, waveProp) = cc0(1, 2);
                cc(sfcToa, waveProp) = roundn(cc(sfcToa, waveProp), -3);
            end

        end

        sfcToalevel = {'sfc', 'toa'}; upperLevel = upper(sfcToalevel);
        Radname = {'dR_{NET} ', 'dR_{LW} ', 'dR_{SW} '};
        set(0, 'defaultfigurecolor', 'w')%设置画布底色为白色

        ss = get(0, 'ScreenSize');
        h = figure('Position', [ss(4) / 2 - 100 ss(3) / 35 ss(3) / 5 * 2.5, (ss(4) - 80) / 5 * 4]); %[离左边缘 离下边缘 自身宽 自身高]

        for sfcToa = 1:2% sfc and toa

            for waveProp = 1:3% total, short,long
                subplot_yc(3, 2, waveProp, sfcToa)
                hold on
                plot(time, squeeze(meanarea(:, sfcToa, waveProp, 1)), 'b')
                plot(time, squeeze(meanarea(:, sfcToa, waveProp, 2)), 'r')
                plot(time, squeeze(meanarea(:, sfcToa, waveProp, 3)), 'k')
                %设置坐标轴
                datetick('x', 'yyyy', 'keepticks'); % 用日期的格式显示横坐标
                xlim([time(1) time(ntime)])

                if areaNum == 1
                    ylim([-5 5])
                else
                    ylim([-10 10])
                end

                xticks([732358 734153 735979])
                xticklabels({'2005', '2010', '2015'})
                ax = gca;
                ax.XMinorTick = 'on'; ax.YMinorTick = 'on'; % 开启次刻度线
                ax.TickLength = [0.02 0.01]; %刻度线长度      set(gca,'ticklength', [0.02 0.01]);
                ax.XColor = 'k'; ax.YColor = 'k'; % 设置刻度线颜色

                ylabel({upperLevel{sfcToa}, ' Rad Anomaly (Wm^{-2})'}, 'Fontsize', 10)
                legend({Radname{waveProp}, ['Diagnosed cc =', num2str(cc(sfcToa, waveProp))], 'Residual'}, 'Fontsize', 6, 'Location', 'northwest'); %'Location', 'best'
                legend('boxoff')%删除图例背景和轮廓
            end

        end
        sgtitle({figureTitle{areaNum}, 'Data: ERAi'})
        figureName=['ERAi_', tLin.time{p_1}, '_',areaName{areaNum}, '_radClos'];
        saveFileName = [outPutPath, '/', figureName, '.png'];
        saveas(gcf, saveFileName)
        % save_png(saveFileName)%high resolution
        close gcf

    end

end
