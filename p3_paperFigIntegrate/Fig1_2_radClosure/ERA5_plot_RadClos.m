%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-07-08 19:48:12
% LastEditTime : 2021-04-01 15:09:23
% LastEditors  : Please set LastEditors
% Description  : plot radClos
% FilePath     : /code/p3_paperFigIntegrate/Fig1_2_radClosure/ERA5_plot_RadClos.m
%
%%---------------------------------------------------------
clc; clear;
lon_k = 0:2.5:357.5; nlonk = length(lon_k); % kernel lat lon
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);

Label.wave = {'net', 'sw', 'lw'};
areaLabel = {'un', 'cn'};
areaName = {'world', 'china'};
figureTitle = {'Radiation closure experiment (90N-90S mean)', 'Radiation closure experiment (China mean)'};

[readme, level, tLin, vars] = obsParameters('ERA5');

for exmNum = 3:3 % only 2000.03-2014.02
    radEffectPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{exmNum}, 'ERA5', level.standVarPath{5}); %radEffect
    outPutPath = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/proj3_PaperFig/v0.3/');
    auto_mkdir(outPutPath)

    % load
    load([radEffectPath, 'RadMean_china.mat']) %'lon_f', 'lat_f', 'time', 'dRclrMean_sfc_cn', 'dRclrMean_toa_cn', 'dRclrMean_sfcKern_cn', 'dRclrMean_toaKern_cn', 'dRclrMean_sfcRes_cn', 'dRclrMean_toaRes_cn'
    load([radEffectPath, 'RadMean_world.mat']) % 'lon_f', 'lat_f', 'time', 'dRclrMean_sfc_un', 'dRclrMean_toa_un', 'dRclrMean_sfcKern_un', 'dRclrMean_toaKern_un', 'dRclrMean_sfcRes_un', 'dRclrMean_toaRes_un'
    ntime = length(time);

    %% plot
    for areaNum = 1:2 % 1 Mean world 2 mean china
        meanarea = zeros(ntime, 2, 3, 3); % (time, sfc/toa, net/long/short, real/kerlCal/Res)

        if areaNum == 1 % world
            meanarea(:, 1, :, 1) = dRclrMean_sfc_un;
            meanarea(:, 1, :, 2) = dRclrMean_sfcKern_un;
            meanarea(:, 1, :, 3) = dRclrMean_sfcRes_un;
            meanarea(:, 2, :, 1) = dRclrMean_toa_un;
            meanarea(:, 2, :, 2) = dRclrMean_toaKern_un;
            meanarea(:, 2, :, 3) = dRclrMean_toaRes_un;
        elseif areaNum == 2 % china
            meanarea(:, 1, :, 1) = dRclrMean_sfc_cn;
            meanarea(:, 1, :, 2) = dRclrMean_sfcKern_cn;
            meanarea(:, 1, :, 3) = dRclrMean_sfcRes_cn;
            meanarea(:, 2, :, 1) = dRclrMean_toa_cn;
            meanarea(:, 2, :, 2) = dRclrMean_toaKern_cn;
            meanarea(:, 2, :, 3) = dRclrMean_toaRes_cn;

        end

        % cal the coef
        cc = zeros(2, 3);

        for sfcToa = 1:2 % sfc and toa

            for waveProp = 1:3 % total, long,short
                temp1 = meanarea(:, sfcToa, waveProp, 1);
                temp2 = meanarea(:, sfcToa, waveProp, 2);
                cc0 = corrcoef(temp1, temp2, 'Rows', 'complete');
                cc(sfcToa, waveProp) = cc0(1, 2);
                cc(sfcToa, waveProp) = roundn(cc(sfcToa, waveProp), -2);
            end

        end

        sfcToalevel = {'sfc', 'toa'}; upperLevel = upper(sfcToalevel);
        Radname = {['d', '\itR', '\rm_{net} '], ['d', '\itR', '\rm_{LW} '], ['d', '\itR', '\rm_{SW} ']};
        set(0, 'defaultfigurecolor', 'w'); %设置画布底色为白色

        ss = get(0, 'ScreenSize');
        h = figure('Position', [-1075 188 934 764]); %[离左边缘 离下边缘 自身宽 自身高][ss(4) / 2 - 100 ss(3) / 35 ss(3) / 5 * 2.5, (ss(4) - 80) / 5 * 4]
        subNum=1;
        for sfcToa = 1:2 % sfc and toa

            for waveProp = 1:3 % total, short,long
                subplot_yc(3, 2, waveProp, sfcToa);
                hold on
                lineWdth = 2.5;
                plot(time, squeeze(meanarea(:, sfcToa, waveProp, 3)), 'color', '#808080', 'LineWidth', 1.5)
                plot(time, squeeze(meanarea(:, sfcToa, waveProp, 1)), 'color', '#BC144A', 'LineWidth', lineWdth) %真实值值
                plot(time, squeeze(meanarea(:, sfcToa, waveProp, 2)), 'color', '#30388B', 'LineWidth', lineWdth) % 计算值

                % x 轴时间label
                timeSer = [1985 1990 1995 2002 2006 2010 2014 2020 2025 2035 2045 2055 2065 2075 2085 2095 2105];
                char_timeSer = cellstr(string(timeSer));
                timeNum = datenum(timeSer, 02, 01);

                timeSerFull = 1985:1:2105;
                char_timeSerFull = cellstr(string(timeSerFull));
                timeNumFull = datenum(timeSerFull, 02, 01);

                % x 轴坐标
                datetick('x', 'yyyy', 'keepticks'); % 用日期的格式显示横坐标
                xlim([time(1) time(ntime)]);
                % 标签
                xticks(timeNum)
                xticklabels([])

                if waveProp == 3
                    xticklabels(char_timeSer)
                end

                % y 轴坐标
                if areaNum == 1
                    ylim([-3 4])
                else
                    ylim([-10 10])
                end

                % 标签
                % if sfcToa==2
                %     yticklabels([])
                % end

                if waveProp == 1
                    text(0.45, 1.15, upperLevel{sfcToa}, 'FontName', 'Microsoft YaHei', 'Fontsize', 16, 'units', 'normalized');
                end

                if sfcToa == 1
                    ylabel({'W m^{-2}'}, 'FontName', 'Microsoft YaHei', 'Fontsize', 16)
                end

                gcaGet = get(gca);
                axPosit = gcaGet.Position;
                axPosit = axPosit + [0.01 0.17 0 -0.18]; %(左边缘, 下边缘, 长度, 宽度)
                legend({'Residual', Radname{waveProp}, [Radname{waveProp}, '(kernel calc) cc =', num2str(cc(sfcToa, waveProp))]}, 'Fontsize', 9, 'Location', 'none', 'position', axPosit, 'NumColumns', 2); %'Location', 'best''FontWeight', 'bold',
                legend('boxoff') %删除图例背景和轮廓

                % 添加序号
                if waveProp == 3
                    text(0.475, -0.2, ['(',char(96+subNum),')'], 'FontName', 'Microsoft YaHei', 'FontWeight', 'bold', 'Fontsize', 14, 'units', 'normalized');
                else
                    text(0.475, -0.1, ['(',char(96+subNum),')'], 'FontName', 'Microsoft YaHei', 'FontWeight', 'bold', 'Fontsize', 14, 'units', 'normalized');
                end
                subNum=subNum+1;
                
                ax = gca;
                ax.FontName = 'Microsoft YaHei'; % Microsoft YaHei 'Time New Roman'
                ax.XMinorTick = 'on'; ax.YMinorTick = 'on'; % 开启次刻度线
                ax.XAxis.MinorTickValues = timeNumFull;
                ax.TickLength = [0.03 0.02]; %刻度线长度      set(gca,'ticklength', [0.02 0.01]);
                ax.XColor = 'k'; ax.YColor = 'k'; % 设置刻度线颜色
                ax.FontSize = 12;
                ax.LineWidth = 1.5;

            end

        end

        sgtitle({figureTitle{areaNum}, 'Data: ERA5'}, 'FontName', 'Microsoft YaHei', 'Fontsize', 18);
        figureName = ['Fig1_ERA5_', tLin.time{exmNum}, '_', areaName{areaNum}, '_radClos'];
        saveFileName = [outPutPath, '/', figureName, '.eps'];
        export_fig(gcf, saveFileName, '-r600', '-cmyk')
        % saveas(gcf, saveFileName)
        % save_png(saveFileName)%high resolution
        % close gcf

    end

end
