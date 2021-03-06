% calculation including Ts, Rheating, cloud, Ta, wv, albedo radiation(ERAi date)
% Ts, Ta, wv, albedo can be caled by kernels
% cloud must be caled by the method in xieyan's paper
% Rheating can be caled by raw date about down lw plus net sw (trend already cal)
%
clc; clear; tic;

sfcToaName{1} = 'sfc'; sfcToaName{2} = 'toa';
allClrName{1} = 'all'; allClrName{2} = 'clr';
dataName{1} = 'ERAi'; dataName{2} = 'CERES';
realRad_source = 'CERES';

%% different time series, 1mean 2000-03 to 2018-02(18*12). 2 mean 200207-201706(15*12)
for p_1 = 1:2
    [readme, tLin] = observeParameters(p_1); % fuction% readme,  tLin,
    %modify path first
    inputpath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, '/');
    outpath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi/radEffect_trend/');
    auto_mkdir(outpath)

    %% ------------------------Read date------------------------
    % radiative effect(144,72,tLin.inter{p_1},2,2)sfc/toa*cld/clr
    s_albEff = zeros(144, 72, tLin.inter{p_1}, 2, 2); l_taEff = s_albEff; l_tsEff = s_albEff; l_wvEff = s_albEff; s_wvEff = s_albEff; totalEff = s_albEff;
    l_rad = s_albEff; s_rad = s_albEff;

    for sfcToa = 1:2% 1.sfc, 2.toa
        var1 = char(sfcToaName{sfcToa});

        for allClr = 1:2% 1.all 2.clr
            var2 = char(allClrName{allClr});
            eff_path = strcat(inputpath, dataName{1}, '/radEffect/dradEffect_', var1, '_', var2, '.nc');
            realRad_path = strcat(inputpath, realRad_source, '/anomaly/drad.mat');
            % grobal vars
            lonPlot = ncread(eff_path, 'longitude'); % 144x72
            latPlot = ncread(eff_path, 'latitude');
            nlon = length(lonPlot); nlat = length(latPlot);
            time = ncread(eff_path, 'time'); ntime = length(time);
            % read effects
            s_albEff(:, :, :, sfcToa, allClr) = ncread(eff_path, 'albRadEff');
            l_taEff(:, :, :, sfcToa, allClr) = ncread(eff_path, 'taRadEff');
            l_tsEff(:, :, :, sfcToa, allClr) = ncread(eff_path, 'tsRadEff');
            l_wvEff(:, :, :, sfcToa, allClr) = ncread(eff_path, 'wvlwRadEff');
            s_wvEff(:, :, :, sfcToa, allClr) = ncread(eff_path, 'wvswRadEff');
            totalEff(:, :, :, sfcToa, allClr) = ncread(eff_path, 'totalRadEff');
            % read real rad Anomly
            l_rad(:, :, :, sfcToa, allClr) = cell2mat(struct2cell(load(realRad_path, ['dR_lw_', var1, '_', var2]))); %net thermal radiation
            s_rad(:, :, :, sfcToa, allClr) = cell2mat(struct2cell(load(realRad_path, ['dR_sw_', var1, '_', var2]))); %net thermal radiation
        end

    end

    % cal allsky vars radEffect(only need allsky)
    % vars:  q, alb, ts, ta
    husEff = l_wvEff + s_wvEff;
    dR_rad = l_rad + s_rad;

    dR_husSfcAll = squeeze(squeeze(husEff(:, :, :, 1, 1))); dR_husToaAll = squeeze(squeeze(husEff(:, :, :, 2, 1))); % q
    dR_albSfcAll = squeeze(squeeze(s_albEff(:, :, :, 1, 1))); dR_albToaAll = squeeze(squeeze(s_albEff(:, :, :, 2, 1))); %albedo
    dR_tsSfcAll = squeeze(squeeze(l_tsEff(:, :, :, 1, 1))); dR_tsToaAll = squeeze(squeeze(l_tsEff(:, :, :, 2, 1))); % ts
    dR_taSfcAll = squeeze(squeeze(l_taEff(:, :, :, 1, 1))); dR_taToaAll = squeeze(squeeze(l_taEff(:, :, :, 2, 1))); % ta
    dR_totalSfcAll = squeeze(squeeze(totalEff(:, :, :, 1, 1))); dR_totalToaAll = squeeze(squeeze(totalEff(:, :, :, 2, 1))); % ta
    dR_realRadSfcAll = squeeze(squeeze(dR_rad(:, :, :, 1, 1))); dR_realRadToaAll = squeeze(squeeze(dR_rad(:, :, :, 2, 1))); % ta

    % cal cloud effect(lw+sw)
    dvarsFeedback = squeeze(totalEff(:, :, :, :, 2) - totalEff(:, :, :, :, 1)); % dRfeedback_clr-dRfeedback_all
    dCRF = squeeze(dR_rad(:, :, :, :, 1) - dR_rad(:, :, :, :, 2)); % dR_all-dR_clr
    dR_cloud = dCRF + dvarsFeedback; % cloud feedback
    dR_cloud_sfc = squeeze(dR_cloud(:, :, :, 1)); dR_cloud_toa = squeeze(dR_cloud(:, :, :, 2));

    % save cloud feedback
    % outfilename = 'E:\Repeat_the_experiment\testdata\radcomponent_data\trendradset_erai_allskysfc';
    outfilename = [inputpath,'ERAi/radEffect/', 'dcldRadEffect_', realRad_source, '.mat'];
    save(outfilename, 'lonPlot', 'latPlot', 'time', ...
        'dR_cloud_sfc', 'dR_cloud_toa')

    % dRsfc_total_c mean varsfeedback plus cloud feedback
    dR_totalSfcAll_c = dR_totalSfcAll + dR_cloud_sfc;
    dR_totalToaAll_c = dR_totalSfcAll + dR_cloud_toa;
    %% ------------------------cal trend------------------------
    startmonth = tLin.startmonth{p_1};
    % sfc trend
    [trendm_dRsfc_ta, trends_dRsfc_ta, trendyr_dRsfc_ta, ~, ~] = autoCalTrend(dR_taSfcAll, nlon, nlat, time, startmonth);
    [trendm_dRsfc_ts, trends_dRsfc_ts, trendyr_dRsfc_ts, ~, ~] = autoCalTrend(dR_tsSfcAll, nlon, nlat, time, startmonth);
    [trendm_dRsfc_q, trends_dRsfc_q, trendyr_dRsfc_q, ~, ~] = autoCalTrend(dR_husSfcAll, nlon, nlat, time, startmonth);
    [trendm_dRsfc_alb, trends_dRsfc_alb, trendyr_dRsfc_alb, ~, ~] = autoCalTrend(dR_albSfcAll, nlon, nlat, time, startmonth);
    [trendm_dRsfc_cloud, trends_dRsfc_cloud, trendyr_dRsfc_cloud, ~, ~] = autoCalTrend(dR_cloud_sfc, nlon, nlat, time, startmonth);

    [trendm_dRsfc_total, trends_dRsfc_total, trendyr_dRsfc_total, ~, ~] = autoCalTrend(dR_totalSfcAll, nlon, nlat, time, startmonth);
    [trendm_dRsfc_total_c, trends_dRsfc_total_c, trendyr_dRsfc_total_c, ~, ~] = autoCalTrend(dR_totalSfcAll_c, nlon, nlat, time, startmonth);
    [trendm_dRsfc_realRad, trends_dRsfc_realRad, trendyr_dRsfc_realRad, ~, ~] = autoCalTrend(dR_realRadSfcAll, nlon, nlat, time, startmonth);
    % toa trend
    [trendm_dRtoa_ta, trends_dRtoa_ta, trendyr_dRtoa_ta, ~, ~] = autoCalTrend(dR_taToaAll, nlon, nlat, time, startmonth);
    [trendm_dRtoa_ts, trends_dRtoa_ts, trendyr_dRtoa_ts, ~, ~] = autoCalTrend(dR_tsToaAll, nlon, nlat, time, startmonth);
    [trendm_dRtoa_q, trends_dRtoa_q, trendyr_dRtoa_q, ~, ~] = autoCalTrend(dR_husToaAll, nlon, nlat, time, startmonth);
    [trendm_dRtoa_alb, trends_dRtoa_alb, trendyr_dRtoa_alb, ~, ~] = autoCalTrend(dR_albToaAll, nlon, nlat, time, startmonth);
    [trendm_dRtoa_cloud, trends_dRtoa_cloud, trendyr_dRtoa_cloud, ~, ~] = autoCalTrend(dR_cloud_toa, nlon, nlat, time, startmonth);

    [trendm_dRtoa_total, trends_dRtoa_total, trendyr_dRtoa_total, ~, ~] = autoCalTrend(dR_totalToaAll, nlon, nlat, time, startmonth);
    [trendm_dRtoa_total_c, trends_dRtoa_total_c, trendyr_dRtoa_total_c, ~, ~] = autoCalTrend(dR_totalToaAll_c, nlon, nlat, time, startmonth);
    [trendm_dRtoa_realRad, trends_dRtoa_realRad, trendyr_dRtoa_realRad, ~, ~] = autoCalTrend(dR_realRadToaAll, nlon, nlat, time, startmonth);

    % save as mat file
    outfilename = [outpath, 'trend_dradEfect_all_', realRad_source, '.mat'];
    save(outfilename, 'lonPlot', 'latPlot', 'time', ...
        'trendm_dRsfc_ta', 'trends_dRsfc_ta', 'trendyr_dRsfc_ta', ...
        'trendm_dRsfc_ts', 'trends_dRsfc_ts', 'trendyr_dRsfc_ts', ...
        'trendm_dRsfc_q', 'trends_dRsfc_q', 'trendyr_dRsfc_q', ...
        'trendm_dRsfc_alb', 'trends_dRsfc_alb', 'trendyr_dRsfc_alb', ...
        'trendm_dRsfc_cloud', 'trends_dRsfc_cloud', 'trendyr_dRsfc_cloud', ...
        'trendm_dRsfc_total', 'trends_dRsfc_total', 'trendyr_dRsfc_total', ...
        'trendm_dRsfc_total_c', 'trends_dRsfc_total_c', 'trendyr_dRsfc_total_c', ...
        'trendm_dRsfc_realRad', 'trends_dRsfc_realRad', 'trendyr_dRsfc_realRad', ...
        'trendm_dRtoa_ta', 'trends_dRtoa_ta', 'trendyr_dRtoa_ta', ...
        'trendm_dRtoa_ts', 'trends_dRtoa_ts', 'trendyr_dRtoa_ts', ...
        'trendm_dRtoa_q', 'trends_dRtoa_q', 'trendyr_dRtoa_q', ...
        'trendm_dRtoa_alb', 'trends_dRtoa_alb', 'trendyr_dRtoa_alb', ...
        'trendm_dRtoa_cloud', 'trends_dRtoa_cloud', 'trendyr_dRtoa_cloud', ...
        'trendm_dRtoa_total', 'trends_dRtoa_total', 'trendyr_dRtoa_total', ...
        'trendm_dRtoa_total_c', 'trends_dRtoa_total_c', 'trendyr_dRtoa_total_c', ...
        'trendm_dRtoa_realRad', 'trends_dRtoa_realRad', 'trendyr_dRtoa_realRad')
end

t = toc; disp(t)
