%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-17 13:16:32
% LastEditTime : 2020-07-23 09:30:50
% LastEditors  : LYC
% Description  : 计算非局地云的辐射效应和温度贡献(用云辐射贡献的ts变化求)
% FilePath     : /code/p2_processCMIP6Data/s3.nonLocalCld/test_radClos.m
% Attention    : both amip and ssp on surface caled
%%---------------------------------------------------------

clear; clc; tic;
latRange = 88.75; % Latitude range
lon1 = [2.5 357.5]; lat1 = [-latRange + 1 latRange - 1]; % world area
toaSfc = {'toa', 'sfc'};
lon_k = 0:2.5:357.5; nlonk = length(lon_k); % kernel lat lon
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);

for exmNum = 4:4%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    nowpath = pwd;
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % exmPath
    exmPath = ['/data1/liuyincheng/cmip6-process/', level.time1{exmNum}]; %/data1/liuyincheng/cmip6-process/amip_2000-2014/

    for mdlNum = 5:5%1:length(level.model2);
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
            esmPath = fullfile(mdlPath, esmName{esmNum, 1});
            % data path
            varsPath = fullfile(esmPath, level.process3{1}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/rawdata_regrid
            dvarsPath = fullfile(esmPath, level.process3{2}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
            dvarsTrendPath = fullfile(esmPath, level.process3{3}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly_trend
            kernelPath = fullfile(esmPath, level.process3{5}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/kernelsCal
            dradEffectPath = fullfile(esmPath, level.process3{6}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
            dnonLocalCldPath = fullfile(esmPath, level.process3{8}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/non_localCld/
            auto_mkdir(dnonLocalCldPath)
            load([dradEffectPath, 'global_vars.mat'])% lat_f lon_f time plev_k readme
            ntime = length(time.date);
            % load dtsg and DeltaTsg
            load([dvarsPath, 'dtsg_nomask.mat'])
            load([dvarsTrendPath, 'DeltaTsg_nomask.mat'])% DeltaTsg, DeltaTsgMeaning

            varUsed = zeros(nlonf, nlatf, ntime, 3); % dR
            varKerUsed = varUsed; % dR/kernels
            %% Step1: calculate Heating(J/m2) and DeltaTsg_cld
            month_second = 30 * 24 * 60 * 60; % sum seconds in 1 month
            % cal total heating
            load([varsPath, 'global_vars.mat'])% lat_k, lon_k, plev_k
            load([varsPath, 'rlut.mat'])% toa_outgoing_longwave_flux
            load([varsPath, 'rsut.mat'])% toa_outgoing_shortwave_flux
            load([varsPath, 'rsdt.mat'])% toa_incoming_shortwave_flux
            % unite define a vector component which is positive when directed downward
            netTOA = rsdt - rlut - rsut; % net radiation in toa (lonxlatxmonth): W/m2
            % regrid to 144x72(unite grids)
            netTOA = autoRegrid3(lon_k, lat_k, time.date, netTOA, lon_f, lat_f, time.date);

            Heat_tot1 = sum(netTOA, 3) * month_second;

            % cal total heating (method 2)
            load([dvarsPath, 'global_vars.mat'])% lat_k, lon_k, plev_k
            load([dvarsPath, 'drlut.mat'])% toa_outgoing_longwave_flux
            load([dvarsPath, 'drsut.mat'])% toa_outgoing_shortwave_flux
            load([dvarsPath, 'drsdt.mat'])% toa_incoming_shortwave_flux
            % unite define a vector component which is positive when directed downward
            dnetTOA = drsdt - drlut - drsut; % net radiation in toa (lonxlatxmonth): W/m2
            dnetTOA = autoRegrid3(lon_k, lat_k, time.date, dnetTOA, lon_f, lat_f, time.date);

            Heat_tot2 = sum(dnetTOA, 3) * month_second;

            % cal cloud heating
            load([dradEffectPath, 'dR_cloud_toa.mat'])% dR_cloud_toa
            % control state = The average state in the first five years
            ctrl_Rcloud = mean(dR_cloud_toa(:, :, 1:5 * 12), 3);
            ctrl_Rcloud = repmat(ctrl_Rcloud, [1 1 length(time)]);
            delta_Rcloud = dR_cloud_toa - ctrl_Rcloud;
            Heat_cld = sum(delta_Rcloud, 3) * month_second;

            % global Zonal weighted average (time, vars)
            jiaquan = cosd(lat_f);
            wei = ones(144, 72); %格点纬度加权

            for latiNum = 1:72
                wei(:, latiNum) = wei(:, latiNum) * jiaquan(latiNum); %格点相对大小
            end

            Heat_tot1Average = nansum(nansum(Heat_tot1 .* wei)) / nansum(nansum(wei));
            Heat_tot2Average = nansum(nansum(Heat_tot2 .* wei)) / nansum(nansum(wei));
            Heat_cldAverage = nansum(nansum(Heat_cld .* wei)) / nansum(nansum(wei));
            % cal DeltaTsg_cld
            DeltaTsg_cld = Heat_cldAverage / Heat_tot1Average * DeltaTsg;
            % save
            % save([dradEffectPath, 'DeltaTsg_cld.mat'], 'DeltaTsg_cld', 'Heat_totAverage', 'Heat_cldAverage', 'DeltaTsg')

            % disp([esmName{esmNum, 1}, ' ensemble is done!'])
        end

        disp([level.model2{mdlNum}, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{exmNum}, ' era is done!'])
    disp(' ')

end

t = toc; disp(t)
