%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-09 15:52:00
% LastEditTime : 2020-07-11 10:10:49
% LastEditors  : LYC
% Description  : cal month mean ps and thickness dps, dp
%                this script mainly include:ps and dps, dp
%                PS:
%                ps_m(144x73x12) denotes the month mean ps
%                dps denotes the thickness (in hPa) of surface layer(month mean): dps((144x73x12))
%                dp denotes the thickness (in hPa) of each layer(month mean): dp((144x73x12))
%                experiment information:
%                time:2000.01-2014.12(interval:15*12);1980.01-2014.12(interval:35*12); 2015.01-2099.12(interval:85*12)
%                initial time in hist(1740 total): 1,561 of 1740(2000.03);1,321 of 1740(1980.01)
%                initial time in futrue(1032 total): 1 of 1032(2015.01);
% FilePath     : /code/p2_processCMIP6Data/s2.radEffTrend/s1_kernelAccumulate.m
%
%%---------------------------------------------------------
clear; clc; tic;
nowpath = pwd;
lon_k = 0:2.5:357.5; nlonk = length(lon_k);
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatk = length(lat_f); % figure lat lon
lon_f = lon_k; nlonk = length(lon_f);
%% kernel data
% kernel path
kernelsPath = '/data1/liuyincheng/y_kernels/kernels_YiH/';
kernelsLevel1_hight = {'surface', 'toa'};
kernelsLevel2_sky = {'cld', 'clr'};
kernelsLevel3_Vars = {'alb', 'ts', 't', 'wv_lw', 'wv_sw'};
kernelsVarsLabel = {{'alb_swkernel', 'ts_lwkernel', 't_lwkernel', 'wv_lwkernel', 'wv_swkernel'}, ...
    {'swkernel', 'lwkernel', 'lwkernel', 'lwkernel', 'swkernel'}}; %kernelsVarsLabel{2} is rawdata name
kernelsHightLabel = {'sfc', 'toa'};
kernelsFileName = cell(1, 4); %filename: 1sfc cld,2sfc clr,3toa cld,4toa clr
kernelsFilePath = cell(1, 4); %file path:~/data/kernels/kernels_YiH/surface/RRTMG_alb_surface_clr_highR.nc
saveName = cell(1, 4);
n = 0;

for sfcToa = 1:2

    for cldClr = 1:2
        n = n + 1;
        tempfile = strcat('RRTMG_', kernelsLevel3_Vars, '_', kernelsLevel1_hight{sfcToa}, '_', kernelsLevel2_sky{cldClr}, '_highR.nc');
        temppath = strcat(kernelsPath, kernelsLevel1_hight{sfcToa}, '/');
        kernelsFileName{n} = tempfile;
        kernelsFilePath{n} = strcat(temppath, tempfile);
        saveName{n} = ['kernels_', kernelsHightLabel{sfcToa}, '_', kernelsLevel2_sky{cldClr}, '.mat'];
    end

end

% read raw kernels
for sfcToa = 1:4%filename: 1sfc cld,2sfc clr,3toa cld,4toa clr
    % read kernel
    alb_swkernel_temp = ncread(kernelsFilePath{sfcToa}{1}, 'swkernel');
    ts_lwkernel_temp = ncread(kernelsFilePath{sfcToa}{2}, 'lwkernel');
    t_lwkernel_temp = ncread(kernelsFilePath{sfcToa}{3}, 'lwkernel');
    wv_lwkernel_temp = ncread(kernelsFilePath{sfcToa}{4}, 'lwkernel');
    wv_swkernel_temp = ncread(kernelsFilePath{sfcToa}{5}, 'swkernel');

    wv_lwkernel_temp(isnan(wv_lwkernel_temp)) = 0; wv_swkernel_temp(isnan(wv_swkernel_temp)) = 0;
    t_lwkernel_temp(isnan(t_lwkernel_temp)) = 0; ts_lwkernel_temp(isnan(ts_lwkernel_temp)) = 0;
    alb_swkernel_temp(isnan(alb_swkernel_temp)) = 0;

    alb_swkernel_copy{sfcToa} = alb_swkernel_temp;
    ts_lwkernel_copy{sfcToa} = ts_lwkernel_temp;
    t_lwkernel_copy{sfcToa} = t_lwkernel_temp;
    wv_swkernel_copy{sfcToa} = wv_swkernel_temp;
    wv_lwkernel_copy{sfcToa} = wv_lwkernel_temp;
end

% read kernels height level
exampKernelPath = '/data1/liuyincheng/y_kernels/kernels_YiH/surface/';
plevel = ncread([exampKernelPath, 'dp.nc'], 'plevel');
player = ncread([exampKernelPath, 'dp.nc'], 'player');
dp_raw = ncread([exampKernelPath, 'dp.nc'], 'dp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment
for p_1 = 1:4%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    % model parameters
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(p_1);
    % exmPath
    exmPath = ['/data1/liuyincheng/cmip6-process/', level.time1{p_1}]; %~/data/cmip6/2000-2014/
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model
    for level1 = 1:length(level.model2)% model numbers
        % model path
        mdlPath = fullfile(exmPath, level.model2{level1});
        eval(['cd ', mdlPath]);
        disp(' ')
        disp([level.model2{level1}, ' model start!'])
        % ensemble member path
        esmName = getPath_fileName(mdlPath, '.');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensemble member
        for esmNum = 1:length(esmName)
            % input and output Path
            esmPath = fullfile(mdlPath, esmName{esmNum, 1});
            rawRegridPath = fullfile(esmPath, level.process3{1}); %~/data/cmip6/2000-2014/MRI-ESM2-0/rawdata_regrids
            kernelsOutPath = fullfile(esmPath, level.process3{5}); % /home/lyc/data/cmip6/2000-2014/MRI-ESM2-0/kernelsCal/
            auto_mkdir(kernelsOutPath)
            %%%%% Part1 month mean ps (ps_m(144x73x12))
            load([rawRegridPath, 'global_vars.mat'])% lat_k lon_k time plev_k readme
            load([rawRegridPath, 'ps.mat'])
            ntime = length(time.date); nlatk = length(lat_k); nlonk = length(lon_k);
            ps = ps ./ 100; % transform hpa unite
            ps_m = zeros(nlonk, nlatk, 12);
            [yy, mm, dd] = datevec(time.date);

            for im = 1:12
                index_mm = mm == im;
                ps_m(:, :, im) = nanmean(ps(:, :, index_mm), 3); %this is the average of each month
            end

            month = 1:12;
            % Part2 cal dps,dp
            dp = zeros(nlonk, nlatk, 24, 12); % revise dp in all layers, under the lowest layer are zeros
            dp_level2 = zeros(nlonk, nlatk, 24, 12); % only contain the near sfc
            dps = zeros(nlonk, nlatk, 12);

            for i = 1:nlonk
                for j = 1:nlatk
                    for nt = 1:12
                        if isnan(ps_m(i, j, nt)) == 1
                            continue
                        end
                        temp = find(player < ps_m(i, j, nt), 1, 'first');
                        dps(i, j, nt) = ps_m(i, j, nt) - plevel(temp + 1);
                        % dps(i, j, nt) = dp_bottom(i,j)123;
                        dp(i, j, temp, nt) = dps(i, j, nt);
                        dp(i, j, temp + 1:24, nt) = dp_raw(temp + 1:24);
                        % near surface level2
                        dp_level2(i, j, temp, nt) = dps(i, j, nt);
                        dp_level2(i, j, temp + 1, nt) = dp_raw(temp + 1);
                    end
                end
            end

            save([kernelsOutPath, 'global_vars.mat'], 'lon_k', 'lat_k', 'time', 'plev_k', 'readme', 'timeseries', 'modelname')
            save([kernelsOutPath, 'kernel_dp.mat'], 'dps', 'dp', 'dp_level2');
            save([kernelsOutPath, 'kernel_ps_m.mat'], 'ps_m');
            disp('ps is done!')
            %%%%% Part2 kernel accumlate
            for sfcToa = 1:4%filename: 1sfc cld,2sfc clr,3toa cld,4toa clr
                % read kernel
                alb_swkernel = alb_swkernel_copy{sfcToa};
                ts_lwkernel = ts_lwkernel_copy{sfcToa};
                t_lwkernel = t_lwkernel_copy{sfcToa};
                t_level2_lwkernel = t_lwkernel_copy{sfcToa}; % only contain near surface 2 levels
                wv_lwkernel = wv_lwkernel_copy{sfcToa};
                wv_swkernel = wv_swkernel_copy{sfcToa};
                numdp = 24;

                if sfcToa <= 2%sfc
                    % handle and unified unit: W/m2, note variable T is special
                    t_lwkernel(:, :, 1, :) = squeeze(t_lwkernel(:, :, 1, :)) .* dps / 100;
                    t_level2_lwkernel(:, :, 1, :) = squeeze(t_level2_lwkernel(:, :, 1, :)) .* dps / 100;

                    for i = 1:numdp
                        wv_lwkernel(:, :, i, :) = wv_lwkernel(:, :, i, :) .* dp(:, :, i, :) / 100;
                        wv_swkernel(:, :, i, :) = wv_swkernel(:, :, i, :) .* dp(:, :, i, :) / 100;
                        t_lwkernel(:, :, i + 1, :) = t_lwkernel(:, :, i + 1, :) .* dp(:, :, i, :) / 100;
                        t_level2_lwkernel(:, :, i + 1, :) = t_level2_lwkernel(:, :, i + 1, :) .* dp_level2(:, :, i, :) / 100;
                    end

                else

                    for i = 1:numdp
                        wv_lwkernel(:, :, i, :) = wv_lwkernel(:, :, i, :) .* dp(:, :, i, :) / 100;
                        wv_swkernel(:, :, i, :) = wv_swkernel(:, :, i, :) .* dp(:, :, i, :) / 100;
                        t_lwkernel(:, :, i, :) = t_lwkernel(:, :, i, :) .* dp(:, :, i, :) / 100;
                        t_level2_lwkernel(:, :, i, :) = t_level2_lwkernel(:, :, i, :) .* dp_level2(:, :, i, :) / 100;
                    end

                end
                % save([kernelsOutPath,'global_vars.mat'], 'lon_k', 'lat_k', 'time','plev_k','readme','timeseries','modelname')
                save([kernelsOutPath, saveName{sfcToa}], 'alb_swkernel', 'ts_lwkernel', 't_lwkernel', 't_level2_lwkernel', 'wv_lwkernel', 'wv_swkernel');
            end
            disp([esmName{esmNum, 1}, ' ensemble is done!'])
        end
        disp([level.model2{level1}, ' model is done!'])
        disp(' ')
    end
    disp([level.time1{p_1}, ' era is done!'])
    disp(' ')
end

t = toc; disp(t)
