%kernel_accumulate_lyc.m
% this program is aimed to calculate every level's accumulate quantitative value of kernels(nodp)
% ------------------------previous tips--------------------------------
% in the file dp.nc
% player denotes the middle pressure of each layer we perturbed
% plevel denotes the boundary pressure of each layer
% dp denotes the thickness (in hPa) of each layer
clc; clear; tic;
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

% read lat lon
lon_k = ncread('/data1/liuyincheng/y_kernels/kernels_YiH/toa/RRTMG_wv_lw_toa_cld_highR.nc', 'lon'); nlonk = length(lon_k);
lat_k = ncread('/data1/liuyincheng/y_kernels/kernels_YiH/toa/RRTMG_wv_lw_toa_cld_highR.nc', 'lat'); nlatk = length(lat_k);

%% different time series, 1mean 2000-03 to 2018-02(18*12). 2 mean 200207-201706(15*12)
[readme, level, tLin, vars] = obsParameters('ERAi');
for p_1 = 1:2
    kernelCalPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi', level.standVarPath{4});
    auto_mkdir(kernelCalPath)
    %%%%% Part1 read dps
    load([kernelCalPath, 'kernel_dp.mat']) % 'dps', 'dp', 'dp_level2'
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

        % save([kernelsOutPath,'global_vars.mat'], 'lonf', 'latf', 'time','plevf','readme','timeseries','modelname')
        save([kernelCalPath, saveName{sfcToa}], 'alb_swkernel', 'ts_lwkernel', 't_lwkernel', 't_level2_lwkernel', 'wv_lwkernel', 'wv_swkernel');
    end
    clear dp dps dp_level2
end
