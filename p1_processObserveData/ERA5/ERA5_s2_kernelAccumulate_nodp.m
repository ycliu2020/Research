%---------------------------------------------------------
% Author       : LYC
% Date         : 2021-04-14 14:03:08
% LastEditors  : Please set LastEditors
% Description  : 
% FilePath     : /code/p1_processObserveData/ERA5/ERA5_s2_kernelAccumulate_nodp.m
%  
%%---------------------------------------------------------
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
kernelsPath = '/data1/liuyincheng/y_kernels/YiH/kernel-nodp-highR/';
toaSfcPath = {'surface', 'toa'};
cldClrLabel = {'cld', 'clr'};
varsKern = {'alb', 'ts', 't', 'wv_lw', 'wv_sw'};
varsKernLabel = {{'alb_swkernel', 'ts_lwkernel', 't_lwkernel', 'wv_lwkernel', 'wv_swkernel'}, ...
    {'swkernel', 'lwkernel', 'lwkernel', 'lwkernel', 'swkernel'}}; %varsKernLabel{2} is rawdata name
sfcToaLabel = {'sfc', 'toa'};
kernelsFileName = cell(1, 4); %filename: 1sfc cld,2sfc clr,3toa cld,4toa clr
% FilePath     : /code/p1_processObserveData/ERA5/ERA5_s2_kernelAccumulate_nodp.m
saveName = cell(1, 4);
n = 0;

for sfcToa = 1:2

    for cldClr = 1:2
        n = n + 1;
        tempName = strcat('RRTMG_', varsKern, '_', toaSfcPath{sfcToa}, '_', cldClrLabel{cldClr}, '_nodp_highR.nc');%RRTMG_alb_surface_cld_nodp_highR.nc
        tempPath = strcat(kernelsPath, toaSfcPath{sfcToa}, '/');
        kernelsFileName{n} = tempName;
        kernelsFilePath{n} = strcat(tempPath, tempName);
        saveName{n} = ['kernels_', sfcToaLabel{sfcToa}, '_', cldClrLabel{cldClr}, '_nodp.mat'];
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

%% different time series, 1mean 2000-03 to 2018-02(18*12). 2 mean 200207-201706(15*12)
[readme, level, tLin, vars] = obsParameters('ERA5');
for p_1 = 1:1
    kernelCalPath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERA5', level.standVarPath{4}); % 'kernelsCal/'
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
        t_level1_lwkernel = t_lwkernel_copy{sfcToa}; % only contain near surface 2 levels
        wv_lwkernel = wv_lwkernel_copy{sfcToa};
        wv_swkernel = wv_swkernel_copy{sfcToa};
        numdp = 24;
        if sfcToa <= 2 % sfc
            % handle and unified unit: W/m2, note variable T is special
            t_level2_lwkernel_temp=t_level2_lwkernel(:,:,2:end,:); %去除最底下一层(ts)
            t_level1_lwkernel_temp=t_level1_lwkernel(:,:,2:end,:); %去除最底下一层(ts)
            temp1=zeros(144,73,24,12);temp2=temp1;
            temp2(dp_level2~=0)=t_level2_lwkernel_temp(dp_level2~=0);
            temp1(dp_level1~=0)=t_level1_lwkernel_temp(dp_level1~=0);
            t_level2_lwkernel(:,:,2:end,:)=temp2;
            t_level1_lwkernel(:,:,2:end,:)=temp1;
        else % toa
            temp1=zeros(144,73,24,12);temp2=temp1;
            temp2(dp_level2~=0)=t_level2_lwkernel(dp_level2~=0);
            temp1(dp_level1~=0)=t_level1_lwkernel(dp_level1~=0);
            t_level2_lwkernel = temp2;
            t_level1_lwkernel = temp1;

        end

        % save([kernelsOutPath,'global_vars.mat'], 'lonf', 'latf', 'time','plevf','readme','timeseries','modelname')
        save([kernelCalPath, saveName{sfcToa}], 'alb_swkernel', 'ts_lwkernel', 't_lwkernel', 'wv_lwkernel', 'wv_swkernel','t_level1_lwkernel','t_level2_lwkernel');
    end
end
