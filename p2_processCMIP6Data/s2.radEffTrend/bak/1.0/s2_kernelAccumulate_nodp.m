%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% mothly data
% read and process kernels for next step
% this script mainly include:'alb_swkernel','ts_lwkernel','t_lwkernel','wv_lwkernel','wv_swkernel'
% PS: note that sfc t_kernels level is different toa
%
% hope to use an unite scrip to read and process all the varies we need
% Compare the pre rad and ts data from 2000-03 to 2018-02(18 years)
%
% experiment information:
% time:2000.01-2014.12(interval:15*12);1980.01-2014.12(interval:35*12); 2015.01-2099.12(interval:85*12)
% initial time in hist(1740 total): 1,561 of 1740(2000.03);1,321 of 1740(1980.01)
% initial time in futrue(1032 total): 1 of 1032(2015.01);
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear; clc; tic;
nowpath = pwd;
% model settings
for p_1 = 5:5%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(p_1);

    % inputPath
    inputPath = ['/data1/liuyincheng/cmip6-process/', level.time1{p_1}]; %~/data/cmip6/2000-2014/

    % kernel path
    % kernelsPath1 = '/data1/liuyincheng/kernels/kernels_YiH/';
    kernelsPath2='/data1/liuyincheng/y_kernels/kernel-nodp-highR/';

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
            tempfile = strcat('RRTMG_', kernelsLevel3_Vars, '_', kernelsLevel1_hight{sfcToa}, '_', kernelsLevel2_sky{cldClr}, '_nodp_highR.nc');
            temppath = strcat(kernelsPath2, kernelsLevel1_hight{sfcToa}, '/');
            kernelsFileName{n} = tempfile;
            kernelsFilePath{n} = strcat(temppath, tempfile);
            saveName{n} = ['kernels_', kernelsHightLabel{sfcToa}, '_', kernelsLevel2_sky{cldClr}, '_nodp.mat'];
        end

    end

    % read raw kernels
    for ii = 1:4%filename: 1sfc cld,2sfc clr,3toa cld,4toa clr
        % read kernel
        alb_swkernel_temp = ncread(kernelsFilePath{ii}{1}, 'swkernel');
        ts_lwkernel_temp = ncread(kernelsFilePath{ii}{2}, 'lwkernel');
        t_lwkernel_temp = ncread(kernelsFilePath{ii}{3}, 'lwkernel');
        wv_lwkernel_temp = ncread(kernelsFilePath{ii}{4}, 'lwkernel');
        wv_swkernel_temp = ncread(kernelsFilePath{ii}{5}, 'swkernel');
        
        wv_lwkernel_temp(isnan(wv_lwkernel_temp)) = 0; wv_swkernel_temp(isnan(wv_swkernel_temp)) = 0;
        t_lwkernel_temp(isnan(t_lwkernel_temp)) = 0; ts_lwkernel_temp(isnan(ts_lwkernel_temp)) = 0;
        alb_swkernel_temp(isnan(alb_swkernel_temp)) = 0;

        alb_swkernel_copy{ii} = alb_swkernel_temp;
        ts_lwkernel_copy{ii} = ts_lwkernel_temp;
        t_lwkernel_copy{ii} = t_lwkernel_temp;
        wv_swkernel_copy{ii} = wv_swkernel_temp;
        wv_lwkernel_copy{ii} = wv_lwkernel_temp;
    end

    % model loop
    for level1 = 1:length(level.model2)
        % read dps
        inputPsPath = [inputPath, level.model2{level1},'/', level.process3{5}]; %~/data/cmip6/2000-2014/MRI-ESM2-0/kernelsCal
        load([inputPsPath, 'global_vars.mat'])% latf lonf time plevf readme
        load([inputPsPath, 'kernel_dp.mat'])%'dps' and 'dp'and  dp_level2
        nlatf = length(latf); nlonf = length(lonf);
        % process
        for ii = 1:4%filename: 1sfc cld,2sfc clr,3toa cld,4toa clr
            % read kernel
            alb_swkernel = alb_swkernel_copy{ii};
            ts_lwkernel = ts_lwkernel_copy{ii};
            t_lwkernel = t_lwkernel_copy{ii};
            wv_lwkernel = wv_swkernel_copy{ii};
            wv_swkernel = wv_lwkernel_copy{ii};
            numdp = 24;


            kernelsOutPath = [inputPath, level.model2{level1}, '/', level.process3{5}]; % /home/lyc/data/cmip6/2000-2014/MRI-ESM2-0/kernelsCal/
            auto_mkdir(kernelsOutPath)
            % save([kernelsOutPath,'global_vars.mat'], 'lonf', 'latf', 'time','plevf','readme','timeseries','modelname')
            save([kernelsOutPath, saveName{ii}], 'alb_swkernel', 'ts_lwkernel', 't_lwkernel',  'wv_lwkernel', 'wv_swkernel');
        end

        disp([level.model2{level1}, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{p_1}, ' era is done!'])
    disp(' ')
end

t = toc; disp(t)
