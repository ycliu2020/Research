%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-10-14 15:37:49
% LastEditTime : 2020-10-15 20:48:22
% LastEditors  : LYC
% Description  : 
% FilePath     : /code/p2_processCMIP6Data/s2.radEffTrend/contribAnalysis/contribAnalysis_eastChina.m
%  
%%---------------------------------------------------------
clc; clear;
nowpath = pwd;

%% calculate mask file 
% load mask map
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/中国干湿以及青藏高原分区地图/main_use/mask720.mat')% high resolution
% 国境线
bou_china=shaperead('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/map_data/中国及各省shp/国界线/bou1_4p.shp'); % load china  boundary
bou_chinaX=[bou_china(:).X];bou_chinaY=[bou_china(:).Y];
% 国境线(带南海线)
bou_china_line=shaperead('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/map_data/中国及各省shp/国界线/bou1_4l.shp'); % load china  boundary
bou_china_lineX=[bou_china_line(:).X];bou_china_lineY=[bou_china_line(:).Y];
% 省界线
bou_chinaProvince=shaperead('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/map_data/中国及各省shp/省界线/bou2_4p.shp'); % load china  boundary
bou_chinaProvinceX=[bou_chinaProvince(:).X];bou_chinaProvinceY=[bou_chinaProvince(:).Y];

% maskFile to area of china east 
maskchina_east=maskchina_cp;
maskchina_east(lonw<112,:)=0;
maskchina_east(:,latw>38)=0;

%% loop and process
[mlabels, areaNum] = cmipPlotParameters('atm', 'land', 'radEffect'); % plot parameters

esm = 'r1i1p1f1';
exmStart = 1; exmEnd = 1;

for exmNum = exmStart:exmEnd
    [readme, Experiment, level, tLin, mPlev, vars] = cmipParameters(exmNum);
    % mPath.input:E:/data/cmip6-process/2000-2014/
    exmPath = fullfile('/data1/liuyincheng/cmip6-process/', level.time1{exmNum});
    % mPath.output:a_research/P02.Ts_change_research/figure/04.cmip6Result/2000-2014/ 根据输出的格式进行修改
    % mPath.uniOutput = fullfile('/home/liuyc/Research/P02.Ts_change_research/figure/02.cmip6Result/1.6/', lower(mlabels.level), level.time1{exmNum});
    % mPath.Output = fullfile(mPath.uniOutput);
    % auto_mkdir(mPath.Output)

    % model loop
    mdlStart = 1; mdlEnd = 1;%length(level.model2); % differnt models%length(level.model2)
    for mdlNum = mdlStart:mdlEnd
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
            dradTrendPath = fullfile(esmPath, level.process3{7}); %/data1/liuyincheng/cmip6-process/amip_1980-2014/CESM2/Effect_trend
            danomTrendPath = fullfile(esmPath, level.process3{3}); %/data1/liuyincheng/cmip6-process/amip_1980-2014/CESM2/anomaly_trend


            
        end
        disp([level.model2{level1}, ' model is done!'])
        disp(' ')
    end
    disp([level.time1{p_1}, ' era is done!'])
    disp(' ')
end
eval(['cd ', nowpath]);





%% Part 2 variance contribution analysis
% load origon data 

% cal mean of time series

% cal variance of time series
 

