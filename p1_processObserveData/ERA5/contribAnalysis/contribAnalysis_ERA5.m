%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-11-07 15:09:57
% LastEditTime : 2020-11-07 22:01:15
% LastEditors  : LYC
% Description  :
% FilePath     : /code/p1_processObserveData/ERA5/contribAnalysis/contribAnalysis_ERA5.m
% symbol_custom_string_obkoro1:
%%---------------------------------------------------------
clear;
nowpath = pwd;
dbstop if error
%% mask file
% load mask map
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_cp144.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/mask/mask_ce72.mat')% load word land mask
load('/home/liuyc/lib/tools/matlab/plot/myMap/02.world_map/mat_file/correct_worldmap.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/mask14472.mat')
load('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/mat_file/中国干湿以及青藏高原分区地图/main_use/mask720.mat')
% 国境线
bou_china = shaperead('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/map_data/中国及各省shp/国界线/bou1_4p.shp'); % load china  boundary
bou_chinaX = [bou_china(:).X]; bou_chinaY = [bou_china(:).Y];
% 国境线(带南海线)
bou_china_line = shaperead('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/map_data/中国及各省shp/国界线/bou1_4l.shp'); % load china  boundary
bou_china_lineX = [bou_china_line(:).X]; bou_china_lineY = [bou_china_line(:).Y];
% 省界线
bou_chinaProvince = shaperead('/home/liuyc/lib/tools/matlab/plot/myMap/01.china_map/map_data/中国及各省shp/省界线/bou2_4p.shp'); % load china  boundary
bou_chinaProvinceX = [bou_chinaProvince(:).X]; bou_chinaProvinceY = [bou_chinaProvince(:).Y];

% lon and lat
lon_k = 0:2.5:357.5; nlonk = length(lon_k); % kernel lat lon
lat_k = 90:-2.5:-90; nlatk = length(lat_k);
lat_f = 88.75:-2.5:-88.75; nlatf = length(lat_f); % figure lat lon
lon_f = lon_k; nlonf = length(lon_f);

% Path
ERA5Path = '/data1/liuyincheng/Observe-process/200003-201802/ERA5';

% Read rhs(RHeating)
anomPath = fullfile(ERA5Path, 'anomaly');
load(fullfile(anomPath, 'global_vars.mat'))
load(fullfile(anomPath, 'drhs.mat'))% dhFlux, drhs, drhsPlus
% regrid
ntime = length(time);
drhs = autoRegrid3(lon_k, lat_k, time, drhs, lon_f, lat_f, time);

% Read sfc values
dradEffectPath = fullfile(ERA5Path, 'radEffect');
load(fullfile(dradEffectPath, 'global_vars.mat'))
load(fullfile(dradEffectPath, 'dR_cloud_sfc'))
load(fullfile(dradEffectPath, 'dradEfect_sfc_cld.mat'))% husEffect, mainEffect, taEffect, taOnlyEffect, taOnlyEffect2, tasEffect, tasEffect2, totalEffect, tsEffect, wvlwEffect, wvswEffect, albEffect
load(fullfile(dradEffectPath, 'dR_residual_cld_sfc.mat'))% dR_resiual_cld_sfc
load(fullfile(dradEffectPath, 'real_dradEfect.mat'))% l_rad, readme_realradEfect, s_rad, dR_clr, dR_allsky

% Add all vars into one var
dTs_x_sfc = zeros(nlonf, nlatf, ntime, 2);
dTs_x_sfc(:, :, :, 1) = drhs;
dTs_x_sfc(:, :, :, 2) = drhs - (mainEffect + dR_cloud_sfc); % or -squeeze(dR_allsky(:,:,:,1))+dR_residual_cld_sfc;
dTs_x_sfc(:, :, :, 3) = albEffect;
dTs_x_sfc(:, :, :, 4) = dR_cloud_sfc;
dTs_x_sfc(:, :, :, 5) = husEffect;
dTs_x_sfc(:, :, :, 6) = taEffect;
size_dTs_x_sfc = size(dTs_x_sfc);
nValue = size_dTs_x_sfc(4);

% Choose Area
latRange = 80;
areaNum = 'china east';

for x_ind = 1:nValue
    [dTs_x_sfc(:, :, :, x_ind), ~, ~] = maskArea(dTs_x_sfc(:, :, :, x_ind), lat_f, latRange, -latRange, areaNum);
end

% Cal the weight average value

glb_dTs_x_sfc = zeros(ntime, nValue);

for timeNum = 1:ntime

    for x_ind = 1:nValue
        glb_dTs_x_sfc(timeNum, x_ind) = areaMeanLatWeight(dTs_x_sfc(:, :, timeNum, x_ind), lat_f);
    end

end
glbMoth_dTs_x_sfc=glb_dTs_x_sfc;

%% cal variance and covariance
cov_glbMoth_dTs_x_sfc = cov(glbMoth_dTs_x_sfc(:, 2:end));
var_glbMoth_dTs_sfc = var(glbMoth_dTs_x_sfc(:, 1));
cov_glbMoth_dTs_x_sfc = cov_glbMoth_dTs_x_sfc ./ var_glbMoth_dTs_sfc;
% transfor into paper used(lower triangular matrix)
cov_triu1 = triu(cov_glbMoth_dTs_x_sfc, 1) .* 2;
cov_diag = diag(diag(cov_glbMoth_dTs_x_sfc));
cov_glbMoth_dTs_x_sfc = cov_triu1 + cov_diag;
cov_glbMoth_dTs_x_sfc = fliplr(cov_glbMoth_dTs_x_sfc);
disp('Rheating decompose: ')
disp(cov_glbMoth_dTs_x_sfc)
% check sum equal to 1
check_sum = sum(sum(cov_glbMoth_dTs_x_sfc));
disp('check sum: ')
disp(check_sum)
