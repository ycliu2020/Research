%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-18 14:06:36
% LastEditTime : 2020-06-24 14:54:38
% LastEditors  : LYC
% Description  : test test
% FilePath     : /code/p2_processCMIP6Data/s99.plot/z_test.m
%  
%%---------------------------------------------------------
inputPath='/data1/liuyincheng/CMIP6-mirror/amip/CMIP/ACCESS-CM2';
File = dir(fullfile(inputPath,'.'));

for i=3:length(File)
    if length(File(i,1).name)~=1&&length(File(i,1).name)~=2
        File1{i-2,1}=fullfile(inputPath,File(i,1).name);
    end
end
