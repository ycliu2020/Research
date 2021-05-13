%%---------------------------------------------------------
% Author       : LYC
% Date         : 2021-04-11 16:37:02
% LastEditors  : Please set LastEditors
% Description  : 
% FilePath     : /code/testExperiment/testWeight.m
%  
%%---------------------------------------------------------

inputPath='/data1/liuyincheng/test/';
lat=ncread([inputPath,'t.kernel.nc'],'lat');
lon=ncread([inputPath,'t.kernel.nc'],'lon');
gw=ncread([inputPath,'t.kernel.nc'],'gw'); %% Gaussian weights for the CESM grid
% Area weight 
weight=repmat(gw(:)',[length(lon) 1]); 
weight2=weight./nansum(nansum(weight));

jiaquan = cosd(lat);
weight1 = ones(length(lon), length(lat));
for latNum = 1:length(lat)
    weight1(:, latNum) = weight1(:, latNum)* jiaquan(latNum); %格点相对大小
end
weight1./nansum(nansum(weight1))
function [varMean] = areaMeanLatWeight(varInput, lat)
    % data must be (lon, lat)
    % varMean: 
    varInput=squeeze(squeeze(varInput));
    jiaquan = cosd(lat);
    sizeVar = size(varInput);
    weight = ones(sizeVar(1), sizeVar(2));
    
    weight(isnan(varInput))=nan;
    for latNum = 1:sizeVar(2)
        weight(:, latNum) = weight(:, latNum)* jiaquan(latNum); %格点相对大小
    end


    varUse = varInput .* weight;
    varMean = nansum(nansum(varUse)) / nansum(nansum(weight));