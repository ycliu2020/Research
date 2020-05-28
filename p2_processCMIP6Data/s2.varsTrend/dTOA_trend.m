%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% mothly data
% cal the anomly trend about any vars in CMIP6 as showed on script name
% this script mainly include:trendm_dts,trends_dts,trendyr_dts; trendm_drhs,trends_drhs,trendyr_drhs...
% PS: m mean month, s mean season, yr mean year
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
% note that TOA use net rad flux
clear; clc; tic;
nowpath = pwd;

for p_1 = 5:5 %1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
% model parameters
[~,modlist_Experiment,level,tLin, mPlev, vars] = modelParameters(p_1);
% chose right month, very important !!! 
startmonth=1;
% input path
inputPath = ['/data1/liuyincheng/cmip6-process/', level.time1{p_1}];%~/data/cmip6/2000-2014/

for level1 =  1:length(level.model2)% model numbers
    inputPath1 = [inputPath, level.model2{level1}, '/', level.process3{2}]; %~/data/cmip6/2000-2014/MRI-ESM2-0/anomaly
    load([inputPath1, 'global_vars.mat'])
    load([inputPath1, 'drlut.mat'])% toa_outgoing_longwave_flux
    load([inputPath1, 'drsut.mat'])% toa_outgoing_shortwave_flux
    load([inputPath1, 'drsdt.mat'])% toa_incoming_shortwave_flux
    % unite define a vector component which is positive when directed downward
    dnetTOA=drsdt-drlut-drsut; % net radiation in toa
    %regrid 144x72(unite grids)
    lat = 88.75:-2.5:-88.75; nlat = length(lat);
    lon=lonf;nlon = length(lon);
    nlonf=length(lonf);nlatf=length(latf);
    [Xlon, Ylat, Ttime] = meshgrid(lat, lon, time);
    [Xlonf, Ylatf, Ttimef] = meshgrid(latf, lonf, time);
    dnetTOA = interp3(Xlonf, Ylatf, Ttimef, dnetTOA, Xlon, Ylat, Ttime);
    % cal the trend
    [trendm_dnetTOA, trends_dnetTOA,trendyr_dnetTOA, p_dnetTOA, cons_dnetTOA] = autoCalTrend(dnetTOA, nlon, nlat, time, startmonth);
    % now we done all the job, now save and output.
    outpathname = [inputPath ,level.model2{level1}, '/', level.process3{3}];%/home/lyc/data/cmip6/2000-2014/MIROC6/anomaly_trend
    auto_mkdir(outpathname)
    save([outpathname,'trend_dnetTOA.mat'],'trendm_dnetTOA',  'cons_dnetTOA',  'p_dnetTOA', 'trends_dnetTOA', 'trendyr_dnetTOA');
    disp([level.model2{level1},' model is done!'])
    disp(' ')
end
disp([level.time1{p_1},' era is done!'])
disp(' ')
end
t = toc; disp(t)
