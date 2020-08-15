%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-08-09 14:58:33
% LastEditTime : 2020-08-15 10:36:17
% LastEditors  : LYC
% Description  : 修正计算非云致温度变化实验可行性分析 feasibility study report 
% FilePath     : /code/p2_processCMIP6Data/s3.nonLocalCld/s2_nonLocalCld_method3_testEBM_FSR.m
%  
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
    for mdlNum = 1:length(level.model2)
        % model path
        mdlName=level.model2{mdlNum};
        mdlPath = fullfile(exmPath, mdlName);
        eval(['cd ', mdlPath]);
        disp(' ')
        disp([mdlName, ' model start!'])
        % ensemble member path
        esmName = getPath_fileName(mdlPath, '.');
        eval(['cd ', nowpath]);

        %% 暂时只看esm实验
        esm = 'r1i1p1f1';
        if sum(strcmp(esmName,esm))==0
            disp(['the ', esm, ' ensemble of ',mdlName, ' didnt exist']);
            continue
        end
        specificNum=find(strcmp(esmName,'r1i1p1f1')==1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensemble member
        for esmNum = specificNum:specificNum%1:1 %length(esmName) % note that r1i1p1 sometime not the first folder
            esmPath = fullfile(mdlPath, esmName{esmNum, 1});
            % data path
            varsPath = fullfile(esmPath,level.process3{1}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/rawdata_regrid
            dvarsPath = fullfile(esmPath,level.process3{2}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly
            dvarsTrendPath = fullfile(esmPath,level.process3{3}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/anomaly_trend
            kernelPath = fullfile(esmPath,level.process3{5}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/kernelsCal
            dradEffectPath = fullfile(esmPath,level.process3{6}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/radEffect/
            dnonLocalCldPath = fullfile(esmPath,level.process3{8}); %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/non_localCld/

    
            load([dradEffectPath, 'global_vars.mat'])% lat_f lon_f time plev_k readme
            ntime = length(time.date);
            % load dtsg and DeltaTsg
            load([dvarsPath, 'dtsg_nomask.mat']) %dtsg
            load([dvarsTrendPath, 'DeltaTsg_nomask.mat'])% DeltaTsg, DeltaTsgMeaning

            varUsed = zeros(nlonf, nlatf, ntime, 3); % dR
            varKerUsed = varUsed; % dR/kernels
            % load model dvars
            load([dvarsPath, 'dts.mat']) %dts
            load([dradEffectPath, 'dR_cloud_toa.mat']) % cloud radRffect
            load([dradEffectPath, 'dradEfect_toa_cld.mat']) % 'totalEffect','wvlwEffect', 'wvswEffect', 'tsEffect', 'albEffect', 'husEffect',  'taEffect', 'tasEffect2', 'taOnlyEffect2', 'totalEffect', 'mainEffect'
            load([dradEffectPath, 'real_dradEfect.mat']) % 'dR_allsky', 'l_rad', 's_rad',  'dR_clr', 'readme_realradEfect'
            dR_nonCloud=totalEffect;
            dR_netTOA=squeeze(dR_allsky(:,:,:,2));
            dts= autoRegrid3(lon_k, lat_k, time.date, dts, lon_f, lat_f, time.date);
            %% Step1: 将月时间序列转换为年时间序列(lon,lat,1020)→(lon,lat,1020/12)并进行纬向加权平均
            gblM_dts = monthlyToYearlyMean(dts,lat_f);
            gblM_dR_netTOA = monthlyToYearlyMean(dR_netTOA,lat_f);

            gblM_dR_cloud_toa = monthlyToYearlyMean(dR_cloud_toa,lat_f);
            gblM_dR_nonCloud = monthlyToYearlyMean(dR_nonCloud,lat_f);

            gblM_dts_ctrl=mean(gblM_dts(1:10));
            dts_ctrl=gblM_dts-gblM_dts_ctrl;
            RnetTOA_ctrl=mean(gblM_dR_netTOA(1:10));
            dRnetTOA_ctrl=gblM_dR_netTOA-RnetTOA_ctrl;
            %% Step2: 划分时间段以求拟合线
            % % 一层EBM模型(拟合效果不佳)
            % startY=1;endY=80;interY=5;
            % [k1,C1, x1, y1,Integral_ts,Integral_RnetTOA,Delta_ts] = EBMLayer1( gblM_dR_netTOA, gblM_dts, startY, endY, interY, mdlName);
            
            % 二层EBM模型
            % model parameters
            D_ref = 55; Dd_ref = 2768;%;% D_ref = 77; Dd_ref = 1105; in Geoffroy and D_ref = 55; Dd_ref = 2768; in Jiménez-de-la-Cuesta,
            cp_water=4180;% heat capacity J/kg/K
            rho_water=1030; % desnsity of saltwater kg/m3
            f0=0.7; % The proportion of ocean cover the earth surface 
            % cal reference C and C0
            C_ref= (D_ref*rho_water*cp_water*f0)/(86400*365.25);% C_ref(upper-ocean) mean the result of formula (22) in reference Geoffroy
            Cd_ref= (Dd_ref*rho_water*cp_water*f0)/(86400*365.25);% Cd_ref(deep-layer ocean) mean the result of formula (22) in reference Geoffroy

            startY=1;endY=85;interY=5; 
            [Seri_k_dts,Seri_k_dtd,gammaa,b_td] = EBMLayer2Test(C_ref, Cd_ref, gblM_dR_netTOA, gblM_dts, startY, endY, interY, mdlName);
            
            
            % %% Step1: calculate Heating(J/m2) and DeltaTsg_cld
            % month_second = 30 * 24 * 60 * 60; % sum seconds in 1 month
            % year_second = 12*month_second;
            % % cal total heating
            % load([varsPath, 'global_vars.mat'])% lat_k, lon_k, plev_k
            % load([varsPath, 'rlut.mat'])% toa_outgoing_longwave_flux
            % load([varsPath, 'rsut.mat'])% toa_outgoing_shortwave_flux
            % load([varsPath, 'rsdt.mat'])% toa_incoming_shortwave_flux
            % % unite define a vector component which is positive when directed downward
            % netTOA = rsdt - rlut - rsut; % net radiation in toa (lonxlatxmonth): W/m2
            % %regrid to 144x72(unite grids)
            % netTOA = autoRegrid3(lon_k, lat_k, time.date, netTOA, lon_f, lat_f, time.date);

            % Heat_tot = sum(netTOA, 3) * month_second;

            % % cal cloud heating
            % load([dradEffectPath, 'dR_cloud_toa.mat'])% dR_cloud_toa
            % % control state = The average state in the first five years
            % ctrl_Rcloud = mean(dR_cloud_toa(:, :, 1:5 * 12), 3);
            % ctrl_Rcloud = repmat(ctrl_Rcloud, [1 1 length(time)]);
            % delta_Rcloud = dR_cloud_toa - ctrl_Rcloud;
            % Heat_cld = sum(delta_Rcloud, 3) * month_second;

            % % global Zonal weighted average (time, vars)
            % jiaquan = cosd(lat_f);
            % wei = ones(144, 72); %格点纬度加权

            % for latiNum = 1:72
            %     wei(:, latiNum) = wei(:, latiNum) * jiaquan(latiNum); %格点相对大小
            % end

            % Heat_totAverage = nansum(nansum(Heat_tot .* wei)) / nansum(nansum(wei));
            % Heat_cldAverage = nansum(nansum(Heat_cld .* wei)) / nansum(nansum(wei));
            % % cal DeltaTsg_cld
            % DeltaTsg_cld = Heat_cldAverage / Heat_totAverage * DeltaTsg;
            % % save
            % save([dradEffectPath, 'DeltaTsg_cld.mat'], 'DeltaTsg_cld', 'Heat_totAverage', 'Heat_cldAverage', 'DeltaTsg')

            % %% Step2: cal dRnonLocalCld2 and trendyr_dRnonLocalCld2(unite: per day)
            % load([dradEffectPath, 'dradEfect_sfc_cld.mat'])%albEffect, husEffect, taEffect, mainEffect, totalEffect, tsEffect, talwEffect, taswEffect
            % varUsed(:, :, :, 1) = husEffect;
            % varUsed(:, :, :, 2) = taEffect;
            % varUsed(:, :, :, 3) = albEffect;
            % varParitial = zeros(nlonf, nlatf, 3);
            % varKerParitial = varParitial;
            % % Part1: cal trendyr_dRnonLocalCld_toa/sfc(unite: per day)
            % for varNum = 1:3
            %     tempTemp = squeeze(varUsed(:, :, :, varNum));

            %     for latNum = 1:nlatf

            %         for lonNum = 1:nlonf
            %             tempUsed = squeeze(squeeze(tempTemp(lonNum, latNum, :)));
            %             k = polyfit(dtsg, tempUsed, 1); % 一元一次拟合
            %             varParitial(lonNum, latNum, varNum) = k(1);
            %         end

            %     end

            % end

            % dRnonLocalCld2_hus = squeeze(varParitial(:, :, 1)) * DeltaTsg_cld;
            % dRnonLocalCld2_ta = squeeze(varParitial(:, :, 2)) * DeltaTsg_cld;
            % dRnonLocalCld2_alb = squeeze(varParitial(:, :, 3)) * DeltaTsg_cld;
            % dRnonLocalCld2 = dRnonLocalCld2_hus + dRnonLocalCld2_ta + dRnonLocalCld2_alb;

            % trendyr_dRnonLocalCld2 = dRnonLocalCld2 ./ (time.date(end) - time.date(1));
            % save([dnonLocalCldPath, 'dRnonLocalCld2_sfc.mat'], 'dRnonLocalCld2_hus', 'dRnonLocalCld2_ta', 'dRnonLocalCld2_alb', 'dRnonLocalCld2')
            % save([dnonLocalCldPath, 'trendyr_dRnonLocalCld2_sfc.mat'], 'trendyr_dRnonLocalCld2')



            disp([esmName{esmNum, 1}, ' ensemble is done!'])
        end

        disp([mdlName, ' model is done!'])
        disp(' ')
    end

    disp([level.time1{exmNum}, ' era is done!'])
    disp(' ')

end

t = toc; disp(t)



function [Seri_k_dts,Seri_k_dtd,gammaa,b_td] = EBMLayer2Test(C_ref, Cd_ref, dR_netTOA, dts, startT, endT, interT, modelName, figName)
    switch nargin
    case 8
        figName=[];
    end

    timeExamp = startT + interT * 6 - 1:interT:endT;
    diffTimeSum = length(timeExamp); % 总共有几组时间段

    Delta_ts = zeros(diffTimeSum, 1);
    Seri_k_dtd=zeros(diffTimeSum, 1);
    
    Seri_k_dts=zeros(diffTimeSum, 1);
    Seri_b_dts=zeros(diffTimeSum, 1);
    
    Seri_k_dRnetTOA=zeros(diffTimeSum, 1);
    Seri_b_dRnetTOA=zeros(diffTimeSum, 1);
    denominator=zeros(diffTimeSum, 1);
    % 计算基于基态的变量
    ts_ctrl=mean(dts(1:10));
    dts_ctrl=dts-ts_ctrl;
    RnetTOA_ctrl=mean(dR_netTOA(1:10));
    dRnetTOA_ctrl=dR_netTOA-RnetTOA_ctrl;

    for diffTime = 1:diffTimeSum%循环每一组
        timeSeries = startT:timeExamp(diffTime);
        Seri_dRnetTOA_ctrl=dRnetTOA_ctrl(startT:timeExamp(diffTime));
        Seri_dts_ctrl=dts_ctrl(startT:timeExamp(diffTime));
        Seri_dts=dts(startT:timeExamp(diffTime));

        dt = length(timeSeries);

        %% cal Delta ts
        k_dts = polyfit(timeSeries', Seri_dts_ctrl, 1);
        Delta_ts(diffTime) = k_dts(1) * dt;
        Seri_k_dts(diffTime) = k_dts(1);
        Seri_b_dts(diffTime) = k_dts(2);
        k_dRnetTOA = polyfit(timeSeries', Seri_dRnetTOA_ctrl, 1);
        Seri_k_dRnetTOA(diffTime)=k_dRnetTOA(1);
        Seri_b_dRnetTOA(diffTime)=k_dRnetTOA(2);
        denominator(diffTime)=timeExamp(diffTime)+startT;
        Seri_k_dtd(diffTime)=(0.5*k_dRnetTOA(1)*denominator(diffTime)+Seri_b_dRnetTOA(diffTime)-C_ref*Seri_k_dts(diffTime))/Cd_ref;
    end

    y=0.5*(Seri_k_dts-Seri_k_dtd).*denominator+Seri_b_dts;
    x=Seri_k_dtd.*Cd_ref;
    xLabel='Cd*ktd';
    yLabel='0.5*(k\_t-k\_td)(t1+t2)+b\_t';

    
    [k_fit, b_fit] = plotLinFit(x, y, xLabel, yLabel, modelName, figName);
    gammaa=1/k_fit;
    b_td=b_fit;
    title({['ssp370 Model: ', modelName]; ['gamma: ', num2str(gammaa), '   初始deep-ocean温度: ', num2str(b_td)]})
    figName=['/home/liuyc/Research/P02.Ts_change_research/figure/02.cmip6Result/1.6/nonLocalCld2_fix/Cktest/ssp370_r1i1p1f1_',modelName,'_Jim.png'];
    saveas(gcf, figName)

end

function [k, C, x, y,Integral_ts,Integral_RnetTOA,Delta_ts] = EBMLayer1Test(dR_netTOA, dts, startT, endT, interT,modelName,figName)
    %
    % description.
    % dts: d均为全球纬向平均
    % startT: 开始的年
    % endT: 结束的年
    % interT: 中间的间隔
    %
    switch nargin
    case 6
        figName=[];
    end
    timeExamp = startT + interT * 6 - 1:interT:endT;
    diffTimeSum = length(timeExamp); % 总共有几组时间段

    Delta_ts = zeros(diffTimeSum, 1);
    Integral_RnetTOA=zeros(diffTimeSum, 1);
    Integral_ts=zeros(diffTimeSum, 1);
    Seri_k_dts=zeros(diffTimeSum, 1);
    Seri_b_dts=zeros(diffTimeSum, 1);
    Seri_k_dRnetTOA=zeros(diffTimeSum, 1);
    Seri_b_dRnetTOA=zeros(diffTimeSum, 1);
    numerator=zeros(diffTimeSum, 1);
    denominator=zeros(diffTimeSum, 1);
    % 计算基于基态的变量
    ts_ctrl=mean(dts(1:10));
    dts_ctrl=dts-ts_ctrl;
    RnetTOA_ctrl=mean(dR_netTOA(1:10));
    dRnetTOA_ctrl=dR_netTOA-RnetTOA_ctrl;

    for diffTime = 1:diffTimeSum%循环每一组
        timeSeries = startT:timeExamp(diffTime);
        Seri_dRnetTOA_ctrl=dRnetTOA_ctrl(startT:timeExamp(diffTime));
        Seri_dts_ctrl=dts_ctrl(startT:timeExamp(diffTime));
        Seri_dts=dts(startT:timeExamp(diffTime));

        dt = length(timeSeries);

        %% cal Delta ts
        k_dts = polyfit(timeSeries', Seri_dts_ctrl, 1);
        Delta_ts(diffTime) = k_dts(1) * dt;
        Seri_k_dts(diffTime) = k_dts(1);
        Seri_b_dts(diffTime) = k_dts(2);
        k_dRnetTOA = polyfit(timeSeries', Seri_dRnetTOA_ctrl, 1);
        Seri_k_dRnetTOA(diffTime)=k_dRnetTOA(1);
        Seri_b_dRnetTOA(diffTime)=k_dRnetTOA(2);
        %% 计算积分
        % old method
        Integral_RnetTOA(diffTime)=sum(Seri_dRnetTOA_ctrl);
        Integral_ts(diffTime)=sum(Seri_dts_ctrl);
        Integral_ts1(diffTime)=sum(Seri_dts);
        % new method
        denominator(diffTime)=timeExamp(diffTime)+startT;
    end

    % polyfit y and x
    %  plot function Newmethod
    y = (Seri_k_dRnetTOA.*denominator+2*Seri_b_dRnetTOA)./(Seri_k_dts.*denominator+2*Seri_b_dts);
    x = 2*Seri_k_dts ./(Seri_k_dts.*denominator+2*Seri_b_dts);
    xLabel='2k\_ts/(k\_ts(t1+t2)+2b\_ts)';
    yLabel='(k\_netTOA(t1+t2)+2b\_netTOA)/(k\_ts(t1+t2)+2b\_ts)';
    [k_fit, b_fit] = plotLinFit(x, y, xLabel, yLabel, modelName, figName);
    C=k_fit;
    k=b_fit;

    % y = Integral_RnetTOA ./ Integral_ts;
    % x = Delta_ts ./ Integral_ts;
    % kb = polyfit(x, y, 1);
    % C = kb(1);
    % k = kb(2);
    % yfit = polyval(kb, x);
    % figure;
    % plot(x, y, 'o', x, yfit, '-')
    % xlabel('DeltaT/T(t)dt积分')
    % ylabel('RnetTOA(t)dt积分/T(t)dt积分')
    % title({['ssp370 Model: ',modelName];['C= ', num2str(C), ' k= ', num2str(k)]})
    % % 计算r2并标注
    % yresid = y - yfit; %将残差值计算为有符号数的向量：
    % SSresid = sum(yresid.^2); %计算残差的平方并相加，以获得残差平方和：
    % SStotal = (length(y) - 1) * var(y); %通过将观测次数减 1(自由度) 再乘以 y 的方差，计算 y 的总平方和：
    % rsq = 1 - SSresid / SStotal; %计算r2
    % text(0.75, 0.8, ['r2= ', num2str(rsq)], 'units', 'normalized');
    % hold on
    % legend('data', 'linear fit')
    % figName=['/home/liuyc/Research/P02.Ts_change_research/figure/02.cmip6Result/1.6/nonLocalCld2_fix/Cktest/ssp370_r1i1p1f1_',modelName,'.png'];
    % saveas(gcf, figName)
    % % close gcf

    % % teacher test
    % x = Integral_RnetTOA;
    % y = Delta_ts;

    % kb = polyfit(x, y, 1);
    % C = kb(1);
    % k = kb(2);
    % yfit = polyval(kb, x);
    % figure;
    % plot(x, y, 'o', x, yfit, '-')
    % xlabel('netTOA*dt积分')
    % ylabel('\DeltaT(相对于基态的变化)')
    % title({['ssp370 Model: ',modelName];['斜率= ', num2str(C), ' 截距= ', num2str(k)]})
    % % title({['ssp370 Model: ',modelName];['C= ', num2str(C), ' k= ', num2str(k)]})
    % % 计算r2并标注
    % yresid = y - yfit; %将残差值计算为有符号数的向量：
    % SSresid = sum(yresid.^2); %计算残差的平方并相加，以获得残差平方和：
    % SStotal = (length(y) - 1) * var(y); %通过将观测次数减 1(自由度) 再乘以 y 的方差，计算 y 的总平方和：
    % rsq = 1 - SSresid / SStotal; %计算r2
    % text(0.75, 0.8, ['r2= ', num2str(rsq)], 'units', 'normalized');
    % hold on
    % legend('data', 'linear fit')
    % figName=['/home/liuyc/Research/P02.Ts_change_research/figure/02.cmip6Result/1.6/nonLocalCld2_fix/RTtest/ssp370_r1i1p1f1_',modelName,'.png'];
    % % saveas(gcf, figName)
    % % close gcf

end
