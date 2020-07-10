%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-09 15:52:00
% LastEditTime : 2020-07-08 19:32:32
% LastEditors  : LYC
% Description  : 
% FilePath     : /code/p1_processObserveData/ERAi_old/s4.radClosure/x_radeffects_cal_world.m
%  
%%---------------------------------------------------------
% c_radeffects_cal.m
% this program is aimed to calculate radiation effects(total,short,long)
% TOA and SFC
% choose china/world as the sample
% chooes clear sky/all sky
clc; clear;tic; %计时


% kind_eff = {'E:\Repeat_the_experiment\testdata\ERAi\eraiRadAnom_2018SFCclrsky_r.nc','E:\Repeat_the_experiment\testdata\ERAi\eraiRadAnom_2018SFCallsky_r.nc'}; % kernel's path
% kind_eff = {'/home/lyc/repeat_Kernel_experiment/testdata/ERAi/eraiRadAnom_2018SFCclrsky_r.nc','/home/lyc/repeat_Kernel_experiment/testdata/ERAi/eraiRadAnom_2018TOAclrsky_r.nc',...
%     '/home/lyc/repeat_Kernel_experiment/testdata/ERAi/eraiRadAnom_2018SFCallsky_r.nc','/home/lyc/repeat_Kernel_experiment/testdata/ERAi/eraiRadAnom_2018TOAallsky_r.nc'}; % clr
kind_eff = {'/home/lyc/repeat_Kernel_experiment/testdata/ERAi/eraiRadAnom_2018SFCclrsky_r.nc','/home/lyc/repeat_Kernel_experiment/testdata/ERAi/eraiRadAnom_2018TOAclrsky_r.nc',...
    '/home/lyc/repeat_Kernel_experiment/testdata/ERAi/eraiRadAnom_2018SFCallsky_r.nc','/home/lyc/repeat_Kernel_experiment/testdata/ERAi/eraiRadAnom_2018TOAallsky_r.nc'}; % clr
f_numb = 1;%1.clr, 3.all
%file_rrad ={'E:\Repeat_the_experiment\testdata\ERAi\ERAi_radAnomclr.nc','E:\Repeat_the_experiment\testdata\ERAi\ERAi_radAnom_all.nc'};
kind_rrad ={'/home/lyc/repeat_Kernel_experiment/testdata/ERAi/ERAi_radAnom_clr.nc','/home/lyc/repeat_Kernel_experiment/testdata/ERAi/ERAi_radAnom_all.nc'} ;%real radiation
fr_numb = 1;%1.clr, 2.all

%load effects
%s_albeff,l_teff,l_tseff,l_wveff,s_wveff,totaleff
%s_sumdff,l_sumdff,totaleff
lone = ncread(kind_eff{f_numb},'longitude');%144x72
late = ncread(kind_eff{f_numb},'latitude');
time = ncread(kind_eff{f_numb},'time');

s_albeff=zeros(144,72,216,4);l_teff=s_albeff; l_tseff=s_albeff;l_wveff=s_albeff;s_wveff=s_albeff;totaleff=s_albeff;
for i=1:4
    s_albeff(:,:,:,i) = ncread(kind_eff{i},'albRadEff');
    l_teff(:,:,:,i) = ncread(kind_eff{i},'tRadEff');
    l_tseff(:,:,:,i) = ncread(kind_eff{i},'tsRadEff');
    l_wveff(:,:,:,i) = ncread(kind_eff{i},'wvlwRadEff');
    s_wveff(:,:,:,i) = ncread(kind_eff{i},'wvswRadEff');
    
    totaleff(:,:,:,i) = ncread(kind_eff{i},'totalRadEff');
end
s_sumdff = s_albeff+s_wveff ;
l_sumdff = l_teff+l_tseff+l_wveff ;
%load real radiation
%l_rad,r_rad,toatalrad
l_rad=zeros(144,72,216,4);s_rad=l_rad;
for i=1:2
    if i==1
        l_rad(:,:,:,1) = ncread(kind_rrad{i},'dR_lw_sfc_clr');
        s_rad(:,:,:,1) = ncread(kind_rrad{i},'dR_sw_sfc_clr');
        l_rad(:,:,:,2) = ncread(kind_rrad{i},'dR_lw_toa_clr');
        s_rad(:,:,:,2) = ncread(kind_rrad{i},'dR_sw_toa_clr');
    else
        l_rad(:,:,:,3) = ncread(kind_rrad{i},'dR_lw_sfc_all');
        s_rad(:,:,:,3) = ncread(kind_rrad{i},'dR_sw_sfc_all');
        l_rad(:,:,:,4) = ncread(kind_rrad{i},'dR_lw_toa_all');
        s_rad(:,:,:,4) = ncread(kind_rrad{i},'dR_sw_toa_all');
    end
end
totalrad = l_rad+s_rad;
%china area mean
meanarea=zeros(216,4,3,2);% 216montime,4toa/sfc,3radiation kinds,2real/effect
[Lone,Late] = ndgrid(lone, late);
wgt = cos(Late./180.*pi);%纬度加权算区域平均值

load('/home/lyc/tools/matlab/map/a_world_map/mask/mask_cp144')% maskworld_cp
mask_t=maskworld_cp;
for i = 1:4
    for j = 1:length(l_teff)
        tmp1=totalrad(:,:,j,i);tmp2=l_rad(:,:,j,i);tmp3=s_rad(:,:,j,i);
        tmp11=totaleff(:,:,j,i);tmp22=l_sumdff(:,:,j,i);tmp33=s_sumdff(:,:,j,i);
        tmp1(~mask_t)=NaN;tmp2(~mask_t)=NaN;tmp3(~mask_t)=NaN;
        tmp11(~mask_t)=NaN;tmp22(~mask_t)=NaN;tmp33(~mask_t)=NaN;
        wgt(~mask_t)=NaN;
        % 只选取-60到60
        for ii = 1:length(late)
            if abs(late(ii))>=60
                tmp11(:,ii)=NaN;tmp22(:,ii)=NaN;tmp33(:,ii)=NaN;
                wgt(:,ii)=NaN;
            end
        end
        meanarea(j,i,1,1)=nansum(nansum(tmp1.*wgt))/nansum(nansum(wgt));
        meanarea(j,i,2,1)=nansum(nansum(tmp2.*wgt))/nansum(nansum(wgt));
        meanarea(j,i,3,1)=nansum(nansum(tmp3.*wgt))/nansum(nansum(wgt));
        meanarea(j,i,1,2)=nansum(nansum(tmp11.*wgt))/nansum(nansum(wgt));
        meanarea(j,i,2,2)=nansum(nansum(tmp22.*wgt))/nansum(nansum(wgt));
        meanarea(j,i,3,2)=nansum(nansum(tmp33.*wgt))/nansum(nansum(wgt));
        % totaleff(:,:,j,i)=tmp1;s_sumdff(:,:,j,i)=tmp2;l_sumdff(:,:,j,i)=tmp3;
        % totalrad(:,:,j,i)=tmp11;s_rad(:,:,j,i)=tmp22;l_rad(:,:,j,i)=tmp33;
    end
end
meanarea(:,:,:,3)=meanarea(:,:,:,1)-meanarea(:,:,:,2);

outfilename='/home/lyc/repeat_Kernel_experiment/testdata/radclosure_date/wordmean_r.mat';
save(outfilename,'lone','late','time','meanarea')





t=toc;
disp(t)