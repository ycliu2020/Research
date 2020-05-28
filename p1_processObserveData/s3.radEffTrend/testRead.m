clc; clear; tic;

varLevel1{1} = 'sfc'; varLevel1{2} = 'toa';
varLevel2{1} = 'all'; varLevel2{2} = 'clr';
varLevel3{1} = 'ERAi'; varLevel3{2} = 'CERES';
realRad_source='ERAi';

inputpath = '/data1/liuyincheng/Observe-process/';


% inputpath='E:\Repeat_the_experiment\testdata\CERES_EBAF\';% ��CERES��������Ʒ���ЧӦ
kind_eff = {[inputpath, 'eraiRadAnom_2018SFCallsky_r.nc'], [inputpath, 'eraiRadAnom_2018TOAallsky_r.nc'], ...
    [inputpath, 'eraiRadAnom_2018SFCclrsky_r.nc'], [inputpath, 'eraiRadAnom_2018TOAclrsky_r.nc']}; % radiative effect
kind_rrad = {[inputpath, 'ERAi_radAnom_clr.nc'], [inputpath, 'ERAi_radAnom_all.nc']}; % real radiation
f_numb = 1; % 1.clr_sfc, 2.clr_toa, 3.all_sfc,4.all_toa,

%% Read date
% radiative effect(144,72,216,2,2)sfc/toa*cld/clr
s_albEff = zeros(144, 72, 216, 2, 2); l_taEff = s_albEff; l_tsEff = s_albEff; l_wvEff = s_albEff; s_wvEff = s_albEff; totalEff = s_albEff;
husEff=s_albEff;
for ii = 1:2 % 1.sfc, 2.toa
    var1 = char(varLevel1{ii}); 
    for jj = 1:2 % 1.cld 2.clr
        var2 = char(varLevel2{jj}); 
        eff_path = strcat(inputpath,varLevel3{1} ,'/radEffect/dradEffect_', var1, '_', var2, '.nc');
        realRad_path = strcat(inputpath,realRad_source,'/anomaly/drad.mat');


        % read real rad Anomly
        l_rad(:,:,:,ii,jj)=cell2mat(struct2cell(load(realRad_path,['dR_lw_',var1,'_',var2])));%mat文件中的某变量

    end

end