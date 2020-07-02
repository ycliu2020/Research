% calculation including Ts, Rheating, cloud, Ta, wv, albedo radiation(ERAi date)
% Ts, Ta, wv, albedo can be caled by kernels
% cloud must be caled by the method in xieyan's paper
% Rheating can be caled by raw date about down lw plus net sw (trend already cal)
% 
clc;clear;tic;

input_path='E:\Repeat_the_experiment\testdata\ERAi\rad_effect\';%用ERA辐射计算云辐射效应
% input_path='E:\Repeat_the_experiment\testdata\CERES_EBAF\';% 用CERES辐射计算云辐射效应
kind_eff = {[input_path,'eraiRadAnom_2018SFCallsky_r.nc'],[input_path,'eraiRadAnom_2018TOAallsky_r.nc'],...
[input_path,'eraiRadAnom_2018SFCclrsky_r.nc'],[input_path,'eraiRadAnom_2018TOAclrsky_r.nc']}; % radiative effect
kind_rrad ={[input_path,'ERAi_radAnom_clr.nc'],[input_path,'ERAi_radAnom_all.nc']} ;% real radiation
f_numb = 1; % 1.clr_sfc, 2.clr_toa, 3.all_sfc,4.all_toa,

% --------------------------read effects-------------------------
lone = ncread(kind_eff{f_numb},'longitude'); % 144x72
late = ncread(kind_eff{f_numb},'latitude');
nlon = length(lone);nlat = length(late);
time = ncread(kind_eff{f_numb},'time');ntime = length(time);
% radiative effect(144,72,216,2)sfc,toa
s_albEffAll=zeros(144,72,216,2);l_taEffAll=s_albEffAll; l_tsEffAll=s_albEffAll;l_wvEffAll=s_albEffAll;s_wvEffAll=s_albEffAll;totalEffAll=s_albEffAll;
s_albEffClr=zeros(144,72,216,2);l_taEffClr=s_albEffClr; l_tsEffClr=s_albEffClr;l_wvEffClr=s_albEffClr;s_wvEffClr=s_albEffClr;totalEffClr=s_albEffClr;
for i=1:4
    if i<=2
        s_albEffAll(:,:,:,i) = ncread(kind_eff{i},'albRadEff');
        l_taEffAll(:,:,:,i) = ncread(kind_eff{i},'tRadEff');
        l_tsEffAll(:,:,:,i) = ncread(kind_eff{i},'tsRadEff');
        l_wvEffAll(:,:,:,i) = ncread(kind_eff{i},'wvlwRadEff');
        s_wvEffAll(:,:,:,i) = ncread(kind_eff{i},'wvswRadEff');
        totalEffAll(:,:,:,i) = ncread(kind_eff{i},'totalRadEff');
    else
        s_albEffClr(:,:,:,i-2) = ncread(kind_eff{i},'albRadEff');
        l_taEffClr(:,:,:,i-2) = ncread(kind_eff{i},'tRadEff');
        l_tsEffClr(:,:,:,i-2) = ncread(kind_eff{i},'tsRadEff');
        l_wvEffClr(:,:,:,i-2) = ncread(kind_eff{i},'wvlwRadEff');
        s_wvEffClr(:,:,:,i-2) = ncread(kind_eff{i},'wvswRadEff');
        totalEffClr(:,:,:,i-2) = ncread(kind_eff{i},'totalRadEff');   
    end
end
% note: only cal the allsky condition, and sfc only used 
% q, alb, ts, t
dR_husEffAll = l_wvEffAll+s_wvEffAll;
dR_husSfcAll = dR_husEffAll(:,:,:,1);dR_qToaAll = dR_husEffAll(:,:,:,2); % q
dR_albSfcAll = s_albEffAll(:,:,:,1);dR_albToaAll =s_albEffAll(:,:,:,2);%albedo
dR_tsSfcAll = l_tsEffAll(:,:,:,1);dR_tsToaAll = l_tsEffAll(:,:,:,2);% ts
dR_taSfcAll = l_taEffAll(:,:,:,1);dR_taToaAll = l_taEffAll(:,:,:,2);% t 

% real radiation anomaly
l_rad=zeros(144,72,216,4);s_rad=l_rad;

for i=1:2
    if i==1
        l_rad(:,:,:,1) = ncread(kind_rrad{i},'dR_lw_sfc_clr');% Surface net thermal radiation, 
        l_rad(:,:,:,2) = ncread(kind_rrad{i},'dR_lw_toa_clr');% Top net thermal radiation
        s_rad(:,:,:,1) = ncread(kind_rrad{i},'dR_sw_sfc_clr');% Surface net solar radiation, 
        s_rad(:,:,:,2) = ncread(kind_rrad{i},'dR_sw_toa_clr');% Top net solar radiation,
    else
        l_rad(:,:,:,3) = ncread(kind_rrad{i},'dR_lw_sfc_all');
        l_rad(:,:,:,4) = ncread(kind_rrad{i},'dR_lw_toa_all');
        s_rad(:,:,:,3) = ncread(kind_rrad{i},'dR_sw_sfc_all');
        s_rad(:,:,:,4) = ncread(kind_rrad{i},'dR_sw_toa_all');
    end
end
dR_allsky = l_rad(:,:,:,[3 4])+s_rad(:,:,:,[3 4]);
dR_clr = l_rad(:,:,:,[1 2])+s_rad(:,:,:,[1 2]);

% ------------------------cal the cloud effect(lw+sw)---------------------------
% dR_cloud(144,72,216,2)(最后一维代表sfc和toa)
% for research only clear and all sky's surface needed (1 and 3 are used)

%% cal the dR_cloud (144,72,216,2)包含地表和大气层顶, 所以最后一维为2, 和地表温度
dR_res = dR_clr-totalEffClr;
dR_cloud = dR_allsky - totalEffAll - dR_res;
dR_cloud_sfc = dR_cloud(:,:,:,1);dR_cloud_toa = dR_cloud(:,:,:,2);%不区分有云无云

% ------------------------trend of variable(原始未乘以10年系数)---------------------------

% trend of 12 months
mon_trend_dR_q = zeros(nlon, nlat, 12, 2); 
mon_trend_dR_alb = mon_trend_dR_q;
mon_trend_dR_ts = mon_trend_dR_q;
mon_trend_dR_t = mon_trend_dR_q;
mon_trend_dR_cloud = mon_trend_dR_q;
% f test
p_dR_q = zeros(nlon, nlat, 12); p_dR_alb = p_dR_q; p_dR_ts = p_dR_q; 
p_dR_t = p_dR_q; p_dR_cloud = p_dR_q; 
% mean
cons_dR_q = zeros(nlon, nlat, 12); cons_dR_alb = cons_dR_q; cons_dR_ts = cons_dR_q; 
cons_dR_t = cons_dR_q; cons_dR_cloud = cons_dR_q; 

% trend of season(im==1 is MAM, im==2 is JJA, im==3 is SON, im==4 is DJF)
sea_trend_dR_q = zeros(nlon, nlat, 4); 
sea_trend_dR_alb = sea_trend_dR_q;
sea_trend_dR_ts = sea_trend_dR_q;
sea_trend_dR_t = sea_trend_dR_q;
sea_trend_dR_cloud = sea_trend_dR_q;

% trend of year
yr_trend_dR_q = zeros(nlon, nlat); 
yr_trend_dR_alb = yr_trend_dR_q;
yr_trend_dR_ts = yr_trend_dR_q;
yr_trend_dR_t = yr_trend_dR_q;
yr_trend_dR_cloud = yr_trend_dR_q;

% sfc
% % months mean trend(unite:per day)
% for i = 1:nlon
%     for j = 1:nlat
%         [~, mon_trend_dR_q(i, j, :, :), cons_dR_q(i, j, :), p_dR_q(i, j, :)] = detrend_yan(dR_husSfcAll(i, j, :), time);
%         [~, mon_trend_dR_alb(i, j, :, :), cons_dR_alb(i, j, :), p_dR_alb(i, j, :)] = detrend_yan(dR_albSfcAll(i, j, :), time);
%         [~, mon_trend_dR_ts(i, j, :, :), cons_dR_ts(i, j, :), p_dR_ts(i, j, :)] = detrend_yan(dR_tsSfcAll(i, j, :), time);
%         [~, mon_trend_dR_t(i, j, :, :), cons_dR_t(i, j, :), p_dR_t(i, j, :)] = detrend_yan(dR_taSfcAll(i, j, :), time);
%         [~, mon_trend_dR_cloud(i, j, :, :), cons_dR_cloud(i, j, :), p_dR_cloud(i, j, :)] = detrend_yan(dR_cloud_sfc(i, j, :), time);
%     end
% end

% toa
% months mean trend(unite:per day)
for i = 1:nlon
    for j = 1:nlat
        [~, mon_trend_dR_q(i, j, :, :), cons_dR_q(i, j, :), p_dR_q(i, j, :)] = detrend_yan(dR_qToaAll(i, j, :), time);
        [~, mon_trend_dR_alb(i, j, :, :), cons_dR_alb(i, j, :), p_dR_alb(i, j, :)] = detrend_yan(dR_albToaAll(i, j, :), time);
        [~, mon_trend_dR_ts(i, j, :, :), cons_dR_ts(i, j, :), p_dR_ts(i, j, :)] = detrend_yan(dR_tsToaAll(i, j, :), time);
        [~, mon_trend_dR_t(i, j, :, :), cons_dR_t(i, j, :), p_dR_t(i, j, :)] = detrend_yan(dR_taToaAll(i, j, :), time);
        [~, mon_trend_dR_cloud(i, j, :, :), cons_dR_cloud(i, j, :), p_dR_cloud(i, j, :)] = detrend_yan(dR_cloud_toa(i, j, :), time);
    end
end

% seasons mean trend (unite:per day)
for ii = 1:4
    if ii ~= 4
     sea_trend_dR_q(:,:,ii) = squeeze(mean(mon_trend_dR_q(:,:,ii*3:ii*3+2,1),3));
     sea_trend_dR_alb(:,:,ii) = squeeze(mean(mon_trend_dR_alb(:,:,ii*3:ii*3+2,1),3));
     sea_trend_dR_ts(:,:,ii) = squeeze(mean(mon_trend_dR_ts(:,:,ii*3:ii*3+2,1),3));
     sea_trend_dR_t(:,:,ii) = squeeze(mean(mon_trend_dR_t(:,:,ii*3:ii*3+2,1),3));
     sea_trend_dR_cloud(:,:,ii) = squeeze(mean(mon_trend_dR_cloud(:,:,ii*3:ii*3+2,1),3));
    elseif ii == 4
       ii1=[1,2,12];
       sea_trend_dR_q(:,:,ii) = squeeze(mean(mon_trend_dR_q(:,:,ii1,1),3));
       sea_trend_dR_alb(:,:,ii) = squeeze(mean(mon_trend_dR_alb(:,:,ii1,1),3));
       sea_trend_dR_ts(:,:,ii) = squeeze(mean(mon_trend_dR_ts(:,:,ii1,1),3));
       sea_trend_dR_t(:,:,ii) = squeeze(mean(mon_trend_dR_t(:,:,ii1,1),3));
       sea_trend_dR_cloud(:,:,ii) = squeeze(mean(mon_trend_dR_cloud(:,:,ii1,1),3));
    end
end

% years mean trend (unite:per day)
yr_trend_dR_q = squeeze(mean(mon_trend_dR_q(:,:,:,1),3));
yr_trend_dR_alb = squeeze(mean(mon_trend_dR_alb(:,:,:,1),3));
yr_trend_dR_ts = squeeze(mean(mon_trend_dR_ts(:,:,:,1),3));
yr_trend_dR_t = squeeze(mean(mon_trend_dR_t(:,:,:,1),3));
yr_trend_dR_cloud = squeeze(mean(mon_trend_dR_cloud(:,:,:,1),3));



%----------------------save varity----------------------------------
% outfilename = 'E:\Repeat_the_experiment\testdata\radcomponent_data\trendradset_erai_allskysfc';
outfilename = 'E:\Repeat_the_experiment\testdata\radcomponent_data\trendradset_erai_allskytoa';
save(outfilename, 'lone', 'late', 'time', ...
'mon_trend_dR_q','mon_trend_dR_alb','mon_trend_dR_ts','mon_trend_dR_t','mon_trend_dR_cloud',...
'sea_trend_dR_q','sea_trend_dR_alb','sea_trend_dR_ts','sea_trend_dR_t','sea_trend_dR_cloud',...
'yr_trend_dR_q','yr_trend_dR_alb','yr_trend_dR_ts','yr_trend_dR_t','yr_trend_dR_cloud')

t = toc; disp(t)