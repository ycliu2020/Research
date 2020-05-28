% calculation including Ts, Rheating, cloud, Ta, wv, albedo radiation(ERAi date)
% Ts, Ta, wv, albedo can be caled by kernels
% cloud must be caled by the method in xieyan's paper
% Rheating can be caled by raw date about down lw plus net sw (trend already cal)




clc;clear;tic;
add_all;
input_path='E:\Repeat_the_experiment\testdata\ERAi\rad_effect\';%用ERA辐射计算云辐射效应
input_path_CERES='E:\Repeat_the_experiment\testdata\CERES_EBAF\';% 用CERES辐射计算云辐射效应

kind_eff = {[input_path,'eraiRadAnom_2018SFCallsky_r.nc'],[input_path,'eraiRadAnom_2018TOAallsky_r.nc'],...
[input_path,'eraiRadAnom_2018SFCclrsky_r.nc'],[input_path,'eraiRadAnom_2018TOAclrsky_r.nc']}; % radiative effect

kind_rrad_sfc ={[input_path_CERES,'ceresAnomclrSFC.nc'],[input_path_CERES,'ceresAnomSFC.nc']} ;% ceres sfc real radiation
kind_rrad_toa={[input_path_CERES,'ceresAnomclrTOA.nc'],[input_path_CERES,'ceresAnomTOA.nc']} ;% ceres toa real radiation
% kind_rrad ={[input_path,'ERAi_radAnom_clr.nc'],[input_path,'ERAi_radAnom_all.nc']} ;% erai real radiation
f_numb = 1; % 1.clr_sfc, 2.clr_toa, 3.all_sfc,4.all_toa,

% --------------------------read effects-------------------------
lone = ncread(kind_eff{f_numb},'longitude'); % 144x72
late = ncread(kind_eff{f_numb},'latitude');
nlon = length(lone);nlat = length(late);
time = ncread(kind_eff{f_numb},'time');ntime = length(time);
% radiative effect(144,72,216,2)sfc1,toa2
s_albeff_all=zeros(144,72,216,2);l_teff_all=s_albeff_all; l_tseff_all=s_albeff_all;l_wveff_all=s_albeff_all;s_wveff_all=s_albeff_all;totaleff_all=s_albeff_all;
s_albeff_clr=zeros(144,72,216,2);l_teff_clr=s_albeff_clr; l_tseff_clr=s_albeff_clr;l_wveff_clr=s_albeff_clr;s_wveff_clr=s_albeff_clr;totaleff_clr=s_albeff_clr;
for i=1:4
    if i<=2
        s_albeff_all(:,:,:,i) = ncread(kind_eff{i},'albRadEff');
        l_teff_all(:,:,:,i) = ncread(kind_eff{i},'tRadEff');
        l_tseff_all(:,:,:,i) = ncread(kind_eff{i},'tsRadEff');
        l_wveff_all(:,:,:,i) = ncread(kind_eff{i},'wvlwRadEff');
        s_wveff_all(:,:,:,i) = ncread(kind_eff{i},'wvswRadEff');
        totaleff_all(:,:,:,i) = ncread(kind_eff{i},'totalRadEff');
    else
        s_albeff_clr(:,:,:,i-2) = ncread(kind_eff{i},'albRadEff');
        l_teff_clr(:,:,:,i-2) = ncread(kind_eff{i},'tRadEff');
        l_tseff_clr(:,:,:,i-2) = ncread(kind_eff{i},'tsRadEff');
        l_wveff_clr(:,:,:,i-2) = ncread(kind_eff{i},'wvlwRadEff');
        s_wveff_clr(:,:,:,i-2) = ncread(kind_eff{i},'wvswRadEff');
        totaleff_clr(:,:,:,i-2) = ncread(kind_eff{i},'totalRadEff');   
    end
end
% note: only cal the allsky condition, and sfc only used 
% q, alb, ts, t
dR_qeff_all = l_wveff_all+s_wveff_all;
dR_q_sfc_all = dR_qeff_all(:,:,:,1);dR_q_toa_all = dR_qeff_all(:,:,:,2); % q
dR_alb_sfc_all = s_albeff_all(:,:,:,1);dR_alb_toa_all =s_albeff_all(:,:,:,2);%albedo
dR_ts_sfc_all = l_tseff_all(:,:,:,1);dR_ts_toa_all = l_tseff_all(:,:,:,2);% ts
dR_t_sfc_all = l_teff_all(:,:,:,1);dR_t_toa_all = l_teff_all(:,:,:,2);% t 

% real radiation anomaly
l_rad=zeros(144,72,216,4);s_rad=l_rad;

% kind_rrad_sfc ={[input_path_CERES,'ceresAnomclrSFC.nc'],[input_path_CERES,'ceresAnomSFC.nc']} ;% ceres sfc real radiation
% kind_rrad_toa={[input_path_CERES,'ceresAnomclrTOA.nc'],[input_path_CERES,'ceresAnomTOA.nc']} ;% ceres toa real radiation
for i=1:2
    if i==1% clr
        dR_clr(:,:,:,1) = ncread(kind_rrad_sfc{i},'dR');% Surface net  radiation, 
        dR_clr(:,:,:,2) = ncread(kind_rrad_toa{i},'dR');% Top net radiation

    else% allsky
        dR_allsky(:,:,:,1) = ncread(kind_rrad_sfc{i},'dR');
        dR_allsky(:,:,:,2) = ncread(kind_rrad_toa{i},'dR');
    end
end


% ------------------------cal the cloud effect(lw+sw)---------------------------
% dR_cloud(144,72,216,2)(最后一维代表sfc和toa)
% for research only clear and all sky's surface needed (1 and 3 are used)

%% cal the dR_cloud (144,72,216,2)包含地表和大气层顶, 所以最后一维为2, 
dR_res = dR_clr-totaleff_clr;
dR_cloud = dR_allsky - totaleff_all - dR_res;
dR_cloud_sfc = dR_cloud(:,:,:,1);dR_cloud_toa = dR_cloud(:,:,:,2);%云的辐射效应只有有云时存在

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

% %sfc
% % months mean trend(unite:per day)
% for i = 1:nlon
%     for j = 1:nlat
%         [~, mon_trend_dR_q(i, j, :, :), cons_dR_q(i, j, :), p_dR_q(i, j, :)] = detrend_yan(dR_q_sfc_all(i, j, :), time);
%         [~, mon_trend_dR_alb(i, j, :, :), cons_dR_alb(i, j, :), p_dR_alb(i, j, :)] = detrend_yan(dR_alb_sfc_all(i, j, :), time);
%         [~, mon_trend_dR_ts(i, j, :, :), cons_dR_ts(i, j, :), p_dR_ts(i, j, :)] = detrend_yan(dR_ts_sfc_all(i, j, :), time);
%         [~, mon_trend_dR_t(i, j, :, :), cons_dR_t(i, j, :), p_dR_t(i, j, :)] = detrend_yan(dR_t_sfc_all(i, j, :), time);
%         [~, mon_trend_dR_cloud(i, j, :, :), cons_dR_cloud(i, j, :), p_dR_cloud(i, j, :)] = detrend_yan(dR_cloud_sfc(i, j, :), time);
%     end
% end

% toa
% months mean trend(unite:per day)
for i = 1:nlon
    for j = 1:nlat
        [~, mon_trend_dR_q(i, j, :, :), cons_dR_q(i, j, :), p_dR_q(i, j, :)] = detrend_yan(dR_q_toa_all(i, j, :), time);
        [~, mon_trend_dR_alb(i, j, :, :), cons_dR_alb(i, j, :), p_dR_alb(i, j, :)] = detrend_yan(dR_alb_toa_all(i, j, :), time);
        [~, mon_trend_dR_ts(i, j, :, :), cons_dR_ts(i, j, :), p_dR_ts(i, j, :)] = detrend_yan(dR_ts_toa_all(i, j, :), time);
        [~, mon_trend_dR_t(i, j, :, :), cons_dR_t(i, j, :), p_dR_t(i, j, :)] = detrend_yan(dR_t_toa_all(i, j, :), time);
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
% outfilename = 'E:\Repeat_the_experiment\testdata\data_radcomponent\trendradset_cerescld__allskysfc';
outfilename = 'E:\Repeat_the_experiment\testdata\data_radcomponent\trend_radset_cerescld_allskytoa';
save(outfilename, 'lone', 'late', 'time', ...
'mon_trend_dR_q','mon_trend_dR_alb','mon_trend_dR_ts','mon_trend_dR_t','mon_trend_dR_cloud',...
'sea_trend_dR_q','sea_trend_dR_alb','sea_trend_dR_ts','sea_trend_dR_t','sea_trend_dR_cloud',...
'yr_trend_dR_q','yr_trend_dR_alb','yr_trend_dR_ts','yr_trend_dR_t','yr_trend_dR_cloud')

t = toc; disp(t)