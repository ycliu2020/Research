clc
% This program is to plot radiation anomalies time series 
% From 2000/03 - 2018/02
clear
vector_var_str{1} = 'net';
vector_var_str{2} = 'lw';
vector_var_str{3} = 'sw';

filepath = 'C:\Yourfilepath\';
filename1 = strcat(filepath,'Rad_Chinamean.nc');
filename2 = strcat(filepath,'Radsupp_Chinamean.nc');
time = ncread(filename1,'time'); ntime = length(time);  

%% Radiation Anomalies
dRx1 = ncread(filename1,'dR1mean');                                      % Radiative feedback anomaly time series
dRx1t = ncread(filename1,'dR2mean');
dRx2 = ncread(filename2,'dR1mean');                                      % x,sky,spectrum,budget,region,time 
dRx2t = ncread(filename2,'dR2mean');
region = 1;%2/3
%% ALL Dimension:6
% dR_sfc0 = zeros(7,ntime);                                                 % x
% dR_sfc0(2,:) = squeeze(dRx2(4,1,2,2,region,:));                           % lw ta
% dR_sfc0(3,:) = squeeze(dRx1(3,1,2,2,region,:));                           % lw wv
% dR_sfc0(4,:) = squeeze(dRx1(3,1,3,2,region,:));                           % sw wv
% dR_sfc0(5,:) = squeeze(dRx1(4,1,3,2,region,:));                           % sw alb
% dR_sfc0(6,:) = squeeze(dRx1(5,1,2,2,region,:));                           % lw cloud
% dR_sfc0(7,:) = squeeze(dRx1(5,1,3,2,region,:));                           % sw cloud
% dR_sfc0(1,:) = sum(dR_sfc0(2:7,:),1);                                     % Surface Heat
% dR_sfct0 = zeros(7,ntime);                                                 % x
% dR_sfct0(2,:) = squeeze(dRx2t(4,1,2,2,region,:));                           % lw ta
% dR_sfct0(3,:) = squeeze(dRx1t(3,1,2,2,region,:));                           % lw wv
% dR_sfct0(4,:) = squeeze(dRx1t(3,1,3,2,region,:));                           % sw wv
% dR_sfct0(5,:) = squeeze(dRx1t(4,1,3,2,region,:));                           % sw alb
% dR_sfct0(6,:) = squeeze(dRx1t(5,1,2,2,region,:));                           % lw cloud
% dR_sfct0(7,:) = squeeze(dRx1t(5,1,3,2,region,:));                           % sw cloud
% dR_sfct0(1,:) = sum(dR_sfct0(2:7,:),1);                                     % Surface Heat
% trend_sfc0 = dR_sfc0 - dR_sfct0;
% 
% dR_toa0 = zeros(7,ntime);                                                 % x
% dR_toa0(2,:) = squeeze(dRx2(4,1,2,1,region,:));                           % lw ta
% dR_toa0(3,:) = squeeze(dRx1(3,1,2,1,region,:));                           % lw wv
% dR_toa0(4,:) = squeeze(dRx1(3,1,3,1,region,:));                           % sw wv
% dR_toa0(5,:) = squeeze(dRx1(4,1,3,1,region,:));                           % sw alb
% dR_toa0(6,:) = squeeze(dRx1(5,1,2,1,region,:));                           % lw cloud
% dR_toa0(7,:) = squeeze(dRx1(5,1,3,1,region,:));                           % sw cloud
% dR_toa0(1,:) = sum(dR_toa0(2:7,:),1);                                     % Surface Heat
% dR_toat0 = zeros(7,ntime);                                                 % x
% dR_toat0(2,:) = squeeze(dRx2t(4,1,2,1,region,:));                           % lw ta
% dR_toat0(3,:) = squeeze(dRx1t(3,1,2,1,region,:));                           % lw wv
% dR_toat0(4,:) = squeeze(dRx1t(3,1,3,1,region,:));                           % sw wv
% dR_toat0(5,:) = squeeze(dRx1t(4,1,3,1,region,:));                           % sw alb
% dR_toat0(6,:) = squeeze(dRx1t(5,1,2,1,region,:));                           % lw cloud
% dR_toat0(7,:) = squeeze(dRx1t(5,1,3,1,region,:));                           % sw cloud
% dR_toat0(1,:) = sum(dR_toat0(2:7,:),1);                                     % Surface Heat
% trend_toa0 = dR_toa0 - dR_toat0;

dR_sfc0 = zeros(7,ntime);                                                 % x
dR_sfc0(2,:) = squeeze(dRx2(4,1,2,2,region,:));                           % lw ta
dR_sfc0(3,:) = squeeze(dRx1(3,1,2,2,region,:));                           % lw wv
dR_sfc0(4,:) = squeeze(dRx1(3,1,3,2,region,:));                           % sw wv
dR_sfc0(5,:) = squeeze(dRx1(4,1,3,2,region,:));                           % sw alb
dR_sfc0(6,:) = squeeze(dRx1(5,1,2,2,region,:));                           % lw cloud
dR_sfc0(7,:) = squeeze(dRx1(5,1,3,2,region,:));                           % sw cloud
dR_sfc0(1,:) = sum(dR_sfc0(2:7,:),1);                                     % Surface Heat
dR_sfct0 = zeros(7,ntime);                                                 % x
dR_sfct0(2,:) = squeeze(dRx2t(4,1,2,2,region,:));                           % lw ta
dR_sfct0(3,:) = squeeze(dRx1t(3,1,2,2,region,:));                           % lw wv
dR_sfct0(4,:) = squeeze(dRx1t(3,1,3,2,region,:));                           % sw wv
dR_sfct0(5,:) = squeeze(dRx1t(4,1,3,2,region,:));                           % sw alb
dR_sfct0(6,:) = squeeze(dRx1t(5,1,2,2,region,:));                           % lw cloud
dR_sfct0(7,:) = squeeze(dRx1t(5,1,3,2,region,:));                           % sw cloud
dR_sfct0(1,:) = sum(dR_sfct0(2:7,:),1);                                     % Surface Heat
trend_sfc0 = dR_sfc0 - dR_sfct0;

dR_toa0 = zeros(7,ntime);                                                 % x
dR_toa0(2,:) = squeeze(dRx2(4,1,2,1,region,:));                           % lw ta
dR_toa0(3,:) = squeeze(dRx1(3,1,2,1,region,:));                           % lw wv
dR_toa0(4,:) = squeeze(dRx1(3,1,3,1,region,:));                           % sw wv
dR_toa0(5,:) = squeeze(dRx1(4,1,3,1,region,:));                           % sw alb
dR_toa0(6,:) = squeeze(dRx1(5,1,2,1,region,:));                           % lw cloud
dR_toa0(7,:) = squeeze(dRx1(5,1,3,1,region,:));                           % sw cloud
dR_toa0(1,:) = sum(dR_toa0(2:7,:),1);                                     % Surface Heat
dR_toat0 = zeros(7,ntime);                                                 % x
dR_toat0(2,:) = squeeze(dRx2t(4,1,2,1,region,:));                           % lw ta
dR_toat0(3,:) = squeeze(dRx1t(3,1,2,1,region,:));                           % lw wv
dR_toat0(4,:) = squeeze(dRx1t(3,1,3,1,region,:));                           % sw wv
dR_toat0(5,:) = squeeze(dRx1t(4,1,3,1,region,:));                           % sw alb
dR_toat0(6,:) = squeeze(dRx1t(5,1,2,1,region,:));                           % lw cloud
dR_toat0(7,:) = squeeze(dRx1t(5,1,3,1,region,:));                           % sw cloud
dR_toat0(1,:) = sum(dR_toat0(2:7,:),1);                                     % Surface Heat
trend_toa0 = dR_toa0 - dR_toat0;
%% NET
% Radiation Anomalies
dRsfc1 = squeeze(dR_sfc0(2:7,:)); dRsfc1 = dRsfc1';
dRtoa1 = squeeze(dR_toa0(2:7,:)); dRtoa1 = dRtoa1';
cov_sfc1= cov(dRsfc1);
cov_toa1= cov(dRtoa1);
contri_sfc1 = cov_sfc1/var(squeeze(dR_sfc0(1,:)));
contri_toa1 = cov_toa1/var(squeeze(dR_toa0(1,:)));
contri_sfc10 = zeros(size(contri_sfc1));
contri_toa10 = zeros(size(contri_toa1));
for jj = 1:6
   for ii = jj:6
       if ii ==  jj
          contri_sfc10(ii,jj) = contri_sfc1(ii,jj);
          contri_toa10(ii,jj) = contri_toa1(ii,jj);
       else
          contri_sfc10(ii,jj) = 2*contri_sfc1(ii,jj);
          contri_toa10(ii,jj) = 2*contri_toa1(ii,jj);
       end
   end
end
contri_sfc20 = flipud(contri_sfc10); contri_toa20 = flipud(contri_toa10);

% Trend Part
dTsfc1 = squeeze(var(trend_sfc0(2:7,:),0,2)); 
dTtoa1 = squeeze(var(trend_toa0(2:7,:),0,2)); 
contri_sfc40 = zeros(1,6);
contri_toa40 = zeros(1,6);
for ii = 1:6
    contri_sfc40(1,ii) = dTsfc1(ii,1)./sum(dTsfc1(1:6,1)); 
    contri_toa40(1,ii) = dTtoa1(ii,1)./sum(dTtoa1(1:6,1));
end
