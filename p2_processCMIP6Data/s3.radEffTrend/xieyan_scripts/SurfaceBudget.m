clc
clear
% This program is to check the surface energy budget balance situation
filepath = 'C:\Yourfilepath\';
filename1 = strcat(filepath,'Rad_Chinamean.nc');
filename2 = strcat(filepath,'Radsupp_Chinamean.nc');
filename3 = strcat(filepath,'ERAiHF_Chinamean.nc');
diag0 = ncread(filename2,'dR1mean'); diag0t = ncread(filename2,'dR2mean');
diag = squeeze(diag0(2,1,1,2,:,:)); diagt = squeeze(diag0t(2,1,1,2,:,:));
res_kernel= squeeze(diag0(3,1,1,2,:,:)); res_kernelt = squeeze(diag0t(3,1,1,2,:,:));
flux0 = ncread(filename3,'dFlux_mean');
flux0t = ncread(filename3,'dFluxmean_detrend');
lahf = squeeze(flux0(1,:,:)); lahft = squeeze(flux0t(1,:,:)); 
sehf = squeeze(flux0(2,:,:)); sehft = squeeze(flux0t(2,:,:));
residual0 = -(diag + lahf + sehf); residual0t = -(diagt+lahft+sehft);
res_budget = residual0 - res_kernel; res_budgett = residual0t - res_kernelt;
time = ncread(filename1,'time'); ntime = length(time);  
%% Radiation Anomalies
dRx1 = ncread(filename1,'dR1mean');                                      % Radiative feedback anomaly time series
dRx1t = ncread(filename1,'dR2mean');
dRx2 = ncread(filename2,'dR1mean');                                      % x,sky,spectrum,budget,region,time 
dRx2t = ncread(filename2,'dR2mean');
region = 3;                                                              % China, QTP, China excluding QTP
%% ALL Dimension:11
dR_sfc0 = zeros(12,ntime);                                                % x
dR_sfc0(1,:) = -squeeze(dRx2(5,1,1,2,region,:));                          % ts  
dR_sfc0(2,:) = squeeze(dRx2(4,1,1,2,region,:));                           % lw ta
dR_sfc0(3,:) = squeeze(dRx1(3,1,2,2,region,:));                           % lw wv
dR_sfc0(4,:) = squeeze(dRx1(3,1,3,2,region,:));                           % sw wv
dR_sfc0(5,:) = squeeze(dRx1(4,1,3,2,region,:));                           % sw alb
dR_sfc0(6,:) = squeeze(dRx1(5,1,2,2,region,:));                           % lw cloud
dR_sfc0(7,:) = squeeze(dRx1(5,1,3,2,region,:));                           % sw cloud
dR_sfc0(8,:) = squeeze(dRx2(3,1,2,2,region,:));                           % lw kernel residual
dR_sfc0(9,:) = squeeze(dRx2(3,1,3,2,region,:));                           % sw kernel residual
dR_sfc0(10,:) = squeeze(lahf(region,:));                                  % latent heat
dR_sfc0(11,:) = squeeze(sehf(region,:));                                  % sensible heat
dR_sfc0(12,:) = squeeze(res_budget(region,:));                            % energy budget residual

dR_sfc0t = zeros(12,ntime);                                                % x
dR_sfc0t(1,:) = -squeeze(dRx2t(5,1,1,2,region,:));                          % ts  
dR_sfc0t(2,:) = squeeze(dRx2t(4,1,1,2,region,:));                           % lw ta
dR_sfc0t(3,:) = squeeze(dRx1t(3,1,2,2,region,:));                           % lw wv
dR_sfc0t(4,:) = squeeze(dRx1t(3,1,3,2,region,:));                           % sw wv
dR_sfc0t(5,:) = squeeze(dRx1t(4,1,3,2,region,:));                           % sw alb
dR_sfc0t(6,:) = squeeze(dRx1t(5,1,2,2,region,:));                           % lw cloud
dR_sfc0t(7,:) = squeeze(dRx1t(5,1,3,2,region,:));                           % sw cloud
dR_sfc0t(8,:) = squeeze(dRx2t(3,1,2,2,region,:));                           % lw kernel residual
dR_sfc0t(9,:) = squeeze(dRx2t(3,1,3,2,region,:));                           % sw kernel residual
dR_sfc0t(10,:) = squeeze(lahft(region,:));                                  % latent heat
dR_sfc0t(11,:) = squeeze(sehft(region,:));                                  % sensible heat
dR_sfc0t(12,:) = squeeze(res_budgett(region,:));                            % energy budget residual
trend_sfc0 = dR_sfc0 - dR_sfc0t;
%% All-sky / NET / Surface
dRsfc1 = squeeze(dR_sfc0(2:12,:)); dRsfc1 = dRsfc1';
cov_sfc1= cov(dRsfc1);
contri_sfc1 = cov_sfc1/var(squeeze(dR_sfc0(1,:)));
contri_sfc10 = zeros(size(contri_sfc1));
for jj = 1:11
   for ii = jj:11
       if ii ==  jj
          contri_sfc10(ii,jj) = contri_sfc1(ii,jj);
       else
          contri_sfc10(ii,jj) = 2*contri_sfc1(ii,jj);
       end
   end
end
contri_sfc20 = flipud(contri_sfc10); 

% Trend Part
dTsfc1 = squeeze(var(trend_sfc0(2:12,:),0,2)); 
contri_sfc40 = zeros(1,11);
for ii = 1:11
    contri_sfc40(1,ii) = dTsfc1(ii,1)./sum(dTsfc1(1:11,1)); 
end