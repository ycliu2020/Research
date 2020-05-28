clc
clear
% This program is to calculate linear trend of time series
filepath = 'C:\Yourfilepath\';
% filename = strcat(filepath,'ERAi_Tsmean.nc');
% filename = strcat(filepath,'HadCRUT4_Tsmean.nc');
filename = strcat(filepath,'MERRA2_Tsmean.nc');
tschm = ncread(filename,'tschm');
tslandm = ncread(filename,'tslandm');
tsmidm = ncread(filename,'tsmidm');
tsnhm = ncread(filename,'tsnhm');
time = ncread(filename,'time');
% tschm = tschm(1:468,1);
% tslandm = tslandm(1:468,1);
% tsmidm = tsmidm(1:468,1);
% tsnhm = tsnhm(1:468,1);
% time = time(1:468,1);
ntime = length(time);

pch = polyfit(time,tschm,1);
plan = polyfit(time,tslandm,1);
pmid = polyfit(time,tsmidm,1);
pnh = polyfit(time,tsnhm,1);
x = time;
ych = pch(1,1)*x+pch(1,2);
ylan = plan(1,1)*x+plan(1,2);
ymid = pmid(1,1)*x+pmid(1,2);
ynh = pnh(1,1)*x+pnh(1,2);
% trend per annual
tch = pch(1,1)*365*10; tlan = plan(1,1)*365*10; 
tmid = pmid(1,1)*365*10; tnh = pnh(1,1)*365*10;
tch = roundn(tch,-3); tlan = roundn(tlan,-3);
tmid = roundn(tmid,-3); tnh = roundn(tnh,-3);
%% Plot
figure(1)
% China area-mean
subplot_yih(4,1,1,1)
hold all
plot(time,tschm,'r')
plot(time,ych,'k')
datetick('x','yyyy','keepticks');
xlim([time(1) time(ntime)])
ylim([-4 4])
ylabel('Ts Anomaly (K)')
legend({['China trend = ', num2str(tch),'°Ê/10 year']},'Fontsize',8,'Location', 'northwest')

% 20 - 50N area mean
subplot_yih(4,1,2,1)
hold all
plot(time,tsmidm,'r')
plot(time,ymid,'k')
datetick('x','yyyy','keepticks');
xlim([time(1) time(ntime)])
ylim([-4 4])
ylabel('Ts Anomaly (K)')
legend({['20°„N-50°„N trend = ', num2str(tmid),'°Ê/10 year']},'Fontsize',8,'Location', 'northwest')

% Northern Hemisphere area mean
subplot_yih(4,1,3,1)
hold all
plot(time,tsnhm,'r')
plot(time,ynh,'k')
datetick('x','yyyy','keepticks');
xlim([time(1) time(ntime)])
ylim([-4 4])
ylabel('Ts Anomaly (K)')
legend({['Northern Hemisphere trend = ', num2str(tnh),'°Ê/10 year']},'Fontsize',8,'Location', 'northwest')

% Global Land area mean
subplot_yih(4,1,4,1)
hold all
plot(time,tslandm,'r')
plot(time,ylan,'k')
datetick('x','yyyy','keepticks');
xlim([time(1) time(ntime)])
ylim([-4 4])
ylabel('Ts Anomaly (K)')
legend({['Global land trend = ', num2str(tlan),'°Ê/10 year']},'Fontsize',8,'Location', 'northwest')

