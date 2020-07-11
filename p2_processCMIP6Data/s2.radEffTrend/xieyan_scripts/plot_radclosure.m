clc
clear
% This program is to test radiation closure of kernel method
% For four regions: China/ 20-50N/ Northern Hemisphere land/ Global land

filepath = 'C:\Yourfilepath\';
filename1 = strcat(filepath,'Rad_Chinamean.nc');
filename2 = strcat(filepath,'Radsupp_Chinamean.nc');
time = ncread(filename1,'time'); ntime = length(time);
dR_all = ncread(filename1,'dR1mean'); dR_all = squeeze(dR_all(1,2,:,:,:,:)); % total, clear-sky, spectrum, budget, region,time
dR_noncloud = ncread(filename2,'dR1mean'); dR_noncloud = squeeze(dR_noncloud(1,2,:,:,:,:));
dR_res = ncread(filename2,'dR1mean'); dR_res = squeeze(dR_res(3,2,:,:,:,:));
dR0 = NaN(3,3,2,3,ntime);                                                 % x,spectrum,budget,region,time
dR0(1,:,:,:,:) = dR_all; dR0(2,:,:,:,:) = dR_noncloud; dR0(3,:,:,:,:) = dR_res;
region = {'China','Qinghai-Tibet Plateau',' China exclude Qinghai-Tibet Plateau'};
for ii = 1:3%3 % region
    reg = region(ii);
    dR_sfc= squeeze(dR0(:,:,2,ii,:)); 
    dR_toa = squeeze(dR0(:,:,1,ii,:));
    cc_sfc = zeros(3,1);    cc_toa = zeros(3,1);                          % LW/ SW/ NET
    for jj = 1:3
       temp1 = squeeze(dR_sfc(1,jj,:)); 
       temp2 = squeeze(dR_sfc(2,jj,:));
       cc0 = corrcoef(temp1,temp2,'Rows','complete');
       cc_sfc(jj,1) = cc0(1,2);
       cc_sfc(jj,1) = roundn(cc_sfc(jj,1),-3);
       
       temp1 = squeeze(dR_toa(1,jj,:)); 
       temp2 = squeeze(dR_toa(2,jj,:));
       cc0 = corrcoef(temp1,temp2,'Rows','complete');
       cc_toa(jj,1) = cc0(1,2);
       cc_toa(jj,1) = roundn(cc_toa(jj,1),-3);
    end
    
    figure(ii)
    %% SFC
    % NET
    subplot_yih(3,2,1,1)
    hold all
    plot(time,squeeze(dR_sfc(1,1,:)),'b')
    plot(time,squeeze(dR_sfc(2,1,:)),'r')
    plot(time,squeeze(dR_sfc(3,1,:)),'k')
    datetick('x','yyyy','keepticks');
    xlim([time(1) time(ntime)])
    xticks([732358 734153 735979])
    xticklabels({'2005','2010','2015'})
    set(gca,'XMinorTick','on'); 
    set(gca,'ticklength',[0.02 0.01]);
    ylim([-15 15])
    set(gca,'YMinorTick','on'); 
    set(gca,'ticklength',[0.02 0.01]);
    set(gca,'ycolor','k')
    ylabel('SFC Rad Anomaly (Wm^{-2})','Fontsize',10)
    legend({'dR_{NET} ', ['Diagnosed cc =',num2str(cc_sfc(1,1))], 'Residual'},'Fontsize',8,'Location', 'northwest')
    legend boxoff
    % LW
    subplot_yih(3,2,2,1)
    hold all
    plot(time,squeeze(dR_sfc(1,2,:)),'b')
    plot(time,squeeze(dR_sfc(2,2,:)),'r')
    plot(time,squeeze(dR_sfc(3,2,:)),'k')
    datetick('x','yyyy','keepticks');
    xlim([time(1) time(ntime)])
    xticks([732358 734153 735979])
    xticklabels({'2005','2010','2015'})
    set(gca,'XMinorTick','on'); 
    set(gca,'ticklength',[0.02 0.01]);
    ylim([-10 10])
    set(gca,'YMinorTick','on'); 
    set(gca,'ticklength',[0.02 0.01]);
    set(gca,'ycolor','k')
    ylabel('SFC Rad Anomaly (Wm^{-2})','Fontsize',10)
    legend({'dR_{LW} ', ['Diagnosed cc =',num2str(cc_sfc(2,1))], 'Residual'},'Fontsize',8,'Location', 'northwest')
    legend boxoff
    % SW
    subplot_yih(3,2,3,1)
    hold all
    plot(time,squeeze(dR_sfc(1,3,:)),'b')
    plot(time,squeeze(dR_sfc(2,3,:)),'r')
    plot(time,squeeze(dR_sfc(3,3,:)),'k')
    datetick('x','yyyy','keepticks');
    xlim([time(1) time(ntime)])
   xticks([732358 734153 735979])
    xticklabels({'2005','2010','2015'})
    set(gca,'XMinorTick','on'); 
    set(gca,'ticklength',[0.02 0.01]);
    ylim([-15 15])
    set(gca,'YMinorTick','on'); 
    set(gca,'ticklength',[0.02 0.01]);
    set(gca,'ycolor','k')
    ylabel('SFC Rad Anomaly (Wm^{-2})','Fontsize',10)
    legend({'dR_{SW} ', ['Diagnosed cc =',num2str(cc_sfc(3,1))], 'Residual'},'Fontsize',8,'Location', 'northwest')
    legend boxoff
    
    %% TOA
    subplot_yih(3,2,1,2)
    hold all
    plot(time,squeeze(dR_toa(1,1,:)),'b')
    plot(time,squeeze(dR_toa(2,1,:)),'r')
    plot(time,squeeze(dR_toa(3,1,:)),'k')
    datetick('x','yyyy','keepticks');
    xlim([time(1) time(ntime)])
    xticks([732358 734153 735979])
    xticklabels({'2005','2010','2015'})
    set(gca,'XMinorTick','on'); 
    set(gca,'ticklength',[0.02 0.01]);
    ylim([-15 15])
    set(gca,'YMinorTick','on'); 
    set(gca,'ticklength',[0.02 0.01]);
    set(gca,'ycolor','k')
    ylabel('TOA Rad Anomaly (Wm^{-2})','Fontsize',10)
    legend({'dR_{NET} ', ['Diagnosed cc =',num2str(cc_toa(1,1))], 'Residual'},'Fontsize',8,'Location', 'northwest')
    legend boxoff
    % LW
    subplot_yih(3,2,2,2)
    hold all
    plot(time,squeeze(dR_toa(1,2,:)),'b')
    plot(time,squeeze(dR_toa(2,2,:)),'r')
    plot(time,squeeze(dR_toa(3,2,:)),'k')
    datetick('x','yyyy','keepticks');
    xlim([time(1) time(ntime)])
    xticks([732358 734153 735979])
    xticklabels({'2005','2010','2015'})
    set(gca,'XMinorTick','on'); 
    set(gca,'ticklength',[0.02 0.01]);
    ylim([-10 10])
    set(gca,'YMinorTick','on'); 
    set(gca,'ticklength',[0.02 0.01]);
    set(gca,'ycolor','k')
    ylabel('TOA Rad Anomaly (Wm^{-2})','Fontsize',10)
    legend({'dR_{LW} ', ['Diagnosed cc =',num2str(cc_toa(2,1))], 'Residual'},'Fontsize',8,'Location', 'northwest')
    legend boxoff
    % SW
    subplot_yih(3,2,3,2)
    hold all
    plot(time,squeeze(dR_toa(1,3,:)),'b')
    plot(time,squeeze(dR_toa(2,3,:)),'r')
    plot(time,squeeze(dR_toa(3,3,:)),'k')
    datetick('x','yyyy','keepticks');
    xlim([time(1) time(ntime)])
    xticks([732358 734153 735979])
    xticklabels({'2005','2010','2015'})
    set(gca,'XMinorTick','on'); 
    set(gca,'ticklength',[0.02 0.01]);
    ylim([-15 15])
    set(gca,'YMinorTick','on'); 
    set(gca,'ticklength',[0.02 0.01]);
    set(gca,'ycolor','k')
    ylabel('TOA Rad Anomaly (Wm^{-2})','Fontsize',10)
    legend({'dR_{SW} ', ['Diagnosed cc =',num2str(cc_toa(3,1))], 'Residual'},'Fontsize',8,'Location', 'northwest')
    legend boxoff
   %  Title
%     set_font_size(14)
%     suptitle(reg)
end 
