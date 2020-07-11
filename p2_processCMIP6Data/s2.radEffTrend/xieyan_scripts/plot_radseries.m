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
    
dRx1 = ncread(filename1,'dR1mean');                                      % Radiative feedback anomaly time series
dRx2 = ncread(filename2,'dR1mean');                                      % x,sky,spectrum,budget,region,time 
region = 3;%2/3
dR_sfc = zeros(7,3,ntime);                                              % x: total,ta,ts,wv,albedo,cloud,residual
dR_sfc(1,:,:) = squeeze(dRx1(1,1,:,2,region,:));
dR_sfc(2:3,:,:) = squeeze(dRx2(4:5,1,:,2,region,:));
dR_sfc(4:6,:,:) = squeeze(dRx1(3:5,1,:,2,region,:));
dR_sfc(7,:,:) = squeeze(dRx2(3,1,:,2,region,:));
dR_toa = zeros(7,3,ntime);                                              % x: total,ta,ts,wv,albedo,cloud,residual
dR_toa(1,:,:) = squeeze(dRx1(1,1,:,1,region,:));
dR_toa(2:3,:,:) = squeeze(dRx2(4:5,1,:,1,region,:));
dR_toa(4:6,:,:) = squeeze(dRx1(3:5,1,:,1,region,:));
dR_toa(7,:,:) = squeeze(dRx2(3,1,:,1,region,:));

ccsfc = zeros(3,1);
temp1 = corrcoef(squeeze(dR_sfc(1,1,:)),squeeze(dR_sfc(6,1,:)),'Rows','complete');
ccsfc(1,1) = temp1(1,2); ccsfc(1,1) = roundn(ccsfc(1,1),-3);
temp1 = corrcoef(squeeze(dR_sfc(1,2,:)),squeeze(dR_sfc(6,2,:)),'Rows','complete');
ccsfc(2,1) = temp1(1,2); ccsfc(2,1) = roundn(ccsfc(2,1),-3);
temp1 = corrcoef(squeeze(dR_sfc(1,3,:)),squeeze(dR_sfc(6,3,:)),'Rows','complete');
ccsfc(3,1) = temp1(1,2); ccsfc(3,1) = roundn(ccsfc(3,1),-3);

cctoa = zeros(3,1);
temp1 = corrcoef(squeeze(dR_toa(1,1,:)),squeeze(dR_toa(6,1,:)),'Rows','complete');
cctoa(1,1) = temp1(1,2); cctoa(1,1) = roundn(cctoa(1,1),-3);
temp1 = corrcoef(squeeze(dR_toa(1,2,:)),squeeze(dR_toa(6,2,:)),'Rows','complete');
cctoa(2,1) = temp1(1,2); cctoa(2,1) = roundn(cctoa(2,1),-3);
temp1 = corrcoef(squeeze(dR_toa(1,3,:)),squeeze(dR_toa(6,3,:)),'Rows','complete');
cctoa(3,1) = temp1(1,2); cctoa(3,1) = roundn(cctoa(3,1),-3);

for ii = 1:3
    var = char(vector_var_str{ii});
    if strcmp(var,'net')
        spectrum = 1;
        nkk = 7; % number of rows
        titles = {'(a) Total','(b) Air Temperature ','(c) Surface Temperature','(d) Water Vapor','(e) Albedo','(f) Cloud','(g) Residual'};
    elseif strcmp(var,'lw')
        spectrum = 2;
        nkk = 6;
        titles = {'(a) Total','(b) Air Temperature','(c) Surface Temperature','(d) Water Vapor','(e) Cloud','(f) Residual'};
    elseif strcmp(var,'sw')
        spectrum = 3;
        nkk = 5;
        titles = {'(a) Total','(b) Water Vapor','(c) Albedo','(d) Cloud','(e) Residual'};
    end
    
    %% read the data
    dRsfc = dR_sfc(:,spectrum,:);                                              % (x, spectrum, time) all-sky Surface
    dRsfc = squeeze(dRsfc);                                                   % (x, time)
    dRtoa = dR_toa(:,spectrum,:);                                              % (x, spectrum, time) all-sky Surface
    dRtoa = squeeze(dRtoa);                                                   % (x, time)
    %% plot
    figure(ii)
    hold all
    if ii == 1
        num = [1,2,3,4,5,6,7];
    elseif ii == 2
        num = [1,2,3,4,6,7];
    elseif ii == 3
        num = [1,4,5,6,7];
    end
    
    for jj = 1:nkk
         xx = num(1,jj);
         subplot_yih(nkk,2,jj,1);
         set_font_size(8);
         hold all
         plot(time,squeeze(dRsfc(xx,:)),'r')
         grid on
         MAX = max(abs(dRsfc(xx,:))); MAX = max(MAX);
            if MAX > roundn(MAX,1)
               ylim([-roundn(MAX,1)-5 roundn(MAX,1)+5]);
            else
               ylim([-roundn(MAX,1) roundn(MAX,1)]);
            end
            set(gca,'YMinorTick','on'); 
            set(gca,'ticklength',[0.02 0.01]);
            datetick('x','yyyy','keepticks');
            xlim([time(1) time(ntime)])
            set(gca,'XMinorTick','on'); 
            set(gca,'ticklength',[0.02 0.01]);
            set(gca,'ycolor','k')
            title(titles(jj),'FontWeight','bold');
            if jj == 1
                legend({'SFC dR'},'Fontsize',8,'Location','Best')
            elseif jj == nkk-1
                legend({['Cld&Tot cc =',num2str(ccsfc(ii,1))]},'Fontsize',8,'Location','Best')
            end
            if jj ~= nkk
                xticks([732358 734153 735979])
                set(gca,'xticklabel',[]);
            else
                xticks([732358 734153 735979])
                xticklabels({'2005','2010','2015'})
            end
         
         subplot_yih(nkk,2,jj,2);
         set_font_size(8);
         hold all
         plot(time,squeeze(dRtoa(xx,:)),'b')
         grid on
         MAX = max(abs(dRtoa(xx,:))); MAX = max(MAX);
            if MAX > roundn(MAX,1)
               ylim([-roundn(MAX,1)-5 roundn(MAX,1)+5]);
            else
               ylim([-roundn(MAX,1) roundn(MAX,1)]);
            end
            set(gca,'YMinorTick','on'); 
            set(gca,'ticklength',[0.02 0.01]);
            datetick('x','yyyy','keepticks');
            xlim([time(1) time(ntime)])
            set(gca,'XMinorTick','on'); 
            set(gca,'ticklength',[0.02 0.01]);
            set(gca,'ycolor','k')
            title(titles(jj),'FontWeight','bold');
            if jj == 1
                legend({'TOA dR'},'Fontsize',8,'Location','Best')
            elseif jj == nkk-1
                legend({['Cld&Tot cc =',num2str(cctoa(ii,1))]},'Fontsize',8,'Location','Best')
            end
            if jj ~= nkk
                xticks([732358 734153 735979])
                set(gca,'xticklabel',[]);
            else
                xticks([732358 734153 735979])
                xticklabels({'2005','2010','2015'})
            end
            hold off
    end
   
end