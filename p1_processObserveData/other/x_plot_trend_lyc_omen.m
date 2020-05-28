% x_plot_month_trend_lyc.m
% plot Ts/rhs trend;ts1:(12,2) k and b, from 1-12 month
% note: units:/10a eg: trend*365*10
% Ts: K/10 a , rhs:W/m^(-2)/10 a
% color map: 16
% contain three database, two resolution(lon,lat, lath,lonh)
%
clc; clear;tic; %��ʱ
add_all;
pre_load;
inputpath_ERA = 'G:\database\homework\ERAi\';
inputpath_Had_CER = 'G:\database\homework\CERES_EBAF\';
load([inputpath_ERA,'trend_tsrad.mat']);% ���ر���: 'lon', 'lat', 'time',cons_rhs, cons_Ts, lat, lon, p_rhs, p_Ts, time, trend_rhs, trend_Ts
load([inputpath_Had_CER,'trend_CERESrad_HadTs.mat']);% ���ر���: mon_trend_CErhs, mon_trend_hadTs, cons_CErhs, cons_hadTs, late, lath, lone, lonh, p_CErhs, p_hadTs, sea_trend_CErhs, sea_trend_hadTs, time, timeh, yr_trend_CErhs, yr_trend_hadTs
load('G:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\tools\plot\map\a_world_map\mask\mask_cp144.mat') % load word land mask
load('G:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\tools\plot\map\a_world_map\mask\mask_ce72.mat') % load word land mask
load('G:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\tools\plot\map\a_world_map\correct_worldmap.mat') % ������������ͼ�ļ�word_mapx(:),word_mapy(:)
load('G:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\tools\plot\map\a_china_map\mask14472.mat')


month = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
season = {'MAM','JJA','SON','DJF'};
action = {'Ts trend(K/10a)','LW��+SW trend(W��m-2/10a)'};
dataname = {'ERAi','CERES','HadCRUT4'};

% multiply 10a
trend_Ts_m = squeeze(mon_trend_Ts(:,:,:,1))*365*10;
trend_rhs_m = squeeze(mon_trend_rhs(:,:,:,1))*365*10;
trend_hadTs_m = squeeze(mon_trend_hadTs(:,:,:,1))*365*10;
trend_CErhs_m = squeeze(mon_trend_CErhs(:,:,:,1))*365*10;
%season
trend_Ts_s = sea_trend_Ts(:,:,:)*365*10;
trend_rhs_s = sea_trend_rhs(:,:,:)*365*10;
trend_hadTs_s = squeeze(sea_trend_hadTs(:,:,:))*365*10;
trend_CErhs_s = squeeze(sea_trend_CErhs(:,:,:))*365*10;
%year
trend_Ts_y = yr_trend_Ts(:,:)*365*10;
trend_rhs_y = yr_trend_rhs(:,:)*365*10;
trend_hadTs_y = squeeze(yr_trend_hadTs(:,:))*365*10;
trend_CErhs_y = squeeze(yr_trend_CErhs(:,:))*365*10;

%----------------------section5:yearͼ(������had)----------------------------------
ts1 = trend_Ts_y; % decide varible to plot
rhs1 = trend_rhs_y; % decide varible to plot
ts2 = trend_hadTs_y; % decide varible to plot
rhs2 = trend_CErhs_y; % decide varible to plot



% ts1 = trend_Ts_y; % decide varible to plot
% rhs1 = trend_rhs_y; % decide varible to plot

ts1(~maskworld_cp)=NaN;%mask land
rhs1(~maskworld_cp)=NaN;%mask land
rhs2(~maskworld_cp)=NaN;


% ts1(:,lat>60|lat<-60,:)=NaN;%����γ���޳���������
% rhs1(:,lat>60|lat<-60,:)=NaN;
cc = zeros(2,1); pp = zeros(2,1);
% pos = zeros(4,2);
% Compute correlation(seasonly)
[cc0,pp0] = corrcoef(ts1(:,:),rhs1(:,:),'Rows','complete');
[cc1,pp1] = corrcoef(ts1(:,:),rhs2(:,:),'Rows','complete');
cc(1,1) = cc0(1,2);pp(1,1) = pp0(1,2);
cc(2,1) = cc1(1,2);pp(2,1) = pp1(1,2);
lon1=[2.5 357.5];lat1=[-75 75];%������γ�ȷ�Χ

%----------------------plot figure about every season and both equation sides(only land)----------------------------------
% MAX = squeeze(max(max(max(abs(trend)))));
% mmin = -roundn(MAX,0); mmax = roundn(MAX,0);
% mmin=mmin';mmax=mmax';
mmin=-1;
mmax=1;
%��ͼǰ����
% set(0,'defaultfigurecolor','w')%���û�����ɫΪ��ɫ
% set(gcf,'outerposition',get(0,'screensize'));%����figureȫ��
figure(1)
hold on
[Xlon,Ylat] = meshgrid(lon,lat);
% ��2*2��year�ֲ���ͼ,iΪ��, jΪ��,kΪts or rhs,p1Ϊdataname
for ii = 1:3
    if ii==1
        i=ii;j=1;kk=1;p1=1;
        trendm = ts1(:,:);
    elseif ii==2
        i=1;j=2;kk=2;p1=1;
        trendm = rhs1(:,:);
        corr = roundn(cc(1,1),-2);% ���ϵ��
    elseif ii==3
        i=2;j=2;kk=2;p1=2;
        trendm = rhs2(:,:);
        corr = roundn(cc(2,1),-2);% ���ϵ��
    end
    
    subplot_yih(2,2,i,j);% ����ͼ
    hold on
    
    m_proj('Mercator','lon',lon1, 'lat',lat1);%ͶӰ�����Ԥ�����Mercator,Equidistant cylindrical,Miller Cylindrical
    m_contourf(Xlon,Ylat,trendm',30,'LineStyle','none')
    colormap(mycolor(16));       %��ת��ɫ��colormap(flipud(mycolor(13)));%colormap(jet(4))%���ٿ����Ӻͼ�С����
    caxis([mmin mmax]);
    hold on
    m_line(world_mapx(:),world_mapy(:),'color',[0 0 0],'LineWidth',0.5);
    m_grid('linestyle','none','tickdir','out','yaxislocation','left','fontsize',10,'color','k');%�߿�������ѡ��
    title({[dataname{p1},' ',action{kk}];season{i}});% cc=',num2str(corr))
    
    if ii==1
        title({[dataname{p1},' ',action{kk}];['year mean ',]});% cc=',num2str(corr))
    else
        title({[dataname{p1},' ',action{kk}];['year mean (cc = ',num2str(corr),')']});% cc=',num2str(corr))
    end
    hold on
    pos = get(gca,'Position');
end
colorbar('Position',[pos(1)+pos(3)+0.005 pos(2) 0.01 pos(4)*1.5])
caxis([mmin mmax])

% %----------------------section4: had yearͼ,ר��Ϊhad�������----------------------------------
% ts1 = trend_Ts_y; % decide varible to plot
% rhs1 = trend_rhs_y; % decide varible to plot
% ts2 = trend_hadTs_y; % decide varible to plot
% rhs2 = trend_CErhs_y; % decide varible to plot

% ts2(~maskworld_ce7236)=NaN;
% ts2(ts2==0)=NaN;

% % note ts2 ���겻��centre ̫ƽ��, ��������ת��
% tpt=ts2;
% ts2(1:36,:)=tpt(37:72,:);
% ts2(37:72,:)=tpt(1:36,:);
% lonh=lonh+180;
% maskhadts=isnan(ts2(:,:)); %nan��1, ��Ϊhadts����Ĥ

% %%������������ֵ��had������
lonf = 0:2.5:357.5; nlonf = length(lonf);
latf = 88.75:-2.5:-88.75; nlatf = length(latf);
lath = double(lath);
lonh = double(lonh);

[Xlonh, Ylath] = meshgrid(lath,lonh);
[Xlonf, Ylatf] = meshgrid(latf,lonf);
% % ��Ϊ��ֵ����
% ts11 = interp2(Xlonf,Ylatf,ts1,Xlonh,Ylath);
% rhs11 = interp2(Xlonf,Ylatf,rhs1,Xlonh,Ylath);
% rhs21 = interp2(Xlonf,Ylatf,rhs2,Xlonh,Ylath);

% ts11(maskhadts)=NaN;%mask land
% rhs11(maskhadts)=NaN;%mask land
% rhs21(maskhadts)=NaN;%mask land


% % ts1(:,lat>60|lat<-60,:)=NaN;%����γ���޳���������
% % rhs1(:,lat>60|lat<-60,:)=NaN;
% cc = zeros(3,1); pp = zeros(3,1);
% % pos = zeros(4,2);
% % Compute correlation(seasonly)

%     [cc0,pp0] = corrcoef(ts2(:,:),ts11(:,:),'Rows','complete');
%     [cc1,pp1] = corrcoef(ts2(:,:),rhs11(:,:),'Rows','complete');
%     [cc2,pp2] = corrcoef(ts2(:,:),rhs21(:,:),'Rows','complete');
%     cc(1,1) = cc0(1,2);pp(1,1) = pp0(1,2);
%     cc(2,1) = cc1(1,2);pp(2,1) = pp1(1,2);
%     cc(3,1) = cc2(1,2);pp(3,1) = pp2(1,2);


% lon1=[2.5 357.5];lat1=[-75 75];%������γ�ȷ�Χ

% %----------------------plot figure about every season and both equation sides(only land)----------------------------------
% % MAX = squeeze(max(max(max(abs(trend)))));
% % mmin = -roundn(MAX,0); mmax = roundn(MAX,0);
% % mmin=mmin';mmax=mmax';
% mmin=-1;
% mmax=1;
% %��ͼǰ����
% % set(0,'defaultfigurecolor','w')%���û�����ɫΪ��ɫ
% % set(gcf,'outerposition',get(0,'screensize'));%����figureȫ��
% figure(1)
% hold on
% [Xlon,Ylat] = meshgrid(lonh,lath);
% % ��2*2��year�ֲ���ͼ,iΪ��, jΪ��,kΪts or rhs,p1Ϊdataname
% for ii = 1:4
%     if ii==1
%         i=ii;j=1;kk=1;p1=3;
%         trendm = ts2(:,:);
%     elseif ii==2
%         i=1;j=2;kk=1;p1=1;
%         trendm = ts11(:,:);
%         corr = roundn(cc(1,1),-2);% ���ϵ��
%     elseif ii==3
%         i=2;j=1;kk=2;p1=1;
%         trendm = rhs11(:,:);
%         corr = roundn(cc(2,1),-2);% ���ϵ��
%     elseif ii==4
%         i=2;j=2;kk=2;p1=2;
%         trendm = rhs21(:,:);
%         corr = roundn(cc(3,1),-2);% ���ϵ��
%     end

%     subplot_yih(2,2,i,j);% ����ͼ
%     hold on

%     m_proj('Mercator','lon',lon1, 'lat',lat1);%ͶӰ�����Ԥ�����Mercator,Equidistant cylindrical,Miller Cylindrical
%     m_contourf(Xlon,Ylat,trendm',30,'LineStyle','none')
%     colormap(mycolor(16));       %��ת��ɫ��colormap(flipud(mycolor(13)));%colormap(jet(4))%���ٿ����Ӻͼ�С����
%     caxis([mmin mmax]);
%     hold on
%     m_line(world_mapx(:),world_mapy(:),'color',[0 0 0],'LineWidth',0.5);
%     m_grid('linestyle','none','tickdir','out','yaxislocation','left','fontsize',10,'color','k');%�߿�������ѡ��
%     title({[dataname{p1},' ',action{kk}];season{i}});% cc=',num2str(corr))

%     if ii==1
%             title({[dataname{p1},' ',action{kk}];['year mean ',]});% cc=',num2str(corr))
%     else
%             title({[dataname{p1},' ',action{kk}];['year mean (cc = ',num2str(corr),')']});% cc=',num2str(corr))
%     end
%     hold on
%     pos = get(gca,'Position');
% end
% colorbar('Position',[pos(1)+pos(3)+0.005 pos(2) 0.01 pos(4)*1.5])
% caxis([mmin mmax])

% %----------------------section3: had seasonͼ,ר��Ϊhad�������----------------------------------
% ts1 = trend_Ts_s; % decide varible to plot
% rhs1 = trend_rhs_s;
% ts2 = trend_hadTs_s;
% rhs2 = trend_CErhs_s;

% for ii=1:4
%     temp4=squeeze(ts2(:,:,ii));
%     temp4(~maskworld_ce7236)=NaN;
%     ts2(:,:,ii)=temp4;
% end
% ts2(ts2==0)=NaN;
% % note ts2 ���겻��centre ̫ƽ��, ��������ת��
% tpt=ts2;
% ts2(1:36,:,:)=tpt(37:72,:,:);
% ts2(37:72,:,:)=tpt(1:36,:,:);
% lonh=lonh+180;
% maskhadts=isnan(ts2(:,:,1)); %nan��1, ��Ϊhadts����Ĥ

% %%������������ֵ��had������
% lonf = 0:2.5:357.5; nlonf = length(lonf);
% latf = 88.75:-2.5:-88.75; nlatf = length(latf);
% lath = double(lath);
% lonh = double(lonh);
% time2 = 1:4;

% [Xlonh, Ylath, Ttimeh] = meshgrid(lath,lonh,time2);
% [Xlonf, Ylatf, Ttimef] = meshgrid(latf,lonf,time2);
% % ��Ϊ��ֵ����
% ts11 = interp3(Xlonf,Ylatf,Ttimef,ts1,Xlonh,Ylath,Ttimeh);
% rhs11 = interp3(Xlonf,Ylatf,Ttimef,rhs1,Xlonh,Ylath,Ttimeh);
% rhs21 = interp3(Xlonf,Ylatf,Ttimef,rhs2,Xlonh,Ylath,Ttimeh);
% for ii=1:4
%     temp1=squeeze(ts11(:,:,ii));%ERA ts
%     temp2=squeeze(rhs11(:,:,ii));%ERA rad
%     temp3=squeeze(rhs21(:,:,ii));% CERES rad
%     temp1(maskhadts)=NaN;%mask land
%     temp2(maskhadts)=NaN;%mask land
%     temp3(maskhadts)=NaN;%mask land
%     ts11(:,:,ii)=temp1;
%     rhs11(:,:,ii)=temp2;
%     rhs21(:,:,ii)=temp3;
% end

% % ts1(:,lat>60|lat<-60,:)=NaN;%����γ���޳���������
% % rhs1(:,lat>60|lat<-60,:)=NaN;
% % rhs2(:,lat>60|lat<-60,:)=NaN;
% lon1=[2.5 357.5];lat1=[-75 75];%������γ�ȷ�Χ

% cc = zeros(4,3); pp = zeros(4,3);
% % Compute correlation(seasonly)
% for ii = 1:4
%     [cc0,pp0] = corrcoef(squeeze(ts2(:,:,ii)),squeeze(rhs11(:,:,ii)),'Rows','complete');
%     [cc1,pp1] = corrcoef(squeeze(ts2(:,:,ii)),squeeze(rhs21(:,:,ii)),'Rows','complete');
%     [cc2,pp2] = corrcoef(squeeze(ts2(:,:,ii)),squeeze(ts11(:,:,ii)),'Rows','complete');
%     cc(ii,1) = cc0(1,2);pp(ii,1) = pp0(1,2);
%     cc(ii,2) = cc1(1,2);pp(ii,2) = pp1(1,2);
%     cc(ii,3) = cc2(1,2);pp(ii,3) = pp2(1,2);
% end


% %----------------------plot figure about every season and both equation sides(only land)----------------------------------
% % MAX = squeeze(max(max(max(abs(trend)))));
% % mmin = -roundn(MAX,0); mmax = roundn(MAX,0);
% % mmin=mmin';mmax=mmax';
% mmin=-1;
% mmax=1;
% %��ͼǰ����
% set(0,'defaultfigurecolor','w')%���û�����ɫΪ��ɫ
% set(gcf,'outerposition',get(0,'screensize'));%����figureȫ��
% figure(1)
% hold on
% [Xlon,Ylat] = meshgrid(lonh,lath);
% % ��4*3��season�ֲ���ͼ, ��Ϊ��ͬseasons,��Ϊ��ͬ���Ϸ�������, iΪ��, jΪ��,kΪts or rhs,p1Ϊdataname
% for ii = 1:16
%     if ii<5
%         i=ii;j=1;kk=1;p1=3;
%         trendm = squeeze(ts2(:,:,i));
%     elseif ii>=5&&ii<9
%         i=ii-4;j=2;kk=1;p1=1;
%         trendm = squeeze(ts11(:,:,i));
%         corr = roundn(cc(i,3),-2);% ���ϵ��
%     elseif ii>=9&&ii<13
%         i=ii-8;j=3;kk=2;p1=1;
%         trendm = squeeze(rhs11(:,:,i));
%         corr = roundn(cc(i,1),-2);% ���ϵ��
%     elseif ii>=13
%         i=ii-12;j=4;kk=2;p1=2;
%         trendm = squeeze(rhs21(:,:,i));
%         corr = roundn(cc(i,2),-2);% ���ϵ��
%     end

%     subplot_yih(4,4,i,j);% ����ͼ
%     hold on
%     m_proj('Mercator','lon',lon1, 'lat',lat1);%ͶӰ�����Ԥ�����Mercator,Equidistant cylindrical,Miller Cylindrical
%     m_contourf(Xlon,Ylat,trendm',15,'LineStyle','none')
%     colormap(mycolor(16));       %��ת��ɫ��colormap(flipud(mycolor(13)));%colormap(jet(4))%���ٿ����Ӻͼ�С����
%     caxis([mmin mmax]);
%     hold on
%     m_line(world_mapx(:),world_mapy(:),'color',[0 0 0],'LineWidth',0.5);


%     m_grid('linestyle','none','tickdir','out','yaxislocation','left','fontsize',10,'color','k');%�߿�������ѡ��
%     if ii<5
%         if ii==1
%             title({[dataname{p1},' ',action{kk}];season{i}});% cc=',num2str(corr))
%         else
%             title(season{i});% cc=',num2str(corr))
%         end
%     elseif ii>=5
%         if ii==5||ii==9||ii==13
%             title({[dataname{p1},' ',action{kk}];[season{i},' (cc = ',num2str(corr),')']});% cc=',num2str(corr))
%         else
%             title([season{i},' (cc = ',num2str(corr),')']);% cc=',num2str(corr))
%         end
%     end
%     hold on
%     pos = get(gca,'Position');

% end
% colorbar('Position',[pos(1)+pos(3)+0.005 pos(2) 0.01 pos(4)*3.73])
% caxis([mmin mmax])


% %----------------------section2: this is every season and both equation sides(only land)----------------------------------
% ts1 = trend_Ts_s; % decide varible to plot
% rhs1 = trend_rhs_s;
% ts2 = trend_hadTs_s;
% rhs2 = trend_CErhs_s;

% for ii=1:4
%     temp1=squeeze(ts1(:,:,ii));%ERA ts
%     temp2=squeeze(rhs1(:,:,ii));%ERA rad
%     temp3=squeeze(rhs2(:,:,ii));% CERES rad
%     temp1(~maskworld_cp)=NaN;%mask land
%     temp2(~maskworld_cp)=NaN;%mask land
%     temp3(~maskworld_cp)=NaN;%mask land
%     ts1(:,:,ii)=temp1;
%     rhs1(:,:,ii)=temp2;
%     rhs2(:,:,ii)=temp3;
% end

% % ts1(:,lat>60|lat<-60,:)=NaN;%����γ���޳���������
% % rhs1(:,lat>60|lat<-60,:)=NaN;
% % rhs2(:,lat>60|lat<-60,:)=NaN;
% lon1=[2.5 357.5];lat1=[-75 75];%������γ�ȷ�Χ

% cc = zeros(4,2); pp = zeros(4,2);
% % Compute correlation(seasonly)
% for ii = 1:4
%     [cc0,pp0] = corrcoef(squeeze(ts1(:,:,ii)),squeeze(rhs1(:,:,ii)),'Rows','complete');
%     [cc1,pp1] = corrcoef(squeeze(ts1(:,:,ii)),squeeze(rhs2(:,:,ii)),'Rows','complete');
%     cc(ii,1) = cc0(1,2);pp(ii,1) = pp0(1,2);
%     cc(ii,2) = cc1(1,2);pp(ii,2) = pp1(1,2);
% end


% %----------------------plot figure about every season and both equation sides(only land)----------------------------------
% % MAX = squeeze(max(max(max(abs(trend)))));
% % mmin = -roundn(MAX,0); mmax = roundn(MAX,0);
% % mmin=mmin';mmax=mmax';
% mmin=-1;
% mmax=1;
% %��ͼǰ����
% set(0,'defaultfigurecolor','w')%���û�����ɫΪ��ɫ
% set(gcf,'outerposition',get(0,'screensize'));%����figureȫ��
% figure(1)
% hold on
% % ��4*3��season�ֲ���ͼ, ��Ϊ��ͬseasons,��Ϊ��ͬ���Ϸ�������, iΪ��, jΪ��,kΪts or rhs,p1Ϊdataname
% for ii = 1:12
%     if ii<5
%         i=ii;j=1;kk=1;p1=1;
%         trendm = squeeze(ts1(:,:,i));
%     elseif ii>=5&&ii<9
%         i=ii-4;j=2;kk=2;p1=1;
%         trendm = squeeze(rhs1(:,:,i));
%         corr = roundn(cc(i,1),-2);% ���ϵ��
%     elseif ii>=9
%         i=ii-8;j=3;kk=2;p1=2;
%         trendm = squeeze(rhs2(:,:,i));
%         corr = roundn(cc(i,2),-2);% ���ϵ��
%     end

%     subplot_yih(4,3,i,j);% ����ͼ
%     hold on

%     m_proj('Mercator','lon',lon1, 'lat',lat1);%ͶӰ�����Ԥ�����Mercator,Equidistant cylindrical,Miller Cylindrical
%     m_contourf(Xlon,Ylat,trendm',30,'LineStyle','none')
%     colormap(mycolor(16));       %��ת��ɫ��colormap(flipud(mycolor(13)));%colormap(jet(4))%���ٿ����Ӻͼ�С����
%     caxis([mmin mmax]);
%     hold on
%     m_line(world_mapx(:),world_mapy(:),'color',[0 0 0],'LineWidth',0.5);


%     m_grid('linestyle','none','tickdir','out','yaxislocation','left','fontsize',10,'color','k');%�߿�������ѡ��
%     if ii<5
%         if ii==1
%             title({[dataname{p1},' ',action{kk}];season{i}});% cc=',num2str(corr))
%         else
%             title(season{i});% cc=',num2str(corr))
%         end
%     elseif ii>=5
%         if ii==5||ii==9
%             title({[dataname{p1},' ',action{kk}];[season{i},' (cc = ',num2str(corr),')']});% cc=',num2str(corr))
%         else
%             title([season{i},' (cc = ',num2str(corr),')']);% cc=',num2str(corr))
%         end
%     end
%     hold on
%     pos = get(gca,'Position');

% end
% colorbar('Position',[pos(1)+pos(3)+0.005 pos(2) 0.01 pos(4)*3.73])
% caxis([mmin mmax])

% %----------------------section1: this is yearly and both equation sides(only land)----------------------------------
% ts1 = trend_Ts_y; % decide varible to plot
% rhs1 = trend_rhs_y; % decide varible to plot

% ts1(~maskworld_cp)=NaN;%mask land
% rhs1(~maskworld_cp)=NaN;%mask land

% ts1(:,lat>60|lat<-60,:)=NaN;%����γ���޳���������
% rhs1(:,lat>60|lat<-60,:)=NaN;
% cc = zeros(1,1); pp = zeros(1,1);
% pos = zeros(4,2);
% % Compute correlation(seasonly)

%     [cc0,pp0] = corrcoef(ts1(:,:),rhs1(:,:),'Rows','complete');
%     cc(1,1) = cc0(1,2);
%     pp(1,1) = pp0(1,2);


% lon1=[2.5 357.5];lat1=[-60 60];%������γ�ȷ�Χ

% %----------------------plot figure about every year and both equation sides(only land)----------------------------------
% % MAX = squeeze(max(max(max(abs(trend)))));
% % mmin = -roundn(MAX,0); mmax = roundn(MAX,0);
% % mmin=mmin';mmax=mmax';
% mmin=-1;
% mmax=1;
% %��ͼǰ����
% % set(0,'defaultfigurecolor','w')%���û�����ɫΪ��ɫ
% % set(gcf,'outerposition',get(0,'screensize'));%����figureȫ��
% figure(1)
% hold on
% % ��4*2��season�ֲ���ͼ, ��Ϊ��ͬseasons,��Ϊ��������, iΪ��, jΪ��
% for ii = 1:2
%     if ii==1
%         i=ii;j=1;kk=1;
%         trendm = ts1(:,:);
%         corr = roundn(cc(1,1),-2);% ���ϵ��
%     elseif ii==2
%         i=1;j=2;kk=2;
%         trendm = rhs1(:,:);
%     end

%     subplot_yih(1,2,i,j);% ����ͼ
%     hold on

%     m_proj('Mercator','lon',lon1, 'lat',lat1);%ͶӰ�����Ԥ�����Mercator,Equidistant cylindrical,Miller Cylindrical
%     m_contourf(Xlon,Ylat,trendm',30,'LineStyle','none')
%     colormap(mycolor(16));       %��ת��ɫ��colormap(flipud(mycolor(13)));%colormap(jet(4))%���ٿ����Ӻͼ�С����
%     caxis([mmin mmax]);
%     hold on
%     m_line(world_mapx(:),world_mapy(:),'color',[0 0 0],'LineWidth',0.5);
%     m_grid('linestyle','none','tickdir','out','yaxislocation','left','fontsize',10,'color','k');%�߿�������ѡ��
%     if ii==1
%             title({action{kk};['year mean (cc = ',num2str(corr),')']});% cc=',num2str(corr))
%     elseif ii==2
%             title({action{kk};['year mean']});% cc=',num2str(corr))
%     end
%     hold on
%     pos(:,ii) = get(gca,'Position');

% end
% colorbar('southoutside','Position',[pos(1,1) pos(2,1)-0.06 pos(1,2)+pos(3,2)-0.06 0.03])
% caxis([mmin mmax])


% %----------------------plot every month figure----------------------------------
% % MAX = squeeze(max(max(max(abs(trend)))));
% % mmin = -roundn(MAX,0); mmax = roundn(MAX,0);
% % mmin=mmin';mmax=mmax';
% mmin=-10;
% mmax=10;
% %��ͼǰ����
% set(0,'defaultfigurecolor','w')%���û�����ɫΪ��ɫ
% set(gcf,'outerposition',get(0,'screensize'));%����figureȫ��
% figure(1)
% hold on
% % 1�¿�ʼ, ����ÿ����������ͼλ��, ��4*3���·ֲ���ͼ, ��Ϊ������,��Ϊ��ͬ����, iΪ��, jΪ��
% for ii = 1:12
%     if ii<11
%         mm = ii+2;
%     else
%         mm = ii-10;
%     end
%     if ii<4
%        i = 1;
%        j = ii;
%     elseif ii<7
%        i = 2;
%        j = ii-3;
%     elseif ii<10
%        i = 3;
%        j = ii - 6;
%     else
%        i = 4;
%        j = ii - 9;
%     end
%     trendm = squeeze(trend(:,:,mm));
%     subplot_yih(4,3,i,j);% ����ͼ
%     hold on
%     m_proj('Mercator','lon',lon1, 'lat',lat1);%ͶӰ�����Ԥ�����Mercator,Equidistant cylindrical,Miller Cylindrical
%     m_contourf(Xlon,Ylat,trendm',30,'LineStyle','none')
%     colormap(mycolor(16));       %��ת��ɫ��colormap(flipud(mycolor(13)));%colormap(jet(4))%���ٿ����Ӻͼ�С����
%     caxis([mmin mmax]);
%     hold on
%     m_line(word_mapx(:),word_mapy(:),'color',[0 0 0],'LineWidth',0.5);


%     m_grid('linestyle','none','tickdir','out','yaxislocation','left','fontsize',10,'color','k');%�߿�������ѡ��
%     title(strcat(month(mm)));% cc=',num2str(corr))
%     hold on
%     pos = get(gca,'Position');

% end
% colorbar('Position',[pos(1)+pos(3)+0.005 pos(2) 0.01 pos(4)*3.73])
% caxis([mmin mmax])



% ------output------
% name=['Ta_',t_name{season+1},'_',plv_name{level},'_',num2str(year(1)),'-',num2str(year(2))];
fig_name='C:\Users\LYC\Desktop\11';
export_fig( gcf , '-png','-r500'  , fig_name );

%----------------------plot figure about every season and both equation sides ----------------------------------

t=toc;
disp(t)