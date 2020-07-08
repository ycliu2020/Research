% x_plot_month_trend_lyc.m
% plot Ts/rhs trend;trend1:(12,2) k and b, from 1-12 month
% note: units:/10a eg: trend*365*10
% Ts: K/10 a , rhs:W/m^(-2)/10 a
% color map: 16
clc; clear;tic; %��ʱ
add_all;
pre_load;
inputpath = 'E:\Repeat_the_experiment\testdata\ERAi\';
load([inputpath,'trend_tsrad.mat']);% ���ر���: 'lon', 'lat', 'time',cons_rhs, cons_Ts, lat, lon, p_rhs, p_Ts, time, trend_rhs, trend_Ts
load('D:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\tools\��ͼ\map\a_world_map\mask144.mat') % load word land mask
load('D:\onedrive\OneDrive - smail.nju.edu.cn\a_research\database_work\programming\matlab\tools\��ͼ\map\a_world_map\correct_worldmap.mat') % ������������ͼ�ļ�word_mapx(:),word_mapy(:)
[Xlon,Ylat] = meshgrid(lon,lat);
month = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
season = {'MAM','JJA','SON','DJF'};
action = {'Ts trend(K/10a)','LW��+SW trend(W��m-2/10a)'};

% multiply 10a
trend_Ts_m = squeeze(mon_trend_Ts(:,:,:,1))*365*10;
trend_rhs_m = squeeze(mon_trend_rhs(:,:,:,1))*365*10;
%season
trend_Ts_s = sea_trend_Ts(:,:,:)*365*10;
trend_rhs_s = sea_trend_rhs(:,:,:)*365*10;
%year
trend_Ts_y = yr_trend_Ts(:,:)*365*10;
trend_rhs_y = yr_trend_rhs(:,:)*365*10;

%----------------------this is every season and both equation sides(only land)----------------------------------
% trend = trend_Ts_m; % decide varible to plot
trend1 = trend_Ts_s; % decide varible to plot
trend2 = trend_rhs_s; % decide varible to plot
for ii=1:4
temp1=squeeze(trend1(:,:,ii));
temp2=squeeze(trend2(:,:,ii));
temp1(~maskworld)=NaN;%mask land
temp2(~maskworld)=NaN;%mask land
trend1(:,:,ii)=temp1;
trend2(:,:,ii)=temp2;
end
trend1(:,lat>60&lat<-60,:)=NaN;%����γ���޳���������
trend2(:,lat>60&lat<-60,:)=NaN;
cc = zeros(4,1); pp = zeros(4,1);
% Compute correlation(seasonly)
for ii = 1:4
    [cc0,pp0] = corrcoef(squeeze(trend1(:,:,ii)),squeeze(trend2(:,:,ii)),'Rows','complete');
    cc(ii,1) = -cc0(1,2);
    pp(ii,1) = -pp0(1,2);
end

lon1=[2.5 357.5];lat1=[-60 60];%������γ�ȷ�Χ


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
% ��4*2��season�ֲ���ͼ, ��Ϊ��ͬseasons,��Ϊ��������, iΪ��, jΪ��
for ii = 1:8
    if ii<5
        i=ii;j=1;kk=1;
        trendm = squeeze(trend1(:,:,ii));   
        corr = roundn(cc(ii,1),-2);% ���ϵ��   
    elseif ii>=5
        i=ii-4;j=2;kk=2;
        trendm = squeeze(trend2(:,:,ii-4));
    end

    subplot_yih(4,2,i,j);% ����ͼ
    hold on
    
    m_proj('Mercator','lon',lon1, 'lat',lat1);%ͶӰ�����Ԥ�����Mercator,Equidistant cylindrical,Miller Cylindrical
    m_contourf(Xlon,Ylat,trendm',30,'LineStyle','none')
    colormap(mycolor(16));       %��ת��ɫ��colormap(flipud(mycolor(13)));%colormap(jet(4))%���ٿ����Ӻͼ�С����
    caxis([mmin mmax]);
    hold on
    m_line(world_mapx(:),world_mapy(:),'color',[0 0 0],'LineWidth',0.5);
    
    
    m_grid('linestyle','none','tickdir','out','yaxislocation','left','fontsize',10,'color','k');%�߿�������ѡ��
    if ii<5
        title(strcat(season(ii),'  cc=',num2str(corr),action(kk)));% cc=',num2str(corr))
    elseif ii>=5
        title(strcat(season(ii-4),action(kk)));% cc=',num2str(corr))
    end
    hold on
    pos = get(gca,'Position');
    
end
colorbar('Position',[pos(1)+pos(3)+0.005 pos(2) 0.01 pos(4)*3.73])
caxis([mmin mmax])


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
% fig_name='C:\Users\LYC\Desktop\��ǰ����ӿ�\222';
% export_fig( gcf , '-png','-r500'  , fig_name );

%----------------------plot figure about every season and both equation sides ----------------------------------

t=toc;
disp(t)