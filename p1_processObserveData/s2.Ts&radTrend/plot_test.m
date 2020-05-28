mean_all(:,:,1,1)=mean_RdwnERA_lw_cld;
mean_all(:,:,1,2)=mean_RdwnCE_lw_cld;
mean_all(:,:,2,1)=mean_RnetERA_sw_cld;
mean_all(:,:,2,2)=mean_RnetCE_sw_cld;
title_name{1,1}={'RdwnERA lw'};title_name{1,2}={'RdwnCE lw'};
title_name{2,1}={'RnetERA sw'};title_name{2,2}={'RnetCE sw'};

[cc0,pp0] = corrcoef(mean_RdwnERA_lw_cld,mean_RdwnCE_lw_cld,'Rows','complete');
[cc1,pp1] = corrcoef(mean_RnetERA_sw_cld,mean_RnetCE_sw_cld,'Rows','complete');
cc(1,1) = cc0(1,2);pp(1,1) = pp0(1,2);
cc(2,1) = cc1(1,2);pp(2,1) = pp1(1,2);
%plot to compare
mmin=-1;
mmax=1;

%��ͼǰ����
set(0,'defaultfigurecolor','w')%���û�����ɫΪ��ɫ
% set(gcf,'outerposition',get(0,'screensize'));%����figureȫ��
figure
hold on
[Xlon,Ylat] = meshgrid(lone,late);
lon1=[2.5 357.5];lat1=[-75 75];%������γ�ȷ�Χ
% ��2*2��year�ֲ���ͼ,iΪ��(���������(down, net_short)), jΪ��(��ͬ������ERA, CE),
for i = 1:2
    for j = 1:2
    corr = roundn(cc(i,1),-2);
    subplot_yc(2,2,i,j);% ����ͼ
    hold on
    temp=squeeze(mean_all(:,:,i,j));
    m_proj('Mercator','lon',lon1, 'lat',lat1);%ͶӰ�����Ԥ�����Mercator,Equidistant cylindrical,Miller Cylindrical
    m_contourf(Xlon,Ylat,temp',30,'LineStyle','none')
    colormap(mycolor(16));       %��ת��ɫ��colormap(flipud(mycolor(13)));%colormap(jet(4))%���ٿ����Ӻͼ�С����
    caxis([mmin mmax]);
    hold on
    m_line(world_mapx(:),world_mapy(:),'color',[0 0 0],'LineWidth',0.5);
    m_grid('linestyle','none','tickdir','out','yaxislocation','left','fontsize',10,'color','k');%�߿�������ѡ��
    title([title_name{i,j},['(cc = ',num2str(corr),')']]);% cc=',num2str(corr))
    hold on
    pos = get(gca,'Position');
    end
end
colorbar('Position',[pos(1)+pos(3)+0.005 pos(2) 0.01 pos(4)*1.5])
caxis([mmin mmax])