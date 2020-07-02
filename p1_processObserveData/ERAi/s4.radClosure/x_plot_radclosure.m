% This program is to test radiation closure of kernel method
% For four regions: China/ 20-50N/ Northern Hemisphere land/ Global land
%����ʱ������ͼ
clc;clear;
pre_load;
inputpath='E:\Repeat_the_experiment\testdata\radclosure_date\';
load([inputpath,'wordmean_r.mat'])% meanarea,late,lone
name_num=1;
name={'Radiation closure experiment (China mean)','Radiation closure experiment (60N-60S mean)'};
ntime = length(time);
% cal the coef
cc=zeros(2,3);
for i = 1:2 % sfc and toa
   for j = 1:3 % total, short,long
      temp1=meanarea(:,i,j,1);
      temp2=meanarea(:,i,j,2);
      cc0 = corrcoef(temp1,temp2,'Rows','complete');
      cc(i,j) = cc0(1,2);
      cc(i,j) = roundn(cc(i,j),-3);
   end
end


level={'sfc','toa'};Level=upper(level);
Radname={'dR_{NET} ','dR_{LW} ','dR_{SW} '};
set(0,'defaultfigurecolor','w')%���û�����ɫΪ��ɫ
figure
for i = 1:2 % sfc and toa
   for j = 1:3 % total, short,long
      subplot_yc(3,2,j,i)
      hold on
      plot(time,squeeze(meanarea(:,i,j,1)),'b')
      plot(time,squeeze(meanarea(:,i,j,2)),'r')
      plot(time,squeeze(meanarea(:,i,j,3)),'k')
      %����������
      datetick('x','yyyy','keepticks');% �����ڵĸ�ʽ��ʾ������
      xlim([time(1) time(ntime)])
      if name_num==1
           ylim([-10 10])
      else
           ylim([-5 5])
      end

      xticks([732358 734153 735979])
      xticklabels({'2005','2010','2015'})
      ax=gca;
      ax.XMinorTick = 'on';ax.YMinorTick = 'on'; % �����ο̶���
      ax.TickLength = [0.02 0.01];  %�̶��߳���      set(gca,'ticklength',[0.02 0.01]);
      ax.XColor = 'k';ax.YColor = 'k'; % ���ÿ̶�����ɫ

      ylabel({Level{i},' Rad Anomaly (Wm^{-2})'},'Fontsize',10)
      legend({Radname{j}, ['Diagnosed cc =',num2str(cc(i,j))], 'Residual'},'Fontsize',8,'Location', 'northwest')
      legend('boxoff') %ɾ��ͼ������������
   end
end
sgtitle(name{name_num})



fig_name='E:\Repeat_the_experiment\a_figure\����պ�\dp';
export_fig( gcf , '-png','-r500'  , fig_name );














