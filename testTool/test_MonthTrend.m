%%---------------------------------------------------------
% Author       : LYC
% Date         : 2020-06-09 15:52:00
% LastEditTime : 2020-06-14 16:47:43
% LastEditors  : LYC
% Description  : 
% FilePath     : /Research/testTool/test_MonthTrend.m
%  
%%---------------------------------------------------------
% aim to test the differnce between detrend month and detrend year
clc; clear;
% % lab1
% time=linspace(1,24,24);
% ntime=length(time);

% vars=20*sin(time*pi/6+pi/6);

% figure;
% % plot(time,vars)
% p=polyfit(time,vars,1);
% f = polyval(p,time);
% plot( time,vars,time, f, 'o')
% % method1
% [~, trendm1, cons_m, p_m] = detrend_yan(vars, time);

% trendyr1 = squeeze(mean(trendm1(:,1),1));

% % method2
% t2=1:ntime/12;
% for im = 1:length(t2)
%     dy_yr(im)=nanmean(vars((im-1)*12+1:im*12));
% end
% [testx, trendm2, cons_m2, p_m2] = detrendyr_yc(dy_yr, t2);
% trendyr2 = trendm2(:,1);
% lab2
filename2 = '/data1/liuyincheng/Observe-rawdata/HadCRUT4/HadCRUT.4.6.0.0.median.nc';
dts0 = ncread(filename2,'temperature_anomaly');
dts0_had = dts0(:,:,1803:2018);    % 1803:2000.03, 2018:2018.02, 72*36*216 surface temperature anomaly, Global land
lonh = ncread(filename2,'longitude'); nlonh = length(lonh);                  % -177.5:5:177.5
lath = ncread(filename2,'latitude'); nlath = length(lath);
timeh = ncread(filename2,'time');  timeh = timeh(1803:2018);                 % 1979/03 - 2018/02
ntimeh = length(timeh);
timeh = timeh + datenum(1850,1,1);
timeh=double(timeh);
%缺测太多将在时间序列中只要缺测一次就剔除
test=dts0_had;
test1=sum(isnan(dts0_had),3);
for ii=1:216
    temp1=test(:,:,ii);
   temp1(test1~=0)=0;
   test(:,:,ii)=temp1;
end
dts0_had=test;
% method 1
trendm_dTsHad = zeros(nlonh, nlath, 12, 2);  % trend, 12 months
p_dTsHad = zeros(nlonh, nlath, 12);  % f test
cons_dTsHad = zeros(nlonh, nlath, 12);  % mean
% months mean trend(unite:per day)
for i = 1:nlonh
    for j = 1:nlath
        [~, trendm_dTsHad(i, j, :, :), cons_dTsHad(i, j, :), p_dTsHad(i, j, :)] = detrend_yan(dts0_had(i, j, :), timeh);
    end
end

% years mean trend (unite:per day)
trendyr1_dTsHad = squeeze(mean(trendm_dTsHad(:,:,:,1),3));

% method 2
time2=1:ntimeh/12;
for im = 1:length(time2)
    dts0_hadyr(:,:,im)=nanmean(dts0_had(:,:,(im-1)*12+1:im*12),3);
end

for i = 1:nlonh
    for j = 1:nlath
        [testx, trendm_dTsHad2(i, j, :), cons_dTsHad2(i, j), p_dTsHad2(i, j)] = detrendyr_yc(dts0_hadyr(i,j,:), time2);
    end
end
trendyr2_dTsHad=trendm_dTsHad2(:,:,1);
