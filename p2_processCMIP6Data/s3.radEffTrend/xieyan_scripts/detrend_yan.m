% function [x, tmp_data] = detrend_yan(x,t)
 function [x, trend1, cons, pvalue] = detrend_yan(x,t)

if ~exist('t','var')
  t = 1:length(x);
end

if length(x) ~= length(t)
  disp('ERROR: mismatching x and t!');
  return;
else
  t = reshape(t,size(x));
end

% trend
ntime = length(t);
nyy = ntime/12;
x1 = reshape(x,12,nyy);
t1 = reshape(t,12,nyy);
x1 = x1';
t1 = t1';
% trend = zeros(12,1);
trend = zeros(12,2);
trend1 = zeros(12,2);
pvalue = zeros(12,1);
cons = zeros(12,1);
x0 = ones(nyy,1);
for i = 1:12
%    x1(i,:) = x1(i,:) - mean(x1(i,:)); 
%    tmp_data = polyfit(t1(i,:),x1(i,:),1);
   cons(i,1) = mean(x1(:,i));
   x1(:,i) = x1(:,i) -cons(i,1);
   X = [x0,squeeze(t1(:,i))];
   [b,~,~,~,stats]= regress(x1(:,i),X);
   trend(i,1) = b(2,1);
   trend(i,2) = b(1,1);
   x1(:,i) = x1(:,i) - t1(:,i).*trend(i,1)-trend(i,2);
   pvalue(i,1) = stats(1,3);
end
x1 = x1';
x = reshape(x1,ntime,1);
trend1(1:2,:) = trend(11:12,:);
trend1(3:12,:) = trend(1:10,:);
