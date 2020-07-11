function [mean_v] = area_mean(v,latitude)
% function [mean_v] = area_mean(v,latitude)

if ~exist('latitude','var')
  latitude = ones(size(v));
elseif (ndims(v) ~= ndims(latitude)) || any(size(v) ~= size(latitude))
  disp('area_mean: Mismatching dimensions!');
  mean_v = [];
  return;
end

idx = ~isnan(v);
wgt = cos(latitude./180.*pi);
if sum(idx(:)) > 0
  mean_v = sum(v(idx).*wgt(idx)) ./ sum(wgt(idx));
else
  mean_v = NaN;
end

