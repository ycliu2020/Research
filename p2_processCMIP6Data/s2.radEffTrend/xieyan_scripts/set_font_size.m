function [] = set_font_size(fontsize,haxis)
%
% This program set font size of the current axis
%
% function [] = set_font_size(fontsize,haxis)
% ARGIN
%   fontsize : [16]
%   haxis    : [gca]
% 

if ~exist('fontsize','var') || isempty(fontsize)
  fontsize = 16;
end

if ~exist('haxis','var') || isempty(haxis)
  haxis = gca;
end

set(haxis,'fontsize',fontsize);

htitle = get(haxis,'title');
set(htitle,'fontsize',fontsize);

hxlabel = get(haxis,'xlabel');
set(hxlabel,'fontsize',fontsize);

hylabel = get(haxis,'ylabel');
set(hylabel,'fontsize',fontsize);

