function [h] = subplot_yih(nrow,ncol,irow,icol)
%
% function [h] = subplot_yih(nrow,ncol,irow,icol)
%

widthplot = 1 ./ ncol .* .8;
widthmargin = (1 - widthplot .* ncol) ./ (ncol+1.05);
% widthmargin = 0;
% heightplot = 1 ./ (nrow+1) .* 1.;
% heightplot = 1 ./ nrow .* .8;
heightplot = 1 ./ (nrow+1) .* .9;
heightmargin = (1 - heightplot .* nrow) ./ (nrow+1);
% heightmargin = 1-heightplot.*nrow;
% heightmargin = 0;

% h = subplot('position',[widthmargin.*icol+widthplot.*(icol-1) 1-heightplot.*irow widthplot heightplot]);
h = subplot('position',[widthmargin.*icol+widthplot.*(icol-1) 1-heightmargin.*irow-heightplot.*irow widthplot heightplot]);

