function [s] = hrtfpath(root,dir_ch,subdir,select,ext,elev,azim)
%
% function [s] = hrtfpath(root,dir_ch,subdir,select,ext,elev,azim)
% Return pathanme for HRTF data file:
%	root is root directory.
%	dir_ch is directory character, '/' (unix) or ':' (mac).
%	subdir is 'compact', 'full', etc.
%	select is 'L', 'R' or 'H'.
%	ext is the filename extension '.dat', etc.
%	elev is elevation.
%	azim is azimuth.
%
s = sprintf('%s%s%s%selev%d%s%s%de%03da%s',...
	root,dir_ch,subdir,dir_ch,round(elev),...
	dir_ch,select,round(elev),round(azim),ext);
