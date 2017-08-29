function [x] = readhrtf(elev,azim,select)
%
% function [x] = readhrtf(elev,azim,select)
%
% elev is elevation from -40 to 90 degrees
% azim is azimuth from 0 to 180 degrees
% select is:
%	'L' use full data from left pinna
%	'R' use full data from right pinna
%	'H' use compact data
% Returns stereo symmetrical hrtf in first two rows of
% x such that left is first row, right is second row.
%
% Bill Gardner
% Copyright 1995 MIT Media Lab. All rights reserved.
%

%
% Root directory for Macintosh or UNIX. Put your own in here.
%
if (length(getenv('MACHTYPE')) == 0)
	root = 'Pokey:hrtf';
	dir_ch = ':';
else
	root = '/grustim/hrtf';
	dir_ch = '/';
end
%
% check arguments
%
azim = round(azim);
% if ((azim < 0) | (azim > 180))
% 	error('azimuth must be between 0 and 180 degrees');
% end
if ((elev < -40) | (elev > 90))
	error('elevation must be between -40 and 90 degrees');
end


%
% format filename
%
flip_azim = 360 - azim;
if (flip_azim == 360)
	flip_azim = 0;
end
ext = '.dat';
if (select == 'L')
	pathname = hrtfpath(root,dir_ch,'full',select,ext,elev,azim);
	x(1,:) = readraw(pathname);
	pathname = hrtfpath(root,dir_ch,'full',select,ext,elev,flip_azim);
	x(2,:) = readraw(pathname);
elseif (select == 'R')
	pathname = hrtfpath(root,dir_ch,'full',select,ext,elev,flip_azim);
	x(1,:) = readraw(pathname);
	pathname = hrtfpath(root,dir_ch,'full',select,ext,elev,azim);
	x(2,:) = readraw(pathname);
elseif (select == 'H')
	pathname = hrtfpath(root,dir_ch,'compact',select,ext,elev,azim);
	tmp = readraw(pathname);
	x(1,:) = tmp(1:2:length(tmp));
	x(2,:) = tmp(2:2:length(tmp));
else
	error(sprintf('%s not a valid selection, use L, R, or H',select));
end
