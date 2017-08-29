%
% Script to equalize the HRTFs according to the diffuse-field
% average (the inverse of the diffuse-field average is used as
% the equalizing filter).
%
% Bill Gardner
% Copyright 1995 MIT Media Lab. All rights reserved.
%

root = '/ti/u/billg/hrtf';
disp(sprintf('root directory is "%s"',root))
dir_ch = '/';

elevs = [-40 -30 -20 -10 0 10 20 30 40 50 60 70 80 90;
	56 60 72 72 72 72 72 60 56 45 36 24 12 1];

%
% Create 'diffuse' directory structure where results will go.
%
disp('creating diffuse directories')
eval(sprintf('!mkdir %s%sdiffuse',root,dir_ch));

%
% Create elev subdirectories.
%
for elev_index = 1 : length(elevs)
	elev = elevs(1,elev_index);
	eval(sprintf('!mkdir %s%sdiffuse%selev%d',root,dir_ch,dir_ch,elev));
end

disp('calculating diffuse-field average')
diffuse_avg

disp('performing diffuse-field equalization')
%
% Normalizing gain - need to run this once with gain = 1, look at
% maximum value, and set gain = 1 / maxval, delete directories, and
% run it again.
%
gain = 1 / 1.08;

num = 0;
maxval = 0;
ycrop = zeros(1,256);

for elev_index = 1 : length(elevs)
	elev = elevs(1,elev_index);
	disp(sprintf('processing elevation %d',elev));
	n_azim = elevs(2,elev_index);
	azim_incr = 360 / n_azim;
	for azim_index = 0 : n_azim - 1
		azim = azim_incr * azim_index;
		%
		% Read symmetrical responses.
		%
		if (azim > 180)
			break;
		end
		flip_azim = 360 - azim;
		if (flip_azim == 360)
			flip_azim = 0;
		end
		pathname = hrtfpath(root,'/','full','L','.dat',elev,azim);
		xl = readraw(pathname);
		pathname = hrtfpath(root,'/','full','L','.dat',elev,flip_azim);
		xr = readraw(pathname);
		%
		% Convolve with equalizing filter.
		%
		yl = conv(xl,df_inv_minph) .* gain;
		yr = conv(xr,df_inv_minph) .* gain;
		%
		% Crop and merge results to stereo file.
		%
		ycrop(1:2:256) = yl(26:(26+127));
		ycrop(2:2:256) = yr(26:(26+127));
		pathname = hrtfpath(root,'/','diffuse','H','.dat',elev,azim);
		writeraw(pathname,ycrop);
		%
		% Keep track of maximum value and file.
		%
		[tmaxval,tmaxi] = max(abs(ycrop));
		if (tmaxval > maxval)
			maxval = tmaxval;
			maxi = tmaxi;
			maxname = pathname;
		end
		num = num + 1;
	end
end

%
% If maxval > 1, then samples will have wrapped around in short files.
%
disp(sprintf('created %d files',num));
disp(sprintf('maximum value %f at sample index %d',maxval,maxi));
disp(sprintf('in file "%s"',maxname));

