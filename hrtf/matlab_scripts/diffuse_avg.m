%
% Script to create the diffuse-field average of the KEMAR HRTFs (y),
% and compute the minimum-phase inverse filter (yminp).
%
% Bill Gardner
% Copyright 1995 MIT Media Lab. All rights reserved.
%

root = '/ti/u/billg/hrtf';

% First, average the magnitude squared spectra of all HRTFs.  In a
% diffuse field of noise excitation, the sound incident from different
% directions will be uncorrelated, and thus averaging the power spectra
% (magnitude squared) over all directions gives you the average power
% spectrum.

elevs = [-40 -30 -20 -10 0 10 20 30 40 50 60 70 80 90;
	56 60 72 72 72 72 72 60 56 45 36 24 12 1];

Y = zeros(1,512);
num = 0;

for elev_index = 1 : length(elevs)
	elev = elevs(1,elev_index);
	disp(sprintf('processing elevation %d',elev));
	n_azim = elevs(2,elev_index);
	azim_incr = 360 / n_azim;
	for azim_index = 0 : n_azim - 1
		azim = azim_incr * azim_index;
		pathname = hrtfpath(root,'/','full','L','.dat',elev,azim);
		x = readraw(pathname);
		X = fft(x')';
		Y = Y + abs(X).^2;
		num = num + 1;
	end
end

disp(sprintf('%d measurements averaged',num));
Y = sqrt(Y / num);
%
% df_avg is zero-phase diffuse-field average. Rotate to center energy
% in impulse response.
%
df_avg = rotate_vect(real(ifft(Y)),256);
figure
plotmaglogf(df_avg);
title('diffuse-field average');

%
% df_inv is the inverse filter, df_inv_minph is the minimum-phase
% inverse filter. Allow 24 dB of dynamic range in inverse filter.
%
df_inv = invert(df_avg,24);
[df_inv_ceps,df_inv_minph] = rceps(df_inv);
figure
plotmaglogf(df_inv_minph);
title('inverse equalizing filter');
