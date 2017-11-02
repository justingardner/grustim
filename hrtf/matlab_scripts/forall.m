%
% Generic script to iterate through all HRTFs
%
%
% Bill Gardner
% Copyright 1995 MIT Media Lab. All rights reserved.
%

elevs = [-40 -30 -20 -10 0 10 20 30 40 50 60 70 80 90;
	56 60 72 72 72 72 72 60 56 45 36 24 12 1];

for elev_index = 1 : length(elevs)
	elev = elevs(1,elev_index);
	n_azim = elevs(2,elev_index);
	azim_incr = 360 / n_azim;
	for azim_index = 0 : n_azim - 1
		azim = azim_incr * azim_index;
		disp(sprintf('elev %f azim %f', elev, azim));
		%
		% processing of (elev,azim) goes here
		%
	end
end

