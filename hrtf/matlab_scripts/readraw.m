function [x] = readraw(pathname)
%
% function [x] = readraw(pathname)
% read raw HRTF data, big-endian (Motorola) format
%
% Bill Gardner
% Copyright 1995 MIT Media Lab. All rights reserved.
%
fid = fopen(pathname,'r','ieee-be');
if (fid == -1)
	error(sprintf('cannot open file %s',pathname));
end
x = fread(fid,inf,'short');
fclose(fid);
%
% return as row vector, +/- 1 max.
%
x = x' / 32768;
