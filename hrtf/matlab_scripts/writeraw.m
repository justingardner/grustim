function writeraw(pathname,x)
%
% function writeraw(pathname,x)
% write raw data  (+/- 1 max) to big-endian (Motorola) short file.
%
% Bill Gardner
% Copyright 1995 MIT Media Lab. All rights reserved.
%

fid = fopen(pathname,'w','ieee-be');
if (fid == -1)
	error(sprintf('cannot create file %s',pathname));
end
fwrite(fid, x .* 32768, 'short');
fclose(fid);
