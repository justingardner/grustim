function y = invert(x,lim,n)
%
% y = invert(x[,lim,n])
%
% creates inverse filter for x using fft, 1/X, ifft
%	lim is maximum range of inverse filter in dB (default 20).
%	n is size of inverse filter (default length(x))
%
% Bill Gardner
% Copyright 1995 MIT Media Lab. All rights reserved.
%
%
if (nargin < 3)
	n = length(x);
end
if (nargin < 2)
	lim = 20;
end

X = fft(x,n);
Xlogmag = 20*log10(abs(X));
%
% Find average log magnitude over 6 octaves, ignoring one octave
% each at top and bottom of usual 8 octave range.
%
lo = ceil(n / (2^8));
hi = floor(n / (2^2));
Xavglogmag = sum(Xlogmag(lo:hi)) / (hi-lo+1);

%
% Subtract average, clip, and restore average.
%
halflim = lim / 2;
for k = 1 : n
	Xlogmag(k) = Xlogmag(k) - Xavglogmag;
	if (Xlogmag(k) > halflim)
		Xlogmag(k) = halflim;
	elseif (Xlogmag(k) < -halflim)
		Xlogmag(k) = -halflim;
	end
	%
	% Commenting out this line gives inverse filter that
	% varies around 0 dB.
	%
	Xlogmag(k) = Xlogmag(k) + Xavglogmag;
end

%
% Y is inverse spectrum (zero phase).
%
Y = zeros(size(X));
Y = 10.^(-Xlogmag/20) .* exp(-i * angle(X));

%
% Take ifft and rotate so that maximum value is centered.
%
y = real(ifft(Y));
[ymax,ymaxi] = max(y);
y = rotate_vect(y, floor(n/2) - ymaxi);

