function plotmag(x,fs,n)
% plotmag(x[,fs,n])
% Plot spectrum magnitude
%
% Bill Gardner
% Copyright 1995 MIT Media Lab. All rights reserved.
%

if (nargin < 3)
	n = max(size(x));
end
if (nargin < 2)
	fs = 44100;
end
fx = fft(x,n);
freq = (0 : n/2) * fs / n;
plot(freq,20*log10(abs(fx(1:(n/2 + 1)))));
grid
title('magnitude');
ylabel('dB');
xlabel('freq, Hz');

