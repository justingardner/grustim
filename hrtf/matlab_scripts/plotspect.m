function plotspect(x,fs,n)
% plotspect(x[,fs,n])
% Plot spectrum magnitude and phase
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
subplot(2,1,1);
freq = (0 : n/2) * fs / n;
plot(freq,20*log10(abs(fx(1:(n/2 + 1)))));
grid
title('magnitude');
ylabel('dB');
xlabel('freq, Hz');
subplot(2,1,2);
plot(freq,angle(fx(1 : (n/2 + 1))));
grid
title('phase');
ylabel('radians');
xlabel('freq, Hz');
