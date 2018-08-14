function [s, waveform, td] = generateITD(theta, duration,samplesPerSecond)
% samplesPerSecond = 2999999;
fs = samplesPerSecond-1;
t = 0:1/fs:duration;
% frequency of signal in hz
% hz = 440;
% amplitude = 0.5;
% wav = amplitude * sin(2*pi*hz*t);
wav = 0.125*randn(1,length(t));

% radius of head in meter (assuming spherical)
r = 0.0875;
% speed of sound at room temperature (m/s)
c = 346;
% distance from monitor
d = 0.56;
% stimulus.displayDistance;
% for low frequency sound
% td = (2*sind(theta)) * r / c;
% left - right
td = (sqrt((d*tand(theta)-r).^2 + (d)^2) - sqrt((d*tand(theta)+r).^2 + (d)^2))./c;
td_a = 0:1/fs:abs(td);
clear waveform  s
if td > 0
    waveform(1,:) = [wav, zeros(1,length(td_a))];
    waveform(2,:) = [zeros(1,length(td_a)), wav];
elseif td < 0
    waveform(2,:) = [wav, zeros(1,length(td_a))];
    waveform(1,:) = [zeros(1,length(td_a)), wav];
else
    waveform(1,:) = wav;
    waveform(2,:) = wav;
end

s = mglInstallSound(waveform,samplesPerSecond);

