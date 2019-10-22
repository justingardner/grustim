% cameraTestCalibrate.m
%
%      usage: cameraTestCalibrate()
%         by: justin gardner
%       date: 10/15/19
%    purpose: tries to calibrate time of acquisition
%
function retval = cameraTestCalibrate(e)

% check arguments
if ~any(nargin == [1])
  help cameraTestCalibrate
  return
end


% get first set of camera images
c = e.stimuli{1}.cameraImages{1};

% bring up figure with some images
mlrSmartfig('cameraTestCalibrate','reuse');clf;

startImage = 50;
skipImage = 50;
nDisp = floor((350-startImage)/skipImage);

for iDisp = 1:nDisp;
  thisImage = startImage+(iDisp-1)*skipImage;
  %subplot(1,nDisp,iDisp);
  imagesc(c.im(:,:,thisImage)');
  colormap(gray);
  imageTimeStamp(iDisp) = c.t(thisImage);
%  title(sprintf('Time: %i',imageTimeStamp(iDisp)));
  title(sprintf('Time: %.4f',imageTimeStamp(iDisp)));
  r = input('What number do you see: ','s');
%  timestamp(iDisp) = str2num(r);
%  disp(sprintf('Frame: %i happened at: %.4f',imageTimeStamp(iDisp),timestamp(iDisp)));
end

imageTimeStampOffset = min(imageTimeStamp);
imageTimeStamp = imageTimeStamp-imageTimeStampOffset;
imageTimeStampStd = std(imageTimeStamp);
imageTimeStamp = imageTimeStamp/imageTimeStampStd;

% get slope and offset
x = imageTimeStamp(:);x(:,2) = 1;
fit = (((x' * x)^-1) * x')*timestamp(:);

clf;subplot(1,2,1);
plot(imageTimeStamp,timestamp,'ko');
hold on
plot(imageTimeStamp,x*fit,'r-');
keyboard


c = e.stimuli{1}.cameraImages{end};
% now convert all time stamps
c.tz = (c.t-imageTimeStampOffset)/imageTimeStampStd;
c.tcalibrated = c.tz*fit(1) + fit(2);


startImage = 1;
skipImage = 30;
nDisp = floor(350/skipImage);
for iDisp = 1:nDisp;
  subplot(1,2,2);cla;
  thisImage = startImage+(iDisp-1)*skipImage;
  imagesc(c.im(:,:,thisImage)');
  colormap(gray);
  imageTimeStamp(iDisp) = c.tcalibrated(thisImage);
  title(sprintf('Time: %f',imageTimeStamp(iDisp)));
  r = input('What number do you see: ','s');
%  timestamp(iDisp) = str2num(r);
%  disp(sprintf('Frame: %i happened at: %.4f',imageTimeStamp(iDisp),timestamp(iDisp)));
end


keyboard