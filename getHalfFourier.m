% getHalfFourier.m
%
%        $Id:$ 
%      usage: d = getHalfFourier(im)
%         by: justin gardner
%       date: 07/08/11
%    purpose: returns the half-fourier representation of the image.
%             (See reconstructFromHalfFourier)
%
function d = getHalfFourier(im)

% check arguments
if ~any(nargin == 1)
  help getHalfFourier
  return
end

% make sure there are an odd number of pixels
if iseven(size(im,1)),im = im(1:end-1,:);end
if iseven(size(im,2)),im = im(:,1:end-1);end

% take fourier transform of image
imf = fft2(im);

% get input dimensions
d.originalDims = size(im);

% get one half of fourier image
imfh = fftshift(imf);
imfh = imfh(1:d.originalDims(1),1:ceil(d.originalDims(2)/2));

% extract dc form half fourier image
d.dc = imfh(ceil(d.originalDims(1)/2),end);
halfFourier = imfh(1:(numel(imfh))-ceil(d.originalDims(1)/2));

d.mag = abs(halfFourier);
d.phase = angle(halfFourier);
d.n = length(d.phase);

