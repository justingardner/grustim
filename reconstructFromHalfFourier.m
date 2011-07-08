% reconstructFromHalfFourier.m
%
%        $Id:$ 
%      usage: im reconstructFromHalfFourier(d)
%         by: justin gardner
%       date: 07/08/11
%    purpose: Reconstructs an image from its half fourier representation
%             (see getHalfFourier)
%
function im = reconstructFromHalfFourier(d)

% check arguments
if ~any(nargin == [1])
  help reconstructFromHalfFourier
  return
end

d.halfFourier = d.mag.*(cos(d.phase)+i*sin(d.phase));

% first make the last column of the half fourier space which includes
% the dc and should have the frequency components replicated corectly
halfFourier = [d.halfFourier d.dc];
halfFourier(end+1:end+floor(d.originalDims(1)/2)) = conj(d.halfFourier(end:-1:end-floor(d.originalDims(1)/2)+1));
halfFourier = reshape(halfFourier,d.originalDims(1),ceil(d.originalDims(2)/2));

% replicate the frequency components to make the negative frequencies which
% are the complex conjugate of the positive frequncies
halfFourier2 = fliplr(flipud(conj(halfFourier(:,1:floor(d.originalDims(2)/2)))));
im = ifft2(ifftshift([halfFourier halfFourier2]));

