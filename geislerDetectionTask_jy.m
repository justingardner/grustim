% geislerDetectionTask_jy.m
% 
%      usage: geislerDetectionTask_jy.m
%         by: jiwon yeon
%       date: 
%  copyright: (c) 2022 Jiwon Yeon
%    purpose: Replicating Najemnik&Geisler's 2005 study
%
%             The study has two experiments: Detection and Search tasks.
%             The current code only replicates the detection task.
%

% making one line change
%%%% test whether the stimulus is correctly displayed

% build pink noise background
image = createPinkNoise(100,100);

% set target locations


mglOpen(0,650,400);
mglVisualAngleCoordinates(57,[51 38]);
mglClearScreen;
mglScreenCoordinates;

% draw noise background


mglTextDraw('Hello World!',[0 0]);

% The above is drawn on the back-buffer of the double-buffered display
% To make it show up you flush the display.
% This function will wait till the end of the screen refresh
mglFlush;
mglClose


%%%%%%%%%% helper functions
function im = createPinkNoise(sz1, sz2)

%   make the odd size of the image
if mod(sz1,2) == 0, sz1 = sz1-1; end
if mod(sz2,2) == 0, sz2 = sz2-1; end
    
white = randn(sz1, sz2);
white = 255 * (white-min(white(:)))./(max(white(:))-min(white(:)));

fwhite = fft2(white);
phase = angle(f_white);

% create a pink noise filter
pink_filter = zeros(ceil(size(f_white,1)/2), ceil(size(f_white,2)/2));
pink_vector = zeros(size(pink_filter,1),1);
for k = 1:length(pink_vector), pink_vector(k) = 1/k; end
pink_filter(1:end-1, end) = flipud(reshape(pink_vector(1:end-1),[],1));
pink_filter(end, 1:end-1) = fliplr(reshape(pink_vector(1:end-10,1,[]));

% fill a quadrant; dc component stays zero
for row = 1:size(pink_filter,1)-1
    for col = 1:size(pink_filter,2)-1
