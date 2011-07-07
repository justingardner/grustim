function makeSound(filename,freq,len,amp)
% makeSound.m
% jan 6 2010
%
%makeSound(filename,freq,len,amp)
%filename should be string such as 'sine400' and extension will be added later by
%the function.  frequency is in hertz, default is 400. length is in
%seconds, default is 1.  amp will adjust volume, defaults to 1.
%
%this script will generate a sine400.wav file in your working directory.
%you can choose to convert the .wav file to .aif format (for apple) using
%itunes.  open itunes, drag the .wav file into library, right click and
%select "create aiff file."  look for the aiff file in your itunes folder.
%
%for use in mgl stimuli, use functions mglInstallSound and mglPlaySound.


if ~any(nargin == [1 2 3 4])
    help makeSound;
    return
end

% some defaults
if ieNotDefined('freq'),freq=400;end
if ieNotDefined('len'),len=1;end
if ieNotDefined('amp'),amp=1;end

% sample rate    
sampleRate = 20000;
x = 0:1/(sampleRate-1):len;
y = amp * sin(2*pi*x*freq);

sound(y,sampleRate);

% ask if we should save
filename = setext(filename,'wav');
if askuser(sprintf('Do you want to save sound %s',filename))
  wavwrite(y,sampleRate,32,filename);   
end