
3Daudio Demonstration Program

Bill Gardner (with help from Mike Casey and Keith Martin)
MIT Media Lab
20 Ames St
Cambridge, MA 02139
billg@media.mit.edu

This program demonstrates how to construct a 3D audio display using
the KEMAR HRTFs provided on our web site. The 3Daudio program will
only run on SGI workstations. The source code is provided so that
interested people can figure out how it works and hack up their own
applications.  For more general information, see the following:

ftp://sound.media.mit.edu/pub/Data/KEMAR/KEMAR-FAQ.txt
ftp://sound.media.mit.edu/pub/Papers/billg-TR-342.ps.Z

****** What does it do?

3Daudio is a realtime "spatializer" that convolves monophonic sound
with a pair of HRTFs to produce a binaural output that should be
listened to via headphones. Using a MIDI interface, you can move the
sound around your head by controlling elevation, azimuth, and
attenuation.  The input sound is taken from the left channel line
input of the SGI, and the binaural output is produced on the stereo
line outputs.

****** What equipment is required?

You need an SGI workstation (Indy or Indigo should work), a serial port
MIDI interface (Macintosh MIDI interfaces work), a MIDI controller capable
of sending continuous controls 1, 2, and 3, a source of audio, and
a pair of headphones.

You also need to have downloaded the 32 kHz diffuse-field equalized
HRTF measurements (ftp://sound.media.mit.edu/pub/Data/KEMAR/diffuse32k.tar.Z)
and the 3Daudio program (ftp://sound.media.mit.edu/pub/Data/KEMAR/3Daudio.tar.Z).

****** How should the hardware be connected?

Connect the serial MIDI interface to one of the SGI serial ports and
plug your MIDI controller into the interface using a MIDI cable. Connect
a source of audio (CD player, for example) to the line audio input. You
can also use the microphone input, but you will need to explicitly set
microphone input using the audio panel (apanel). Connect your headphones
to the headphone output of the SGI, or use a headphone amplifier connected
to the line audio output.

****** How is the program run?

First you need to set the "HRTFROOT" environment variable to point to
the root directory of the 32kHz diffuse-field equalized  HRTF data. If your
data is in "/usr/pat/hrtf/diffuse32k", then you should type:

setenv HRTFROOT /usr/pat/hrtf/diffuse32k

Don't put a trailing '/' in the directory path.

Then you need to start the MIDI driver. If you connected the MIDI
interface to serial port 2, then type:

startmidi -d /dev/ttyd2

If you get a cryptic error message, it's probably already running, and you
can ignore the error.

Now run the 3Daudio program by typing "3Daudio" with no arguments.
After 3Daudio starts up, you can use the audio panel to set microphone
input. 3Daudio will start up by loading the HRTFs which takes a few
seconds.  The screen should display something like:

buffer length = 256
queue frames  = 10
num channels  = 2
Loading HRTFs from '/usr/pat/hrtf/diffuse32k'
elevation -40...
elevation -30...
elevation -20...
elevation -10...
elevation 0...
elevation 10...
elevation 20...
elevation 30...
elevation 40...
elevation 50...
elevation 60...
elevation 70...
elevation 80...
elevation 90...
elev -40 azim 180 atten 0

3Daudio active...
......

At this point, if 3Daudio starts to display dots ".....", it means the audio
processing is not keeping up. You must kill active processes or use a faster
machine. 3Daudio always displays a few dots right as it starts, as shown above.

****** How do you control the position of the sound?

The MIDI controls are mapped as follows:

MIDI control 1: elevation from -40 to +90 degrees.
MIDI control 2: azimuth from -180 to +180 degrees.
MIDI control 3: attenuation from 0 to 20 dB.

As you move the controls, 3Daudio will display the current source position
in polar coordinates:

elev -10 azim -85 atten 0
elev -10 azim -95 atten 0
elev -10 azim -90 atten 0
elev -10 azim -80 atten 0
elev 0 azim -70 atten 0
elev 10 azim -70 atten 0
elev 0 azim -70 atten 0
elev 0 azim -60 atten 0
elev 0 azim -45 atten 0
elev 0 azim -40 atten 0
...etc

If nothing happens when you move the MIDI controls, then something is
wrong with your MIDI setup.

****** How can changes be made?

The source code and Makefile are provided. You will need the C and C++
compilers, and also the various SGI libraries. Chances are good the
libraries aren't all installed on your machine, so you'll have to dig
around for them. Also, if you're running an old version of the system,
the libraries are different.

****** Specifications?

This spatializer only handles a single source at 32 kHz sampling rate
with no reverberation. It is therefore primitive compared to
commercial spatializers. See the technical report TR-342 for more
information about the inner workings of the spatializer. The software
was written fairly quickly, and there is a lot that could be done to
improve the performance of the system, but it is not our intention to
improve or support this software in any way.

****** Questions?

I'm not going to field endless questions about how to build 3D audio
systems, nor am I going to update this software to increase performance.
There are many things that could be done to increase performance, but
I leave it to the reader to implement them.
