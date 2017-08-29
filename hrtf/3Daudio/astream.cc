
// astream.cc -- definitions for astream class members
// language: C++
// author  : Michael A. Casey
// copyright: (c) 1994 MIT Media Lab, All Rights Reserved

#include <iostream.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "astream.h"

extern float fabsf(float);

#define DEBUG 1
#define DCOUNT (500)

// local prototypes
void error(const char*, const char*);
int _audio_inited=0;

// ::astream -- default constructor
astream::astream(len_t sr, len_t length, len_t channels){

  if(!_audio_inited)
    initAudio(sr);
  else
    astream::srate=rateToSrate(sr);
  astream::length=length; // set the buffer length
  astream::channels=channels; // set the number of channels
  newconfig();
#ifdef DEBUG
  cout << "\nbuffer length = " << astream::length;
  cout << "\nqueue frames  = " << QUEUE_FRAMES;
  cout << "\nnum channels  = " << astream::channels << "\n";
  cout.flush();
#endif
}

// ::astream -- sample rate parameter constructor
astream::astream(len_t sr, len_t length, len_t channels,
		 const char fname[]="snd.out"){

  if(!_audio_inited)
    initAudio(sr);
  else
    srate=rateToSrate(sr);  
  astream::length=length;
  astream::channels=channels;
  
  // do all the audio library stuff
  newconfig();
  
  // do the stream library stuff (for output file only, file of ASCII chars)
  outFile = new ofstream(fname,ios::out);
  if(!outFile)
    error("Couldn't open sound output file","<astream::astream>");
  astream::fname = fname;
#ifdef DEBUG
  cout << "\nbuffer length = " << astream::length;
  cout << "\nqueue frames  = " << QUEUE_FRAMES;
  cout << "\nnum channels  = " << astream::channels << "\n";
  cout.flush();
#endif  
}

// ::newconfig -- build the audio configuration struct for this stream
void astream::newconfig(){
  inconfig = ALnewconfig();
  ALsetqueuesize(inconfig, length*QUEUE_FRAMES);  /* Set Queue Size */

  if(channels==1)
    ALsetchannels(inconfig, AL_MONO);            /* 1channel */
  else
    ALsetchannels(inconfig, AL_STEREO);          /* 2 channels */
  
  ALsetsampfmt(inconfig, AL_SAMPFMT_TWOSCOMP);   /* two's compliment */
  ALsetwidth(inconfig, AL_SAMPLE_16);            /* 16 bit ints (shorts) */

  outconfig = ALnewconfig();
  ALsetqueuesize(outconfig, length*QUEUE_FRAMES);  /* Set Queue Size */

  if(channels==1)
    ALsetchannels(outconfig, AL_MONO);            /* 1 channels */
  else
    ALsetchannels(outconfig, AL_STEREO);          /* 2 channels */
  
  ALsetsampfmt(outconfig, AL_SAMPFMT_TWOSCOMP);   /* two's compliment */
  ALsetwidth(outconfig, AL_SAMPLE_16);            /* 16 bit ints (shorts) */
  sbuf = new short[length*QUEUE_FRAMES];   /* Allocate shorts space */
  
  if(!sbuf) error("Can't allocate for sbuf","<::astream()>");
  
  /* Open the audio ports */
  inport = ALopenport("input", "r", inconfig); /* read port */
  if (!inport) error("couldn't open audio input port","");
  outport = ALopenport("output", "w", outconfig); /* write port */
  if (!outport) error("couldn't open audio output port","");

}  

// ::~astream -- audio class destructor
astream::~astream(){
  ALfreeconfig(inconfig);  // free config structures
  ALfreeconfig(outconfig);  
  ALcloseport(inport);     // free audio ports
  ALcloseport(outport);
}

void astream::initAudio(len_t sr){
  len_t pvbuf[2];

  
  srate=rateToSrate(sr);       //  get the sample rate identifier

  /* Configure Audio Hardware */
  pvbuf[0] = AL_MONITOR_CTL; pvbuf[1] = AL_MONITOR_OFF; /* Input Monitor */
  ALsetparams(AL_DEFAULT_DEVICE, pvbuf, 2);  

  pvbuf[0] = AL_LEFT_INPUT_ATTEN; pvbuf[1] = ATTEN;   /* Input Attenuation L */
  ALsetparams(AL_DEFAULT_DEVICE, pvbuf, 2);

  pvbuf[0] = AL_RIGHT_INPUT_ATTEN; pvbuf[1] = ATTEN;  /* Input Attenuation R */
  ALsetparams(AL_DEFAULT_DEVICE, pvbuf, 2);

  pvbuf[0] = AL_INPUT_RATE; pvbuf[1] = srate;         /* Input Sample Rate */
  ALsetparams(AL_DEFAULT_DEVICE, pvbuf, 2);

  pvbuf[0] = AL_INPUT_SOURCE; pvbuf[1] = AL_INPUT_LINE;/* Input Device */
  ALsetparams(AL_DEFAULT_DEVICE, pvbuf, 2);

  pvbuf[0] = AL_OUTPUT_RATE;pvbuf[1]=srate;          /* Output Sample Rate */
  ALsetparams(AL_DEFAULT_DEVICE, pvbuf, 2);

  pvbuf[0] = AL_LEFT_SPEAKER_GAIN; pvbuf[1] = GAIN;  /* Output Gain L */
  ALsetparams(AL_DEFAULT_DEVICE, pvbuf, 2);

  pvbuf[0] = AL_RIGHT_SPEAKER_GAIN; pvbuf[1] = GAIN; /* Output Gain R */
  ALsetparams(AL_DEFAULT_DEVICE, pvbuf, 2);

}



// ::rateToSrate -- convert sample rate into AL's srate identifier
long astream::rateToSrate(long rate)
{
  long srate;

  switch(rate){
  case 8000:
    srate = AL_RATE_8000;
    break;
  case 11025:
    srate = AL_RATE_11025;
    break;
  case 16000:
    srate = AL_RATE_16000;
    break;
  case 22050:
    srate = AL_RATE_22050;
    break;
  case 32000:
    srate = AL_RATE_32000;
    break;
  case 44100:
    srate = AL_RATE_44100;
    break;
  case 48000:
    srate = AL_RATE_48000;
    break;
  default:
    error("unrecognised sample rate","<rateToSrate>");
  }
  return srate;
}


// ::operator<< -- send the samples in the sample struct to the audio port
//
//
//
void astream::operator<<(sampPtr s)
{
  // Write the samples to the output buffer
  for(int i=0;i<length;i++)
    sbuf[i]=(short) s[i];
  ALwritesamps(outport, sbuf, length);
}


// ::operator>> -- read ::length samples into the sample buffer
// It is up to the user to make sure that there are "length"
// slots available for sample input in the array.
//
void astream::operator>>(sampPtr s)
{
len_t i;
len_t l;
	if ((l = ALgetfilled(inport)) > length * (QUEUE_FRAMES - 1)) {
		printf(".");
		fflush(stdout);
	}
  
  // read the samples from the input port
    ALreadsamps(inport, sbuf, length);
  
  // Transfer samples to the sampT buffer re-casting as necessary
  // STEREO input is interleaved L R L R alternately
  for (i=0;i<length;i++)
    s[i]=(sampT) sbuf[i];
} 

/* This code already exists in mstream.cc
// report errors
void error(const char*a, const char*b=""){
  cout << a << "," << b << "\n";
  exit(1);
}
*/
