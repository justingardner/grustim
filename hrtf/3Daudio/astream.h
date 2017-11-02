
// astream.h -- header for audio stream implementation (IRIX 4.0.5F)
// language: C++
// author  : Michael A. Casey
// copyright: (c) 1994 MIT Media Lab, All Rights Reserved

#ifndef _astream_h
#define _astream_h

#define SHORTSCALE (32767.0)
#define QUEUE_FRAMES (10)
#define FRAME_OFFSET (0)
#define MAX_VAL (1.0)
#define GAIN (64)
#define ATTEN (32)



#include "audio.h"
#include "fstream.h"
#include "iostream.h"

typedef float* sampPtr;
typedef float  sampT;
typedef long len_t;

// declare audio stream class, derived from ALport
class astream{
protected:
  ALconfig inconfig,outconfig;
  ALport   inport,outport;
  len_t srate;
  len_t length;
  len_t channels;
  short* sbuf;

  void initAudio(len_t);     // initialize the audio hardware
  void newconfig(void);     // make audio configuration struct
  len_t rateToSrate(len_t);   // convert SRATE to srate

public:
  const char* fname;
  ofstream* outFile;

  astream(len_t,len_t,len_t);     // default contructor
  astream(len_t,len_t,len_t,const char*); // sample rate, output file name constructor
  ~astream(); // destructor

  void operator<<(sampPtr);  // stream samples to output
  void operator>>(sampPtr);  // stream samples to input
};


#endif

