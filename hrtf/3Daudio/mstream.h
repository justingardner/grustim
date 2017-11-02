
// mstream.h -- header for MIDI stream implementation (IRIX 4.0.5F)
// language: C++
// author  : Michael A. Casey
// copyright: (c) 1994 MIT Media Lab, All Rights Reserved

#ifndef _mstream_h
#define _mstream_h


// implement MIDI Library Routines as C calls, not C++
#include "midi.h"
#include "midiio.h"

// declare structure to hold MIDI information
typedef struct {
float freq;
float amp;
int valid;
int newNote;
int endNote;
int channel;
int notenum;
int velnum;
int ctrl;
int val;
} midiStruct;

class mstream{
protected:
  MIconfig* config;        // MIDI port configuration
  MIport*   port;          // MIDI port structure
  
  void midiError(long,const char*); // MIDI error handler
  float noteToHertz(u_char);
  float velToVol(u_char);
public:
  mstream();     // default contructor
  ~mstream(); // destructor
  midiStruct& operator>>(midiStruct&);  // stream samples to output
};
  

#endif

