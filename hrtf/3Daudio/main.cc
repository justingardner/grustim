
// hrtf main.cc -- 3Daudio top-level program
//                 senses MIDI, inputs and outputs audio
// Author: Michael A. Casey
// Language: C++
// Copyright (c) 1994 MIT Media Lab All Rights Reserved

#include <stdlib.h>

#include "astream.h"       // audio stream classes and types
#include "mstream.h"       // MIDI  stream classes and types

extern "C"{
#include "complex.h"	   // header for complex data type
#include "hrtf.h"          // header for do_init and do_ctrl
#include "procbuf.h"       // header for do_buf
#include "arraylib.h"	   // header for f_alloc, etc.
#include "misc.h"          // misc stuff
};

#define SAMPLERATE	32000
#define CHANNELS 2 // number of audio channels of output
#define BUFLEN 256 // For  stereo make sure BUFLEN is twice actual size
#define COUNT 20   // number of buffers to calculate before sampling MIDI

void main(void){
  float *zbuf;
  int i;

#ifdef PROFILE
  long loop_count=0;			// Profiling counter
#endif
  
  // construct audio stream object
  astream AUDIO(SAMPLERATE,BUFLEN,CHANNELS);
  
  // Construct MIDI objects
  
  mstream MIDI;
  midiStruct M;
  
  
  // allocate output buffer
  sampPtr ibuf; // 3Daudio allocates input buffer
  sampPtr obuf; // 3Daudio allocates output buffer
  
  do_init(&ibuf);		// 3Daudio initialization
  
  // Control Loop
  int doing_the_thang=1;
  
  // Output QUEUE_FRAMES / 2 buffers of zeros
  zbuf = f_alloc(BUFLEN);
  for (i = 0; i < (QUEUE_FRAMES / 2); i++)
    AUDIO << zbuf;

  cout << "\n3Daudio active...\n";
  cout.flush();
  while (doing_the_thang) {
    
    
    // Sense all buffered MIDI control messages
    MIDI >> M; // look at MIDI port
    while(M.valid){ // if it's a control message
      if(M.ctrl)
	do_ctl(M.ctrl,M.val);	// handle control event
      MIDI >> M;		// look at MIDI port and do it all over again
    }
    
    int flag = TRUE;		// set control flag
    int n = COUNT;		// set loop count
    
    // Audio Loop
    do {
      AUDIO >> ibuf;		//read input signal into ibuf;
      obuf = (sampPtr) do_buf(&ibuf, flag); // process signal
      flag = FALSE;		// unset control flag
      AUDIO << obuf;		// write to output
    } while (n-- > 0);		// audio loop
    
#ifdef PROFILE
    if(loop_count++>100)  break;	// Break hook for profiling
#endif
    
  }				// control loop
}				// end of program

