
// mstream.cc -- definitions for MIDI stream implementation (IRIX 4.0.5F)
// language: C++
// author  : Michael A. Casey
// copyright: (c) 1994 MIT Media Lab, All Rights Reserved


#include "iostream.h"
#include "mstream.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// Prototypes
extern float powf(float,float);
void error(const char*, const char*);

// Macros
#define MIDI_NOTE_ZERO (8.1757989)   // MIDI Note 0 in Hz

// global initialization flag (for all objects of this class)
int _midi_inited=0;

// ::mstream -- constructor for mstream class
mstream::mstream(){
  u_int pbuf[3];

  if(_midi_inited) // error scope should be resolved by argument types
    ::error("You should only open one MIDI stream object","<mstream::mstream>");
  else
    _midi_inited=1;
  
  // Set up the configuration for real-time MIDI control
  config = new MIconfig;
  pbuf[0]=MI_TIMEOUT;pbuf[1]=0;pbuf[2]=0;    /* Don't wait for input */
  config->setparams(pbuf,3);
  pbuf[0]=MI_BLOCKING;pbuf[1]=MINONBLOCKING; /* Don't block if no input */
  config->setparams(pbuf,2);
  pbuf[0]=MI_STAMPING;pbuf[1]=MINOSTAMP;     /* Don't use time stamping */
  config->setparams(pbuf,2);
  

  // Open the MIDI port
  // 
  // If "startmidi" has not been called opening the port will fail.
  // Call start MIDI before running this program.
  //  
  port = new MIport;
  if (port->open("rw", config) < 0) 
    ::error("couldn't open MIDI input port","\nHave you executed startmidi?");

  /* Set the MIDI error handler */
  /*MIseterrorhandler((MIerrfunc) mstream::midiError());*/
}

// Destructor 
mstream::~mstream(){
  // insert destructor code here
}

// operator>> -- get MIDI message from input stream, parse and make into
//               a midiStruct message
midiStruct& mstream::operator>>(midiStruct& M){
  MIevent e;    // declare a MIDI event structure
  MImessage m;  // declare a MIDI message structure
  int retVal;

  M.ctrl=0; // reset continuous controllers
  if((retVal=port->receive(&e,1))<0)
    ::error("MIDI receive error:","[operator>>]");
  if(!retVal){ // if there is no new MIDI message
    M.valid=0;
    return M;
  }
  
  if(e.t!=MIMESSAGE){
    printf("Wrong Message Type: SYSEX\n");
	M.valid=0;
    return M;
  }
  
  m=MESSAGE(e); // Get the message part of the MIDI event
  M.valid=1;
  M.channel=m.channel();  
  switch(m.status()){  
  case MIDI_NoteOn: case MIDI_NoteOff:
    if(m.byte2() && (m.status()==MIDI_NoteOn)){
      M.newNote=1;
      M.freq=noteToHertz(m.byte1());
      M.amp =velToVol(m.byte2());
      M.ctrl=0;
    }
    else // It\'s a note off
      if(M.freq==noteToHertz(m.byte1())){
	M.endNote=1;
	M.ctrl=0;
      }
    break;

  case MIDI_ControlChange:
    M.ctrl = m.byte1();
    M.val  = m.byte2();
    break;

  default:
    ;
  }

  return M;
}


// cnvrtNoteToHz 
// turn a MIDI note-on event into Hz
//
float mstream::noteToHertz(u_char MIDInote)
{
  float freqval;
  freqval=MIDI_NOTE_ZERO*powf(2.0,(float)MIDInote/12.0);
  return freqval;
}


// cnvrtVelToVol
// scale volume messages to lie in range 0..1
//
float mstream::velToVol(u_char vel)
{
  float vol;

  vol=(float)vel/127.0;
  return vol;
}


// midiError -- error handler for MIDI errors
void mstream::midiError(long err, const char *s)
{
  printf("MIDI error%d, %s\n",err,s);
  exit(1);
}



void error(const char* s, const char* t){
  cout << s << ":" << t << "\n";
  exit(-1);
}
  
