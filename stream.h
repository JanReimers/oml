// File: StreamableObject.H  mixin class for streamable objects.
#ifndef _StreamableObject_H_
#define _StreamableObject_H_

// Copyright (1994-2003), Jan N. Reimers

#include <iosfwd>


typedef const char* c_str;

/*! \class StreamableObject stream.h oml/stream.h
  \brief Common IO featrues for all OML containers.

  The class stores a static flag indicating the type of output. THere are three options:
  -binary Condensed/fast  and mostly unreadable.
  -ascii  Condensed/slow but readable.
  -pretty Formatted for reading by humans.
 */
class StreamableObject
{
 public:
  enum Mode {binary, ascii, pretty};

  void WriteHeader(std::ostream&,c_str type) const;
  Mode ReadHeader (std::istream&,c_str type)      ;

  static  c_str PeekAtName(std::istream&);
  static  void  CheckName (std::istream&,c_str);

  //! Set the output mode.  Returns current mode incase you need to restore it later.
  static  Mode SetOutputMode(Mode);
  //! Is mode set to binary?
  static  bool      Binary() {return theMode==binary;}
  //! Is mode set to ascii?
  static  bool      Ascii () {return theMode==ascii ;}
  //! Is mode set to pretty?
  static  bool      Pretty() {return theMode==pretty;}
  //! Set output mode to binary. Returns current mode incase you need to restore it later.
  static  Mode SetToBinary() {return SetOutputMode(binary);}
  //! Set output mode to ascii. Returns current mode incase you need to restore it later.
  static  Mode SetToAscii () {return SetOutputMode(ascii );}
  //! Set output mode to ascii. Returns current mode incase you need to restore it later.
  static  Mode SetToPretty() {return SetOutputMode(pretty);}
  
 private:
  static  Mode theMode;
};

void OMLArrayIndexError(int i, int n);
void OMLListIndexError(int i, int n);

#endif //_StreamabelObject_H_
