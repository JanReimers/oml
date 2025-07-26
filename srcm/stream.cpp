module;
#include <iostream>
#include <cassert>

export module oml.StreamableObject;

// stream.h
typedef const char* c_str;

/*! \class StreamableObject stream.h oml/stream.h
  \brief Common IO featrues for all OML containers.

  The class stores a static flag indicating the type of output. THere are three options:
  -binary Condensed/fast  and mostly unreadable.
  -ascii  Condensed/slow but readable.
  -pretty Formatted for reading by humans.
 */
export class StreamableObject
{
 public:
  enum Mode {binary, ascii, pretty};

  //! Outputs the IO mode and type name.
  void WriteHeader(std::ostream&,c_str type) const;
  //! Reads in IO mode and class type name and compares it ith the expected type name
  Mode ReadHeader (std::istream&,c_str expected_type)      ;

  static  c_str PeekAtName(std::istream&);
  static  void  CheckName (std::istream&,c_str);

  static  Mode GetOutputMode() {return theMode;}
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

  static  bool CheckArrayTypes(bool);

 private:
  static  Mode theMode;
  static  bool theCheckArrayTypes;
};

// binio.h
export template <class T> inline void BinaryWrite(const T& t,std::ostream& os) {os << t;}
export template <class T> inline void BinaryRead (      T& t,std::istream& is) {is >> t;}

#define INTRINSIC_IO(Type) \
   template <> inline void BinaryWrite( Type const &t,std::ostream& os)  {os.write((const char*)&t,sizeof( Type));} \
   template <> inline void BinaryRead ( Type       &t,std::istream& is)  {is.read ((      char*)&t,sizeof( Type));} \

  INTRINSIC_IO(int)
  INTRINSIC_IO(unsigned int)
  INTRINSIC_IO(long int)
  INTRINSIC_IO(unsigned long int)
  INTRINSIC_IO(float)
  INTRINSIC_IO(double)
  INTRINSIC_IO(long double)
  INTRINSIC_IO(short)
  INTRINSIC_IO(char)
  INTRINSIC_IO(bool)
  INTRINSIC_IO(void*)

inline std::istream& operator>>(std::istream& is,void* t) 
{
	long l;
	is >> l >> std::ws ;
	t=reinterpret_cast<void*>(l);
	return is;
} 


export template <typename T> inline void BinaryWrite(const std::complex<T>& t,std::ostream& os) 
{
	BinaryWrite(t.real(),os);
	BinaryWrite(t.imag(),os);
}

export template <typename T> inline void BinaryRead(std::complex<T>& t,std::istream& is) 
{
	double re,im;
	BinaryRead(re,is);
	BinaryRead(im,is);
	t=std::complex<T>(re,im);
}

std::string ReadName  (std::istream&);

StreamableObject::Mode StreamableObject::theMode = StreamableObject::ascii;
bool StreamableObject::theCheckArrayTypes=true;

StreamableObject::Mode StreamableObject::SetOutputMode(Mode newMode)
{
  Mode ret=theMode;
  theMode=newMode;
  return ret;
}

bool StreamableObject::CheckArrayTypes(bool check)
{
    bool ret=theCheckArrayTypes;
    theCheckArrayTypes=check;
    return ret;
}

void StreamableObject::WriteHeader(std::ostream& os,c_str type) const
{
  assert(os);
  assert(type);
  if (!Pretty())  os << (int)theMode << " " << type << " ";
  assert(os);
}

StreamableObject::Mode StreamableObject::ReadHeader(std::istream& is,c_str type)
{
  assert(is);
  assert(type);
  Mode current,temp;
  int itemp;
  is >> itemp;
  temp=static_cast<Mode>(itemp);
  assert(temp!=pretty);
  current=SetOutputMode(temp);
  CheckName(is,type);
  assert(is);
  return current;
}

//
//  Check tests if the object type on the stream match this's type.
//  CheckName extracts the name from the stream.
//
void StreamableObject::CheckName(std::istream& is,c_str type)
{
  assert(is);
  assert(type);
  std::string Name=ReadName(is);
  if (theCheckArrayTypes && Name!=type)
  {
    std::cerr << "Read '" << Name << "' from stream, but actual object type is '" << type << "'" << std::endl;
    exit(-1);
  }
}

//
//  ReadName extracts the name from the stream.
//

std::string ReadName(std::istream& is)
{
  assert(is);
  is >> std::ws;
  std::string Name;
  is >>  Name; //Apparently this also reads in the space after the name!
  if (StreamableObject::Binary()) is.get();
  return Name;
}

//
//  Peek at name reads the name from the stream without
//  extracting it.
//
c_str StreamableObject::PeekAtName(std::istream& is)
{

  assert(is);
  is >> std::ws;
  static std::string Name;
  std::streampos beforeName=is.tellg();
  Mode temp;
  int itemp;
  is >> itemp;            //Get the mode flag out of the way.
  (void)temp; //Avoid unused warning
  temp=static_cast<Mode>(itemp);
  is >>  Name;
  is.seekg(beforeName); //Reset stream before start of name and mode flag.
  assert(is);
  return Name.c_str();
}

