// File: vectorm.cpp  Make a module for vector
module;
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstdint>
#include <complex>
#include <vector>

export module oml.vector;
export
{

typedef int64_t        index_t;


// veclimits.h
//-----------------------------------------------------------------------------
//
//  Default upper and lower bounds.  This where the Fortan lower bound
//  of 1 gets decided.  Just make LOW=0 if want the C default indexing.
//
const index_t LOW=1;
const index_t HIGH=LOW-1;
//-----------------------------------------------------------------------------
/*! \class VecLimits veclimit.h oml/veclimit.h
  \brief Encapsulate all effects of non-zero lower subscript bounds in this class.

  This is just a data structure so the data is public and there is no
  attempt at any encapsulation stuff here, don't need it.
  \nosubgrouping
*/
class VecLimits
{
 public:
  /*! \name Constructors*/
  //@{
  VecLimits(               ); //!<Limits for a null Vector.
  VecLimits(size_t         ); //!<Construct from size, use default lower bound.
  VecLimits(index_t,index_t); //!<Construct from lower and upper bounds.
  //@}
 ~VecLimits();

  void ReBase(int low);

  static size_t  size(index_t,index_t)      ;
  //! Returns number of elements.
  size_t  size(               ) const;

  index_t Offset         (index_t) const;
  bool    Check          (       ) const;
  bool    CheckIndex     (index_t) const;
  //! Comparison.
  bool operator==(const VecLimits&) const;
  //! Comparison.
  bool operator!=(const VecLimits&) const;

  /*! \name IO
    See StreamableObject for details on output formats. You must include std::vector_io.h to get
    definitions.  op>> and op<< are also defined.
  */
  //@{
  std::ostream& Write(std::ostream&) const; //!<Write to stream.
  std::istream& Read (std::istream&)      ; //!<Read from stream.
  //@}

  index_t Low;   //!<lower std::vector index limit.
  index_t High;  //!<upper std::vector index limit.
};

#ifdef DEBUG
#define CHECK Check()
#else
#define CHECK
#endif

//-----------------------------------------------------------------------------
//
//  Constructors.
//
inline VecLimits::VecLimits() :                 //Limits for a null std::vector.
  Low(LOW),
  High(HIGH)
  {
    CHECK;
  }

inline VecLimits::VecLimits(size_t size) :     //Construct from size, use default lower bound.
  Low(LOW),
  High(LOW+size-1)
  {
    CHECK;
  }

inline VecLimits::VecLimits(index_t low, index_t high) : //COnstruct from lower and upper bounds.
  Low(low),
  High(high)
  {
    CHECK;
  }

inline VecLimits::~VecLimits() {}

inline void VecLimits::ReBase(int low)
{
    High=High-Low+low;
    Low=low;
}

//-----------------------------------------------------------------------------
//
//  Static size calculator
//
inline size_t  VecLimits::size(index_t low, index_t high)
{
  return (index_t)(high-low+1);
}

inline size_t  VecLimits::size() const
{
  return size(Low,High);
}

//-----------------------------------------------------------------------------
//
//  Allow == operators for comparison (what else?).
//
inline bool VecLimits::operator==(const VecLimits& lim) const
{
  return (Low==lim.Low)&&(High==lim.High);
}

inline bool VecLimits::operator!=(const VecLimits& l) const
{
  return !((*this)==l);
}

inline index_t VecLimits::Offset(index_t i) const
{
  return i-Low;
}

inline std::ostream& operator<<(std::ostream& os,const VecLimits& v)
{
  return v.Write(os);
}

inline std::istream& operator>>(std::istream& is,VecLimits& v)
{
  return v.Read (is);
}

// Use this to get limits of a tensor product
VecLimits operator*(const VecLimits& a, const VecLimits& b);


//-----------------------------------------------------------------------------
//
//  Index checking with descriptive erros.
//
bool VecLimits::Check() const
{
  return High+1 >= Low;
}

bool VecLimits::CheckIndex(index_t i) const
{
  bool ok=true;
  if (i< Low )
  {
    std::cout << "Index " << i << " too low in VecLimits, Low=" << Low << std::endl;
    ok=false;
  }
  if (i>High)
  {
    std::cout << "Index " << i << " too high in VecLimits, High=" << High << std::endl;
    ok=false;
  }
#ifdef DEBUG
  assert(ok);
#endif
  return ok;
}

VecLimits operator*(const VecLimits& a, const VecLimits& b)
{
    assert(a.Low==b.Low);
    int size=a.size()*b.size();
    return VecLimits(a.Low,a.Low+size-1);
}


#undef CHECK

// shape.h
//--------------------------------------------------
//
//  Containers will have the following traits:
//                       Shape        Store   Data
//  ------------------------------------------------
//  Array               ArrayShape    Full    Real.
//  Vector              VectorShape   Full    Real.
//  Matrix              MatrixShape   Full    Real.
//  SMatrix             MatrixShape   Sym     Real.
//  HMatrix             MatrixShape   Sym     Real.
//  GetColumn(Matrix)   VectorShape   Full    Abstract.
//  Transpose(Matrix)   MatrixShape   Full    Abstract.
//  OuterProduct(V1,V2) MatrixShape   Sym     Abstract.
//

enum Shape   {ArrayShape,VectorShape,MatrixShape};
enum Store   {Full,Symmetric,Diagonal,Sparse,Upper,Lower,Row,Column};
enum Data    {Real,Abstract};


// Define how different storage layouts and Data type interact
template <Store S1, Store S2> class ReturnStore{};
template <Store S> class ReturnStore<S   ,S   >{public: static const Store RetType=S;};
template <Store S> class ReturnStore<Full,S   >{public: static const Store RetType=Full;};
template <Store S> class ReturnStore<S   ,Full>{public: static const Store RetType=Full;};
template <       > class ReturnStore<Full,Full>{public: static const Store RetType=Full;};
template <>        class ReturnStore<Symmetric,Diagonal >{public: static const Store RetType=Symmetric;};
template <>        class ReturnStore<Diagonal ,Symmetric>{public: static const Store RetType=Symmetric;};

template <Data D1, Data D2> class ReturnData;
template <Data D> class ReturnData<D       ,D       > {public: static const Data RetType=D;};
template <Data D> class ReturnData<D       ,Abstract> {public: static const Data RetType=Abstract;};
template <Data D> class ReturnData<Abstract,D       > {public: static const Data RetType=Abstract;};
template <>       class ReturnData<Abstract,Abstract> {public: static const Data RetType=Abstract;};

// stream.h
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
template <class T> inline void BinaryWrite(const T& t,std::ostream& os) {os << t;}
template <class T> inline void BinaryRead (      T& t,std::istream& is) {is >> t;}

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


template <typename T> inline void BinaryWrite(const std::complex<T>& t,std::ostream& os) 
{
	BinaryWrite(t.real(),os);
	BinaryWrite(t.imag(),os);
}

template <typename T> inline void BinaryRead(std::complex<T>& t,std::istream& is) 
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


void OMLArrayIndexError(index_t i, index_t n)
{
  std::cerr << "Array index " << i
	    << " exceeds the array bounds (0:" << n
	    << ")" << std::endl;
}

void OMLListIndexError(index_t i, index_t n)
{
  std::cerr << "List index " << i
	    << " exceeds the list bounds (0:" << n
	    << ")" << std::endl;
}

//-----------------------------------------------------------------------------
//
//  IO
//
std::ostream& VecLimits::Write(std::ostream& os) const
{
  assert(os);
  if(StreamableObject::Binary()) {BinaryWrite(Low,os);BinaryWrite(High,os);}
  if(StreamableObject::Ascii ()) os << Low << " " << High << " ";
  if(StreamableObject::Pretty()) os << "(" << Low << ":" << High << ")";
  assert(os);
  return os;
}

std::istream& VecLimits::Read(std::istream& is)
{
  assert(is);
  if(StreamableObject::Binary()) {BinaryRead(Low,is);BinaryRead(High,is);}
  if(StreamableObject::Ascii ()) is >> Low >> High;
  assert(is);
  return is;
}


template <class T> class OpLT
{
 public:
   typedef bool RetType;
   static inline bool apply(const T a, const T b) {return a<b;}
};

template <class T> class OpLE
{
 public:
   typedef bool RetType;
   static inline bool apply(const T a, const T b) {return a<=b;}
};

template <class T> class OpGT
{
 public:
   typedef bool RetType;
   static inline bool apply(const T a, const T b) {return a>b;}
};

template <class T> class OpGE
{
 public:
   typedef bool RetType;
   static inline bool apply(const T a, const T b) {return a>=b;}
};

//-----------------------------------------------------------
//
//  Specify how different types are allowed to mix in binary operators.
//  Main thing is to define a return type that does not lose any information.
//
template <class T1, class T2> struct ReturnType;
//
//  Allow case where both types are the same.  Assumes no units.
//
template <class T> struct ReturnType<T,T> {typedef T RetType;};
//
//  complex/scalar mixing is allowed.
//
template <class T> struct ReturnType<             T ,std::complex<T>> {typedef std::complex<T> RetType;};
template <class T> struct ReturnType<std::complex<T>,             T > {typedef std::complex<T> RetType;};
//
//  double int is allowed
//
template <> struct ReturnType<double,int   > {typedef double RetType;};
template <> struct ReturnType<int   ,double> {typedef double RetType;};
//---------------------------------------------------------------
//
//  Temporary unary operation holder.
//
template <class T, class TR, class A, Shape S> class XprUnary
{};


template <class T, class TR, class A> class XprUnary<T,TR,A,VectorShape>
{
 public:
  XprUnary(const A& a,TR(*f)(const T&)) : itsA(a), itsF(f) {};
  ~XprUnary() {};

  TR        operator()(index_t n) const {return itsF(itsA(n));}
  size_t    size      (         ) const {return itsA.size();}
  VecLimits GetLimits (         ) const {return itsA.GetLimits();}
 private:
   A itsA;
   TR(*itsF)(const T&) ;
};

//----------------------------------------------------
//
//  Hold a reference to terminals
//
template <class T, class R,Shape> class Ref {};

template <class T, class R> class Ref<T,R,VectorShape>
{
 public:
  Ref(const R& r) : itsRef(r) {};
  T         operator()(index_t n) const {return itsRef(n);}
  size_t    size      (         ) const {return itsRef.size();}
  VecLimits GetLimits (         ) const {return itsRef.GetLimits();}
 private:
  const R& itsRef;
};
template <class Derived, Shape S> class IndexableBase;

//----------------------------------------------------
//
//  Make a scalar or constant look like it is indexalble.
//
template <class T, class A, Shape> class Val {};

template <class T, class A> class Val<T,A, VectorShape>
{
 public:
  explicit Val(const T& v, const A& a) : itsVal(v), itsA(a), itsLimits(a.GetLimits()) {};
  T         operator()(index_t) const {return itsVal;}
  size_t    size      (       ) const {return itsLimits.size();} //Apparently itsA is gone so why store it?
  VecLimits GetLimits (       ) const {return itsLimits;} //Apparently itsA.GetLimits() fails at runtime
 private:
  T        itsVal;
  const A& itsA;
  VecLimits itsLimits;
};

//---------------------------------------------------------------------------
//
//  Primary template for the IndexableBase class which provide index iterators.
//  Use partial specialization for each  container shape.
//
template <class Derived, Shape S> class IndexableBase {};

//---------------------------------------------------------------------------
//
//  Primary template for the indexable class.  Use partial specialization
//  for each  container shape.
//
template <class T, class Derived,Store M,Data D,Shape S> class Indexable;




//---------------------------------------------------
//
//  Hold a temporary expression.
//
template <class T, class Expression,Store M,Data D,Shape S> class Xpr;

template <class T, class Expression,Store M,Data D> class Xpr<T,Expression,M,D,VectorShape>
: public Indexable<T,Xpr<T,Expression,M,D,VectorShape>,M,Abstract,VectorShape>
{
 public:
  typedef Indexable<T,Xpr<T,Expression,M,D,VectorShape>,M,Abstract,VectorShape> IndexableT;
  typedef Ref<T,IndexableT,VectorShape> RefT;

  Xpr(Expression e) : itsExp(e) {};
  Xpr(const Xpr& x) : itsExp(x.itsExp) {};
  ~Xpr() {};

  T         operator()(index_t i) const {return itsExp(i);}
  index_t   size      (         ) const {return itsExp.size();}
  VecLimits GetLimits (         ) const {return itsExp.GetLimits();}
 private:
  Expression itsExp;
};

//-----------------------------------------------------------------------------
//
//  Overload lots of Unary functions
//
template <class T, class TR, class A, Store M, Data D, Shape S> inline
auto UnaryFunction(const Indexable <T,A,M,D,S>& a,TR(*f)(const T&))
{
  typedef typename A::RefT RefA;
  typedef XprUnary<T,TR,RefA,S> ExprTLambda;
  return Xpr<TR,ExprTLambda,M,D,S>(ExprTLambda(RefA(a),f));
}
template <class T, class A, Store M, Data D, Shape S> inline
auto operator-(const Indexable <T,A,M,D,S>& a)
{return UnaryFunction<T,T,A,M,D,S>(a,[](const T &x) { return -x; });}

template <class T, class A, Store M, Data D, Shape S> inline
auto operator+(const Indexable <T,A,M,D,S>& a)
{return UnaryFunction<T,T,A,M,D,S>(a,[](const T &x) { return x; });}


#define Op(f) \
template <class T, class A, Store M, Data D, Shape S> inline \
auto f(const Indexable <T,A,M,D,S>& a) \
{return UnaryFunction<T,T,A,M,D,S>(a,[](const T &x) { return f(x); });}

#define OpR(f,TR) \
template <class T, class A, Store M, Data D, Shape S> inline \
auto f(const Indexable <T,A,M,D,S>& a) \
{return UnaryFunction<T,TR,A,M,D,S>(a,[](const T &x) { return f(x); });}

Op(sin  )
Op(cos  )
Op(tan  )
Op(asin  )
Op(acos  )
Op(atan  )
Op(sinh  )
Op(cosh  )
Op(tanh  )
Op(exp  )
Op(log  )
Op(log10)
//Op(pow10)
Op(sqrt)
Op(conj)
OpR(real,double)
OpR(imag,double)
OpR(norm,double)
OpR(arg ,double)
OpR(fabs,double)

#undef Op
#undef OpR

template <class TR, class TA, class TB, class A, class B, Shape S> class XprBinary
{};


template <class TR, class TA, class TB, class A, class B> class XprBinary<TR,TA,TB,A,B,VectorShape>
{
 public:
   XprBinary(const A& a,const B& b,TR(*f)(const TA&,const TB&)) : itsA(a), itsB(b), itsF(f) {};
  ~XprBinary() {};

  TR        operator[](index_t n) const {return itsF(itsA[n],itsB[n]);}
  TR        operator()(index_t n) const {return itsF(itsA(n),itsB(n));}
  size_t    size      (         ) const {return itsA.size();}
  VecLimits GetLimits (         ) const {return itsA.GetLimits();}
 private:
   A itsA;
   B itsB;
   TR(*itsF)(const TA&,const TB&);
};


//----------------------------------------------------------------------------------
//
//  Overload lots of binary operators.  The combinatorics starts to explode here
//      Ob Ob binary
//
template <class TA, class TB, class A, class B, Store MA, Store MB, Data DA, Data DB, Shape S> inline
auto BinaryFunction1(const Indexable<TA,A,MA,DA,S>& a,
                     const Indexable<TB,B,MB,DB,S>& b,
                     typename ReturnType<TA,TB>::RetType (*f)(const TA&, const TB&))
{
//  if (a.size()!=b.size())
//  {
//    std::cout << "indexable.h BinaryFunction1 a=" << a.size() << " b=" << b.size() << std::endl;
//  }
  //assert(a.size()==b.size());
  typedef typename A::RefT RefA;
  typedef typename B::RefT RefB;
  constexpr Data  DR=ReturnData <DA,DB>::RetType ;
  constexpr Store MR=ReturnStore<MA,MB>::RetType;
  typedef typename ReturnType<TA,TB>::RetType TR;
  typedef XprBinary<TR,TA,TB,RefA,RefB,S> ExprTLambda;
  return Xpr<TR,ExprTLambda,MR,DR,S>(ExprTLambda(RefA(a),RefB(b),f));
}

template <class TA, class TB, class A, class B, Store MA, Store MB, Data DA, Data DB, Shape S> inline
auto operator+  (const Indexable<TA,A,MA,DA,S>& a, const Indexable<TB,B,MB,DB,S>& b)
{
  return BinaryFunction1<TA,TB,A,B,MA,MB,DA,DB,S>(a,b,[](const TA &xa,const TB& xb) { return xa+xb; });
}

template <class TA, class TB, class A, class B, Store MA, Store MB, Data DA, Data DB, Shape S> inline
auto operator-  (const Indexable<TA,A,MA,DA,S>& a, const Indexable<TB,B,MB,DB,S>& b)
{
  return BinaryFunction1<TA,TB,A,B,MA,MB,DA,DB,S>(a,b,[](const TA &xa,const TB& xb) { return xa-xb; });
}


template <class TA, class TB, class A, class B, Store MA,Store MB, Data DA, Data DB, Shape S> inline
auto DirectMultiply (const Indexable<TA,A,MA,DA,S>& a, const Indexable<TB,B,MB,DB,S>& b)
{
  return BinaryFunction1<TA,TB,A,B,MA,MB,DA,DB,S>(a,b,[](const TA &xa,const TB& xb) { return xa*xb; });
}
template <class TA, class TB, class A, class B, Store MA,Store MB, Data DA, Data DB, Shape S> inline
auto DirectDivide (const Indexable<TA,A,MA,DA,S>& a, const Indexable<TB,B,MB,DB,S>& b)
{
  return BinaryFunction1<TA,TB,A,B,MA,MB,DA,DB,S>(a,b,[](const TA &xa,const TB& xb) { return xa/xb; });
}


//
//  Ob Scalar binary ops
//
template <class TA, class TB, class A, Store M, Data DA, Shape S> inline
auto BinaryFunction(const Indexable<TA,A,M,DA,S>& a,
                    const TB& b,
                    typename ReturnType<TA,TB>::RetType (*f)(const TA&, const TB&))
{
  typedef typename A::RefT RefA;
  typedef  Val<TB,RefA,S> ValB;
  typedef typename ReturnType<TA,TB>::RetType TR;
  typedef XprBinary<TR,TA,TB,RefA,ValB,S> ExprTLambda;
  return Xpr<TR,ExprTLambda,M,DA,S>(ExprTLambda(RefA(a),ValB(b,RefA(a)),f));
}

// Return type is deduced
template <class TA, class TB, class B, Store M, Data DB, Shape S> inline
auto BinaryFunction(const TA& a,
                    const Indexable<TB,B,M,DB,S>& b,
                    typename ReturnType<TA,TB>::RetType (*f)(const TA&, const TB&))
{
  typedef typename B::RefT RefB;
  typedef  Val<TA,RefB,S> ValA;
  typedef typename ReturnType<TA,TB>::RetType TR;
  typedef XprBinary<TR,TA,TB,ValA,RefB,S> ExprTLambda;
  return Xpr<TR,ExprTLambda,M,DB,S>(ExprTLambda(ValA(a,RefB(b)),RefB(b),f));
}

// Return type is explicit
template <class TR,class TA, class TB, class A, Store M, Data DA, Shape S> inline
auto BinaryFunctionR(const Indexable<TA,A,M,DA,S>& a,
                    const TB& b,
                    TR (*f)(const TA&, const TB&))
{
  typedef typename A::RefT RefA;
  typedef  Val<TB,RefA,S> ValB;
  typedef XprBinary<TR,TA,TB,RefA,ValB,S> ExprTLambda;
  return Xpr<TR,ExprTLambda,M,DA,S>(ExprTLambda(RefA(a),ValB(b,RefA(a)),f));
}
template <class TR, class TA, class TB, class B, Store M, Data DB, Shape S> inline
auto BinaryFunctionR(const TA& a,
                    const Indexable<TB,B,M,DB,S>& b,
                    TR (*f)(const TA&, const TB&))
{
  typedef typename B::RefT RefB;
  typedef  Val<TA,RefB,S> ValA;
  typedef XprBinary<TR,TA,TB,ValA,RefB,S> ExprTLambda;
  return Xpr<TR,ExprTLambda,M,DB,S>(ExprTLambda(ValA(a,RefB(b)),RefB(b),f));
}


#define ObSc1(func,op)\
template <class TA, class A, Store M, Data DA, Shape S> \
inline auto func (const Indexable<TA,A,M,DA,S>& a, const TA& b)\
{  return BinaryFunction<TA,TA,A,M,DA,S>(a,b,[](const TA &xa,const TA& xb) { return op; });}\
template <class TB, class B, Store M, Data DB, Shape S> \
inline auto func  (const TB& a,const Indexable<TB,B,M,DB,S>& b)\
{  return BinaryFunction<TB,TB,B,M,DB,S>(a,b,[](const TB &xa,const TB& xb) { return op; });}

template <class T, class A, class B,Store MA,Store MB, Data DA, Data DB, class L>
inline bool LogicalIII(const Indexable<T,A,MA,DA,VectorShape>& a, const Indexable<T,B,MB,DB,VectorShape>& b,const L& lambda)
{
  assert(a.GetLimits()==b.GetLimits());
  bool ret(true);
    for (index_t i:a.indices())
        ret = ret && lambda(a(i),b(i));
  return ret;
}
template <class T, class A, class B,Store MA,Store MB, Data DA, Data DB, class L>
inline bool LogicalIII(const Indexable<T,A,MA,DA,MatrixShape>& a, const Indexable<T,B,MB,DB,MatrixShape>& b,const L& lambda)
{
  assert(a.GetLimits()==b.GetLimits());
  bool ret(true);
    for (index_t i:a.rows())
        for (index_t j:a.cols())
        ret = ret && lambda(a(i,j),b(i,j));
  return ret;
}

#define ObScBool(func,op)\
template <class T, class A,class B,Store MA,Store MB, Data DA, Data DB,Shape S> \
inline bool func(const Indexable<T,A,MA,DA,S>& a, const Indexable<T,B,MB,DB,S>& b) \
{return LogicalIII(a,b,[](const T& xa,const T&xb){return op;});}

ObScBool(operator==,xa==xb)
ObScBool(operator!=,xa!=xb)

#undef ObScBool


#define ObScMix1(func,op,T1,T2)\
template <class A, Store M, Data DA, Shape S> \
inline auto func (const Indexable<T1,A,M,DA,S>& a, const T2& b) \
{return BinaryFunction<T1,T2,A,M,DA,S>(a,b,[](const T1 &xa,const T2& xb) { return op; });} \
template <class B, Store M, Data DB, Shape S> \
inline auto func  (const T1& a,const Indexable<T2,B,M,DB,S>& b) \
{return BinaryFunction<T1,T2,B,M,DB,S>(a,b,[](const T1 &xa,const T2& xb) { return op; });} \
template <class A, Store M, Data DA, Shape S> \
inline auto func (const Indexable<T2,A,M,DA,S>& a, const T1& b) \
{return BinaryFunction<T2,T1,A,M,DA,S>(a,b,[](const T2 &xa,const T1& xb) { return op; });} \
template <class B, Store M, Data DB, Shape S> \
inline auto func  (const T2& a,const Indexable<T1,B,M,DB,S>& b) \
{return BinaryFunction<T2,T1,B,M,DB,S>(a,b,[](const T2 &xa,const T1& xb) { return op; });}



//----------------------------------------------------------------------------
//
//  Generate lots of expression template functions.
//

ObSc1(operator+ ,xa+xb )
ObSc1(operator- ,xa-xb )
ObSc1(operator* ,xa*xb )
ObSc1(operator/ ,xa/xb )


ObScMix1(operator+ ,xa+xb  ,std::complex<double>,double)
ObScMix1(operator- ,xa-xb  ,std::complex<double>,double)
ObScMix1(operator* ,xa*xb  ,std::complex<double>,double)
ObScMix1(operator/ ,xa/xb  ,std::complex<double>,double)
//ObScMix1(operator+ ,xa+xb  ,std::complex<double>,int)
//ObScMix1(operator- ,xa-xb  ,std::complex<double>,int) //stdcomplex does not support int
//ObScMix1(operator* ,xa*xb  ,std::complex<double>,int)
//ObScMix1(operator/ ,xa/xb  ,std::complex<double>,int)
ObScMix1(operator+ ,xa+xb  ,double,int)
ObScMix1(operator- ,xa-xb  ,double,int)
ObScMix1(operator* ,xa*xb  ,double,int)
ObScMix1(operator/ ,xa/xb  ,double,int)


template <class T, class A, class B, Store MA, Data DA, Store MB, Data DB, Shape S> inline
T Dot(const Indexable<T,A,MA,DA,S>& a,const Indexable<T,B,MB,DB,S>& b)
{
  return Sum(DirectMultiply(a,b));
}


// arrindex.h
//-----------------------------------------------------------------------------
//
//  These macros invoke array index bounds checking if DEBUG is on.
//
#if DEBUG
  #include <cassert>
  #define CHECK(i)\
  assert(i>=0);\
  assert(static_cast<size_t>(i)<size());
#else
  #define CHECK(i)
#endif

template <class T, class Derived,Store M, Shape S> class ArrayIndexable
{
protected: //Can only copy and construct the derived class.
    ArrayIndexable() {};
    ~ArrayIndexable() {};
    ArrayIndexable(const ArrayIndexable&) {};
    ArrayIndexable& operator=(const ArrayIndexable&) {return *this;}

public:
    typedef       T*       iterator;
    typedef const T* const_iterator;

    T operator[](index_t i) const {CHECK(i);return begin()[i];}

    const_iterator begin() const {return static_cast<const Derived*>(this)->priv_begin();}
          iterator begin()       {return static_cast<      Derived*>(this)->priv_begin();}
    const_iterator end  () const {return static_cast<const Derived*>(this)->priv_begin()+size();}
          iterator end  ()       {return static_cast<      Derived*>(this)->priv_begin()+size();}

    size_t  size() const {return static_cast<const Derived*>(this)->size();}

    //  Some overloaded operators
    Derived& operator+=(const T& scalar) {return ArrayAdd(*this,scalar);}
    Derived& operator-=(const T& scalar) {return ArraySub(*this,scalar);}
    Derived& operator*=(const T& scalar) {return ArrayMul(*this,scalar);}
    Derived& operator/=(const T& scalar) {return ArrayDiv(*this,scalar);}

    Derived& operator+=(const ArrayIndexable<T,Derived,M,S>& b) {return ArrayAdd(*this,b);}
    Derived& operator-=(const ArrayIndexable<T,Derived,M,S>& b) {return ArraySub(*this,b);}

    //
    //  Support (index_t i:arr) range iterators over indices
    //
    class index_iterator
    {
    public:
        index_iterator(index_t i) : current{i} {};
        index_iterator operator++(){current++;return (*this);}
        const index_t operator*() const {return current;}
              index_t operator*() {return current;}
        bool operator!=(const index_iterator& b) {return current!=b.current;}
    private:
        index_t current;
    };

    class iterator_proxy
    {
    public:
        iterator_proxy(index_t l, index_t h) : low(l), high(h) {};
        index_iterator begin() const {return low;}
        index_iterator end  () const {return high+1;}
    private:
        index_t low;
        index_t high;
    };

    iterator_proxy arr_indices() const {return iterator_proxy(0,size()-1);}


};

#undef CHECK

template <class T, class A, Shape S> inline T Sum(const ArrayIndexable<T,A,Full,S>& a)
{
  T ret(0);
  for (const T& ai:a) ret+=ai;
  return ret;
}

template <class T,class A, Store M, Shape S> inline void Fill(ArrayIndexable<T,A,M,S>& arr,T value)
{
  for (T& a:arr) a = value; //Sill fails for complex
}

template <class T,class A, Store M, Shape S> inline void FillLinear(ArrayIndexable<T,A,M,S>& arr,T start, T stop)
{
  T del = (stop-start)/(double)(static_cast<A&>(arr).size()-1);
  T val=start;
  for (T& a:arr)
  {
      a = val;
      val+=del;
  }
}

//
//  FE calculus
//
template <class T,class A, Store M, Shape S> inline A Integrate(const ArrayIndexable<T,A,M,S>& arr,T y0=0)
{
  size_t  n=arr.size();
  A ret(n);
  for (index_t i:arr.indices())
  {
    y0+=arr[i];
    ret[i]=y0;
  }
  return ret;
}

template <class T,class A, Store M, Shape S> inline A Differentiate(const ArrayIndexable<T,A,M,S>& arr)
{
  index_t n=arr.size();
  A ret(n);
  auto ri=ret.begin();
  *ri=arr[0]; //Save integration constant in case caller needs it.
  for (index_t i=1;i<n;i++,ri++) *ri=arr[i]-arr[i-1];
  return ret;
}

//------------------------------------------------------------------------------
//
//  IO stuff.  Just dumps the data to the stream, ... thats it.
//
template <class T, class A, Store M, Shape S> inline std::ostream& Write(std::ostream& os,const ArrayIndexable<T,A,M,S>& arr)
{
  assert(os);
  if (StreamableObject::Binary()) for(T b:arr) BinaryWrite(b,os);
  if (StreamableObject::Ascii ()) for(T b:arr) os << b << " ";
  if (StreamableObject::Pretty())
  {
    std::streamsize prec=os.precision();
    std::streamsize wid =os.width();
    os << "{ ";
    for(T b:arr) os << std::setw(wid) << std::setprecision(prec) << b << " ";
    os << "}";
  }
  assert(os);
  return os;
}

template <class T, class A, Store M, Shape S> inline std::istream& Read(std::istream& is,ArrayIndexable<T,A,M,S>& arr)
{
  assert(is);
  if(StreamableObject::Binary())
    for(T& i:arr) BinaryRead(i,is);
  else
    for(T& i:arr) is >> i;

  assert(is);
  return is;
}
//
//  Create assign functions
//
#define OP(NAME,OP) \
template <class T, class Derived,Store M,Shape S> inline \
Derived& Array##NAME (ArrayIndexable<T,Derived,M,S>& a,const T& scalar)\
{\
  for (T& ai:a) ai OP##=scalar;\
  return static_cast<Derived&>(a);\
}\
template <class T,class Derived,Store M,Shape S,class B> inline \
Derived& Array##NAME (ArrayIndexable<T,Derived,M,S>& a, const ArrayIndexable<T,B,M,S>& b)\
{\
    auto bi=b.begin();\
	for (T& ai:a) ai OP##= *bi++;\
	return static_cast<Derived&>(a);\
}\

OP(Add,+)
OP(Sub,-)
OP(Mul,*)
OP(Div,/)

#undef OP


//
//  Logical operators mapped over iterable arrays
//
template <class T, class A, class B, class L, Store M, Shape S>
inline bool LogicalII(const ArrayIndexable<T,A,M,S>& a, const ArrayIndexable<T,B,M,S>& b,const L& lambda)
{
  assert(a.size()==b.size());
  bool ret(true);
  auto bi=b.begin();
  for (const T& ai:a)
  {
    ret = ret && lambda(ai,*bi);
    bi++;
    if (!ret) break;
  }
  return ret;
}

template <class T, class A, class L, Store M, Shape S>
inline bool Logical(const ArrayIndexable<T,A,M,S>& a, const T& b,const L& lambda)
{
  bool ret(true);
  for (const T& i:a) ret = ret && lambda(i,b);
  return ret;
}

template <class T, class A, class L, Store M, Shape S>
inline bool Logical(const T & a, const ArrayIndexable<T,A,M,S>& b,const L& lambda)
{
  bool ret(true);
  for (const T& i:b) ret = ret && lambda(a,i);
  return ret;
}

#define ObScBool(func,op)\
template <class T, class A,class B,Store M, Shape S> \
inline bool func(const ArrayIndexable<T,A,M,S>& a, const ArrayIndexable<T,A,M,S>& b) \
{return LogicalII(static_cast<const A&>(a),static_cast<const B&>(b),[](const T& xa,const T&xb){return op;});}\
template <class T, class A,Store M, Shape S> inline bool func(const ArrayIndexable<T,A,M,S>& a, const T& b)\
{return Logical(a,b,[](const T& xa,const T&xb){return op;});}\
template <class T, class A,Store M, Shape S> inline bool func(const T& a, const ArrayIndexable<T,A,M,S>& b)\
{return Logical(a,b,[](const T& xa,const T&xb){return op;});}\

ObScBool(operator==,xa==xb)
ObScBool(operator!=,xa!=xb)

#undef ObScBool
//
//  isnan isinf logical functions
//
inline constexpr bool isnan(const std::complex<double>& c)
{
    return std::isnan(c.real()) || std::isnan(c.imag());
}
inline bool isinf(const std::complex<double>& c)
{
    return std::isinf(c.real()) || std::isinf(c.imag());
}

template <class T, class D, Store M, Shape S> inline bool isnan(const ArrayIndexable<T,D,M,S>& arr)
{
    bool ret=false;
    for (const T& a : arr) if (isnan(a)) {ret=true;break;}
    return ret;
}

template <class T, class D, Store M, Shape S> inline bool isinf(const ArrayIndexable<T,D,M,S>& arr)
{
    bool ret=false;
    for (const T& a : arr) if (isinf(a)) {ret=true;break;}
    return ret;
}

//------------------------------------------------------------------
//
//  Max/Min functions.
//  TODO: Try using partial specialization of template functions.
//
template <class T, class A, class Op, Store M, Data D,Shape S> class MinMax;

template <class T, class A, class Op, Store M, Shape S> class MinMax<T,A,Op,M,Real,S>
{
public:
    static T apply(const ArrayIndexable<T,A,M,S>& a)
    {
        T ret=a.size()>0 ? a[0] : T(0); // Don't try and read a[0] if there is no data in a!
        for (index_t i:a.arr_indices())
        {
            T ai=a[i];
            if (Op::apply(ai,ret)) ret=ai;
        }
        return ret;
    }
};

template <class T, class A, Store M, Shape S> inline T Min(const ArrayIndexable<T,A,M,S>& a)
{
	return MinMax<T,A,OpLT<T>,M,Real,S>::apply(a);
}

template <class T, class A, Store M, Shape S> inline T Max(const ArrayIndexable<T,A,M,S>& a)
{
	return MinMax<T,A,OpGT<T>,M,Real,S>::apply(a);
}

// vecindex.h
//
// Template specialization provides index iterators for vector shape
//
template <class Derived> class IndexableBase<Derived,VectorShape>
{
    public:
    //
//  Support range based iteration for rows and columns so client code and do
//     for (index_t i:V)
//          {do something with V(i)
//     for (index_t i:V.all())
//          {do something with V[i]
//
  class index_iterator
  {
    public:
        index_iterator(index_t i) : current{i} {};
        index_iterator operator++(){current++;return (*this);}
        const index_t operator*() const {return current;}
        index_t operator*() {return current;}
        friend bool operator!=(const index_iterator& a, const index_iterator& b) {return a.current!=b.current;}
    private:
        index_t current;
  };

    class iterator_proxy
    {
    public:
        iterator_proxy(const VecLimits& lim) : low(lim.Low), high(lim.High) {};
        iterator_proxy(index_t l, index_t h) : low(l), high(h) {};
        index_iterator begin() const {return low;}
        index_iterator end  () const {return high+1;}
    private:
        index_t low;
        index_t high;
    };

    iterator_proxy indices() const {return iterator_proxy(static_cast<const Derived*>(this)->GetLimits());}
    iterator_proxy indices(index_t i) const 
    {
        return iterator_proxy(i,static_cast<const Derived*>(this)->GetLimits().High);
    }
};

//-------------------------------------------------
//
//  template specialization for Vectors's.
//
template <class T, class Derived, Store M, Data D> class Indexable<T,Derived,M,D,VectorShape>
 : public IndexableBase<Derived,VectorShape>
{
 public:

  T  operator()(index_t n) const {return static_cast<const Derived*>(this)->operator()(n);}
  T& operator()(index_t n)       {return static_cast<      Derived*>(this)->operator()(n);}

  size_t    size  () const {return static_cast<const Derived*>(this)->size();}
  VecLimits GetLimits() const {return static_cast<const Derived*>(this)->GetLimits();}

 protected:
  template <class B> void AssignFrom(const Indexable<T,B,Full,Abstract,VectorShape>& b) {VectorAssign(*this,b);}

  explicit Indexable() {};
  ~Indexable() {};
  Indexable& operator=(const Indexable&) {return *this;}
  Indexable(const Indexable&) {};
};

//--------------------------------------------------------------
//
//  Template specialization for abstract vectors.
//
template <class T, class Derived, Store M> class Indexable<T,Derived,M,Abstract,VectorShape>
 : public IndexableBase<Derived,VectorShape>
{
 public:
  explicit Indexable() {};
  ~Indexable() {};

  T operator()(index_t n) const {return static_cast<const Derived*>(this)->operator()(n);}

  size_t    size  () const {return static_cast<const Derived*>(this)->size();}
  VecLimits GetLimits() const {return static_cast<const Derived*>(this)->GetLimits();}

 private:
  Indexable& operator=(const Indexable&);
  Indexable(const Indexable&);
};

//
//  Create assign functions
//
template <class T, class Derived,Store M,Data D,class B,Data DB> inline
void VectorAssign(Indexable<T,Derived,M,D,VectorShape>& a,const Indexable<T,B,M,DB,VectorShape>& b)
{
#ifdef WARN_DEEP_COPY
  std::cerr << "Doing abstract VectorAssign n=" << a.size() << std::endl;
#endif
  assert(a.GetLimits()==b.GetLimits());
  typename Derived::Subscriptor s(a);
  for (index_t i:a.indices()) s(i)=b(i);
}

template <class T, class Derived,Store M,Data D> inline
void VectorAssign(Indexable<T,Derived,M,D,VectorShape>& a,T scalar)
{
#ifdef WARN_DEEP_COPY
  std::cerr << "Doing abstract VectorAssign n=" << a.size() << std::endl;
#endif
  typename Derived::Subscriptor s(a);
  for (index_t i:a.indices()) s(i)=scalar;
}

#define OP(NAME,OP) \
template <class T, class Derived,Store M,Data D> inline \
Derived& Vector##NAME (Indexable<T,Derived,M,D,VectorShape>& a,const T& scalar)\
{\
  typename Derived::Subscriptor s(a); \
  for (index_t i:a) s(i) OP##=scalar;\
  return static_cast<Derived&>(a);\
}\
template <class T,class Derived,Store M,Data D,class B,Data DB> inline \
Derived& Vector##NAME (Indexable<T,Derived,M,D,VectorShape>& a,\
                  const Indexable<T,B,M,DB,VectorShape>& b)\
{\
    typename Derived::Subscriptor s(a); \
    for (index_t i:a) s(i) OP##=b(i);\
	return static_cast<Derived&>(a);\
}\

OP(Add,+)
OP(Sub,-)
OP(Mul,*)
OP(Div,/)

#undef OP


//----------------------------------------------------------------
//
//  Abstract vector specializations for some helper functions.
//

template <class T, class A, Store M> inline
T Sum(const Indexable<T,A,M,Abstract,VectorShape>& a)
{
  T ret(0);
  for (index_t i:a.indices()) ret+=a(i);
  return ret;
}

template <class T, class A, class Op, Store M> class MinMax<T,A,Op,M,Abstract,VectorShape>
{
 public:
  static T apply(const Indexable<T,A,M,Abstract,VectorShape>& a)
  {
    index_t low=a.GetLimits().Low;
    index_t hi =a.GetLimits().High;
    T ret=a(low);
    for (index_t i=low+1;i<=hi;i++)
    {
      T ai=a(i);
      if (Op::apply(ai,ret)) ret=ai;
    }
    return ret;
  }
};

template <class T, class A, Store M> inline T Min(const Indexable<T,A,M,Abstract,VectorShape>& a)
{
	return MinMax<T,A,OpLT<T>,M,Abstract,VectorShape>::apply(a);
}

template <class T, class A, Store M> inline T Max(const Indexable<T,A,M,Abstract,VectorShape>& a)
{
	return MinMax<T,A,OpGT<T>,M,Abstract,VectorShape>::apply(a);
}

// minmax.h
template <class T> inline       T& Max(      T& a,       T& b)  {return a > b ? a : b;}
template <class T> inline const T& Max(const T& a, const T& b)  {return a > b ? a : b;}
template <class T> inline       T& Min(      T& a,       T& b)  {return a < b ? a : b;}
template <class T> inline const T& Min(const T& a, const T& b)  {return a < b ? a : b;}

// tstram.h
template <class A> class TStreamableObject;

template <class A> std::ostream& operator<<(std::ostream& os, const TStreamableObject<A>& o);
template <class A> std::istream& operator>>(std::istream& os,       TStreamableObject<A>& o);

// Allows template based streaming.  IO methods for the parent type A are called directly.
//  ***No virtual dispatch***!  If you need virtual dispatch see the PMStreamableObject class.
template <class A> class TStreamableObject
: public StreamableObject
{
 public:
  friend std::ostream& operator<< <>(std::ostream& os, const TStreamableObject& o);
  friend std::istream& operator>> <>(std::istream& os,       TStreamableObject& o);
//  friend std::ostream& operator<<(std::ostream& os, const TStreamableObject* o) {return os << *o;}
//  friend std::istream& operator>>(std::istream& is,       TStreamableObject* o) {return is >> *o;}
  static A* Factory(std::istream& is);
};


template <class A> inline std::ostream& operator<<(std::ostream& os, const TStreamableObject<A>& o)
{
  o.WriteHeader(os,typeid(A).name());
  return static_cast<const A&>(o).Write(os);
}

template <class A> inline std::istream& operator>>(std::istream& is,TStreamableObject<A>& o)
{
  StreamableObject::Mode current=o.ReadHeader(is,typeid(A).name());
  static_cast<A&>(o).Read (is);
  o.SetOutputMode(current); //Restore to previous state.
  return is;
}

template <class A> inline A* TStreamableObject<A>::Factory(std::istream& is)
{
  //  CheckName(is,typeid(A).name()); Read operation should do this.
  A* ret=new A;
  is >> ret;
  return ret;
}


// cow.h
#ifndef OML_USE_STDVEC
template <class T> class cow_array
{
public:
    cow_array(size_t theSize     );
    cow_array(const cow_array& ca);
    ~cow_array(                   );
    cow_array& operator=(const cow_array& ca);

    const T* begin   () const {                         return itsData;}
          T* begin   ()       {if (*itsOwners>1) COW(); return itsData;}
    const T& operator[](index_t i) const {return itsData[i];}
          T& operator[](index_t i)       {if (*itsOwners>1) COW(); return itsData[i];}

    size_t   size() const {return itsSize;}
    int      GetNumOwners() const {return *itsOwners;}

private:
    void Release();
    void COW    ();

          size_t  itsSize;
          T*      itsData;
    mutable int*  itsOwners;
};


#ifdef DEBUG
  #define CHECK \
  assert(*itsOwners >0); \
  assert(itsData);       \
  assert(itsSize>=0);
#else
  #define CHECK
#endif
#ifdef WARN_DEEP_COPY
  #include <iostream>
#endif

template <class T> inline cow_array<T>::cow_array(size_t theSize)
  : itsSize  (theSize       )
  , itsData  (new T[itsSize])
  , itsOwners(new int(1)    )
  {
#ifdef WARN_DEEP_COPY
  std::cerr << "*** Copy-on-write array allocating size=" << itsSize << " ***" << std::endl;
#endif
    CHECK
  }

template <class T> inline cow_array<T>::cow_array(const cow_array& ca)
  : itsSize  (ca.itsSize  )
  , itsData  (ca.itsData  )
  , itsOwners(ca.itsOwners)
  {
    assert(itsOwners);
    (*itsOwners)++;
    CHECK
  }

template <class T> inline cow_array<T>::~cow_array()
{
  CHECK
  Release();
}

template <class T> cow_array<T>& cow_array<T>::operator=(const cow_array<T>& ca)
{
  CHECK
  if (this!=&ca)
  {
    Release();
    itsSize  =ca.itsSize;
    itsData  =ca.itsData;
    itsOwners=ca.itsOwners;
    (*itsOwners)++;
  }
  CHECK
  return *this;
}

//
//  Release link to current data.
//
template <class T> inline void cow_array<T>::Release()
{
  CHECK
  if (--(*itsOwners)==0)
  {
    delete [] itsData;
    delete    itsOwners;
#ifdef WARN_DEEP_COPY
  std::cerr << "*** Copy-on-write array de-allocating size=" << itsSize << " ***" << std::endl;
#endif
  }
}

//
//  Deep copy for COW (Copy On Write) symantics.
//
template <class T> void cow_array<T>::COW()
{
  CHECK
#ifdef WARN_DEEP_COPY
  std::cerr << "*** Copy-on-write array doing deep copy, size=" << itsSize << " ***" << std::endl;
#endif

  T* newData=new T[itsSize];
  T* source =itsData;
  T* dest   =newData;
  for (;dest<newData+itsSize; dest++,source++) *dest=*source;
  Release();
  itsData=newData;
  itsOwners=new int(1);
  CHECK
}

#undef CHECK

#endif //USE_STD_VEC

// ran250.h
#include <complex>
#include <cassert>
#include <climits>
#ifdef HAVE_STDINT_H
  #include <stdint.h>
#else
  #ifdef HAVE_INTTYPES_H
    #include <inttypes.h>
  #else
    #ifdef  HAVE_SYS_TYPES_H
      #include <sys/types.h>
      typedef  __uint32_t uint32_t;
    #endif
  #endif
#endif

template <class T> union Int2RealUnion;

template <> union Int2RealUnion<float>
{		   	// used to access floats as unsigneds
  float rep;
  uint32_t u;
};

template <> union Int2RealUnion<double>
{		   	// used to access doubles as unsigneds
  double rep;
  uint32_t u[2];
};

template <class T> class Int2Real;

template <> class Int2Real<float>
{
 public:
  Int2Real();
  float convert(uint32_t l) const;
 private:
  Int2RealUnion<float> itsMantissa;
};

template <> class Int2Real<double>
{
 public:
  Int2Real();
  double convert(uint32_t l1,uint32_t l2) const;
 private:
  Int2RealUnion<double> itsMantissa;
};


inline float Int2Real<float>::convert(uint32_t l) const
{
  Int2RealUnion<float> result;
  result.rep = 1.0;
  result.u   |= (l & itsMantissa.u);
  result.rep -= 1.0;
  assert( result.rep < 1.0 && result.rep >= 0);
  return result.rep;
}

inline double Int2Real<double>::convert(uint32_t l1,uint32_t l2) const
{
  Int2RealUnion<double> result;
  result.rep = 1.0;
  result.u[0] |= (l1 & itsMantissa.u[0]);
  result.u[1] |= (l1 & itsMantissa.u[1]);
  result.rep -= 1.0;
  assert( result.rep < 1.0 && result.rep >= 0);
  return result.rep;
}

class TwoTap
{
  unsigned int   Next,Tap1,Tap2,Mask;
  long*          itsArray;

  Int2Real<float > floatConverter;
  Int2Real<double> doubleConverter;

public:
  TwoTap(unsigned int tap1,unsigned int tap2);
 ~TwoTap();
  const char* Name();

  long GetNext()
  {
    ++Next;
    return itsArray[Next&Mask]=
//      (itsArray[(Next-Tap1)&Mask]%(SHRT_MAX-2))
//       *
//      (itsArray[(Next-Tap2)&Mask]%(SHRT_MAX-2));
      itsArray[(Next-Tap1)&Mask]^itsArray[(Next-Tap2)&Mask];
  }

  float  GetNextFloat () {return  floatConverter.convert(GetNext());}
  double GetNextDouble() {return doubleConverter.convert(GetNext(),GetNext());}


};


class FourTap
{
  unsigned int   Next,Tap1,Tap2,Tap3,Tap4,Mask;
  long*          itsArray;

  Int2Real<float > floatConverter;
  Int2Real<double> doubleConverter;


public:
  FourTap(unsigned int tap1,unsigned int tap2,unsigned int tap3,unsigned int tap4);
 ~FourTap();
  const char* Name();

  long GetNext()
  {
    ++Next;
    return itsArray[Next&Mask]=
      itsArray[(Next-Tap1)&Mask]^
      itsArray[(Next-Tap2)&Mask]^
      itsArray[(Next-Tap3)&Mask]^
      itsArray[(Next-Tap4)&Mask];
  }
  float  GetNextFloat () {return  floatConverter.convert(GetNext());}
  double GetNextDouble() {return doubleConverter.convert(GetNext(),GetNext());}
};

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <time.h>
#ifdef HAVE_LIMITS_H
#include <limits.h>
#endif
#include <stdio.h>
#include <cmath>
#include <cassert>
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

//
//	Additive number generator. This method is presented in Volume II
//	of The Art of Computer Programming by Knuth. I've coded the algorithm
//	and have added the extensions by Andres Nowatzyk of CMU to randomize
//	the result of algorithm M a bit	by using an LCG & a spatial
//	permutation table.
//
//	The version presented uses the same constants for the LCG that Andres
//	uses (chosen by trial & error). The spatial permutation table is
//	the same size (it's based on word size). This is for 32-bit words.
//
//	The ``auxillary table'' used by the LCG table varies in size, and
//	is chosen to be the the smallest power of two which is larger than
//	twice the size of the state table.
//

class ACG
{
  uint32_t initialSeed;	// used to reset generator
  int initialTableEntry;

  uint32_t *state;
  uint32_t *auxState;
  short stateSize;
  short auxSize;
  uint32_t lcgRecurr;
  short j;
  short k;
public:
  ACG(uint32_t seed = 0, int size = 55);
  ~ACG();
  uint32_t asLong();
  void reset();
};

double INorm=1.0/LONG_MAX;
double LNorm=1.0/LONG_MAX;

TwoTap::TwoTap(unsigned int tap1,unsigned int tap2)
  : Next(0)
  , Tap1(tap1)
  , Tap2(tap2)
  , Mask(0)
  , itsArray(0)
{
  unsigned int max=Max(Tap1,Tap2);
  unsigned int n=(unsigned int)pow(2.0,(int)floor(log(max)/log(2.0)+1.0));
  itsArray=new long[n];
  Mask=n-1;
  ACG Rand(time(0)); //Use this guy to boot strap R(250,103).
  for (unsigned int i=0;i<n;i++) itsArray[i]=Rand.asLong();
}

TwoTap::~TwoTap()
{
  delete [] itsArray;
}

const char* TwoTap::Name()
{
    std::ostringstream os;
    os << "R(" << Tap1 << "," << Tap2 << ",*)";
    return os.str().c_str();
}

FourTap::FourTap(unsigned int tap1,unsigned int tap2,unsigned int tap3,unsigned int tap4)
  : Next(0)
  , Tap1(tap1)
  , Tap2(tap2)
  , Tap3(tap3)
  , Tap4(tap4)
  , Mask(0)
  , itsArray(0)
{
  unsigned int max=Max(Max(Tap1,Tap2),Max(Tap3,Tap4));
  unsigned int n=(unsigned int)pow(2.0,(int)floor(log(max)/log(2.0)+1.0));
  itsArray=new long[n];
  Mask=n-1;
  ACG Rand(time(0)); //Use this guy to boot strap R(250,103).
  for (unsigned int i=0;i<n;i++) itsArray[i]=Rand.asLong();
}

FourTap::~FourTap()
{
  delete [] itsArray;
}

const char* FourTap::Name()
{
    std::ostringstream os;
    os << "R(" << Tap1 << "," << Tap2 << "," << Tap3 << "," << Tap4 << ")";
    return os.str().c_str();
}


FourTap GlobalRandomNumberGenerator(471,1586,6988,9869);
//TwoTap GlobalRandomNumberGenerator(1063,1279);
//TwoTap GlobalRandomNumberGenerator(103,250);

//
//	The following is a hack that I attribute to
//	Andres Nowatzyk at CMU. The intent of the loop
//	is to form the smallest number 0 <= x < 1.0,
//	which is then used as a mask for two longwords.
//	this gives us a fast way way to produce double
//	precision numbers from longwords.
//
//	I know that this works for IEEE and VAX floating
//	point representations.
//
//	A further complication is that gnu C will blow
//	the following loop, unless compiled with -ffloat-store,
//	because it uses extended representations for some of
//	of the comparisons. Thus, we have the following hack.
//	If we could specify #pragma optimize, we wouldn't need this.
//
Int2Real<double>::Int2Real()
{
  assert (sizeof(double) == 2 * sizeof(uint32_t));
  Int2RealUnion<double> t;
#if _IEEE == 1

  t.rep = 1.5;
  if ( t.u[1] == 0 )
  {		// sun word order?
    t.u[0] = 0x3fffffff;
    t.u[1] = 0xffffffff;
  }
  else
  {
    t.u[0] = 0xffffffff;	// encore word order?
    t.u[1] = 0x3fffffff;
  }

#else
  volatile double x = 1.0; // volatile needed when fp hardware used,
  // and has greater precision than memory doubles
  double y = 0.5;
  do
  {			    // find largest fp-number < 2.0
    t.rep = x;
    x += y;
    y *= 0.5;
  }
  while (x != t.rep && x < 2.0);

#endif
  // set doubleMantissa to 1 for each doubleMantissa bit
  itsMantissa.rep = 1.0;
  itsMantissa.u[0] ^= t.u[0];
  itsMantissa.u[1] ^= t.u[1];
}

Int2Real<float>::Int2Real()
{
  assert (sizeof(float) == sizeof(uint32_t));
  Int2RealUnion<float> t;
#if _IEEE == 1
  t.u = 0x3fffffff;
#else
  t.rep=0.0; //avoid valgrind warning.
  volatile float x = 1.0; // volatile needed when fp hardware used,
  // and has greater precision than memory floats
  float y = 0.5;
  do
  {			    // find largest fp-number < 2.0
    t.rep = x;
    x += y;
    y *= 0.5;
  }
  while (x != t.rep && x < 2.0);
#endif
  // set singleMantissa to 1 for each singleMantissa bit
  itsMantissa.rep = 1.0;
  itsMantissa.u ^= t.u;

}




extern double  INorm;
extern double  LNorm;

template <class T> T      OMLRand();
template <class T> T      OMLRandPos();
template <class T> double OMLRandScale(T max);

template <> inline int    OMLRand<int>   () {return GlobalRandomNumberGenerator.GetNext();}
template <> inline long   OMLRand<long>  () {return GlobalRandomNumberGenerator.GetNext();}
template <> inline float  OMLRand<float> ()
{
  return GlobalRandomNumberGenerator.GetNextFloat();
}
template <> inline double OMLRand<double>()
{
  return GlobalRandomNumberGenerator.GetNextDouble();
}

template <> inline std::complex<double> OMLRand<std::complex<double> >()
{
  return std::complex<double>(OMLRand<double>(),OMLRand<double>());
}

template <> inline int    OMLRandPos<int>   () {return OMLRand<int> ()&0x7fffffff;}
template <> inline long   OMLRandPos<long>  () {return OMLRand<long>()&0x7fffffff;}
template <> inline float  OMLRandPos<float> () {return OMLRand<float>();}
template <> inline double OMLRandPos<double>() {return OMLRand<double>();}
template <> inline std::complex<double> OMLRandPos<std::complex<double> >()
{
  return std::complex<double>(OMLRandPos<double>(),OMLRandPos<double>());
}


template <> inline double OMLRandScale<int>   (int    max) {return INorm*max;}
template <> inline double OMLRandScale<long>  (long   max) {return INorm*max;}
template <> inline double OMLRandScale<float> (float  max) {return max;}
template <> inline double OMLRandScale<double>(double max) {return max;}
template <> inline double OMLRandScale<std::complex<double> >(std::complex<double> max) {return real(max);}

// random.h
/*! \file random.h
  \brief Routines for filling OML containers with random numbers.

  Random number are supplied by a 2 tap xor random number generator. Very fast!
 */
//! Fill with random numbers, range = (-1,1).
template <class T,class A,Store M, Shape S> void FillRandom(ArrayIndexable<T,A,M,S>& arr)
{
  for (T& i:arr) i=OMLRand<T>();
}

//! Fill with positive random numbers, range = [0,1).
template <class T,class A,Store M, Shape S> void FillRandomPositive(ArrayIndexable<T,A,M,S>& arr)
{
  for (T& i:arr) i=OMLRandPos<T>();
}

//! Fill with positive random numbers, range = [-max,max).
template <class T,class A,Store M, Shape S> void FillRandom(ArrayIndexable<T,A,M,S>& arr,T max)
{
  double scale=OMLRandScale<T>(max);
  for (T& i:arr) i = (T)(OMLRand<T>() * scale);
}

//! Fill with positive random numbers, range = [0,max).
template <class T,class A,Store M, Shape S> void FillRandomPositive(ArrayIndexable<T,A,M,S>& arr,T max)
{
  double scale=OMLRandScale<T>(max);
  for (T& i:arr) i = (T)(OMLRandPos<T>() * scale);
}

// vector3d.h
const double pi=acos(-1.0);

//-----------------------------------------------------------------------------
/*! \class Vector3D vector3d.h oml/vector3d.h
  \brief Very light weight 3D %Vector class with lots of overloaded operators.

  Structure for real space three vectors.  Most algeabraic operators have
  been overloaded.  So Vector3Ds should behave like any intrinsic data type.
  Lots of overloaded operators form 3 element vectors.
  - +, -, +=, -=, ==, !=
  - \a a \c % \a b is the vector cross product.
  - \c * is a dot product.
  - /, *=, /= for scalers.
  - <, >, <=, >= all compare magnitudes.
  - Unary + and -
  - \a a| \a b returns the angle in radians between \a a and \a b.
  - \a a|| \a b returns the angle in degrees between \a a and \a b.
  - ! \a a returns the magnitude of \a a.
  - ~ \a a returns a unit vector in the direction of \a a.
  - For complex values \c conj, \c real, \c imag and \c norm are defined.

  This class is much more efficient that \c Vector<double>(3) would be.
  You will need to include \c io3d.h to get \c op<< and \c op>> for IO.
  \nosubgrouping
*/
template <class T> class Vector3D
{
 public:
  /*! \name Constructors/Assignment*/
  //@{

  //! Default contructor, \c x=y=z=0.
  Vector3D(                 ): x(  0),y(  0),z(  0) {};
  //! Contruct from individual components.
  Vector3D(T _x,T _y,T _z   ): x( _x),y( _y),z( _z) {};
  //! Copy constructor.
  Vector3D(const Vector3D& v): x(v.x),y(v.y),z(v.z) {};
  //! Construct from another data type.
  template <class T1> Vector3D(T1 _x,T1 _y,T1 _z    ) : x( _x),y( _y),z( _z) {}
  //! Construct form another Vector3D type.
  template <class T1> Vector3D(const Vector3D<T1>& v) : x(v.x),y(v.y),z(v.z) {}

  //! Assign
  Vector3D& operator =(const Vector3D& v) {x=v.x;y=v.y;z=v.z;return *this;}
  //! Assign form another Vector3D type.
  template <class T1> Vector3D& operator =(const Vector3D<T1>& v) {x=v.x;y=v.y;z=v.z;return *this;}
  //@}

 ~Vector3D() {};

  //! Element access
  const T& operator()(index_t i) const {return (&x)[i-1];}
  //! Element access
  T& operator()(index_t i)       {return (&x)[i-1];}

 


  /*! \name Coordinates*/
  //@{
  T x; //!< \a x coordinate.
  T y; //!< \a y coordinate.
  T z; //!< \a z coordinate.
  //@}
};


//-----------------------------------------------------------------------------
//
//  Binary algeabra.
//


template <class T1, class T2> inline
auto operator +(const Vector3D<T1>& a,const Vector3D<T2>& b)
{
  typedef typename ReturnType<T1,T2>::RetType TR;
  return Vector3D<TR>(a.x+b.x,a.y+b.y,a.z+b.z);
}

template <class T1, class T2> inline
auto operator -(const Vector3D<T1>& a,const Vector3D<T2>& b)
{
  typedef typename ReturnType<T1,T2>::RetType TR;
  return Vector3D<TR>(a.x-b.x,a.y-b.y,a.z-b.z);
}

template <class T1, class T2> inline
auto operator *(const Vector3D<T1>& a,const Vector3D<T2>& b)
{
  //typedef typename ReturnType<T1,T2>::RetType TR;
  return (a.x*b.x+a.y*b.y+a.z*b.z);
}

template <class T> inline
Vector3D<T> operator *(const Vector3D<T>& a,const T F)
{
  return Vector3D<T>(a.x*F,a.y*F,a.z*F);
}

template <class T> inline
Vector3D<T> operator *(const T F,const Vector3D<T>& a)
{
  return Vector3D<T>(a.x*F,a.y*F,a.z*F);
}

template <class T> inline
Vector3D<T> operator /(const Vector3D<T>& a,const T F)
{
  return Vector3D<T>(a.x/F,a.y/F,a.z/F);
}

template <class T1, class T2> inline
auto operator %(const Vector3D<T1>& a,const Vector3D<T2>& b)
{
  typedef typename ReturnType<T1,T2>::RetType TR;
  return Vector3D<TR>( a.x%b.x, a.y%b.y, a.z%b.z );
}

template <class T1, class T2> inline
auto Cross(const Vector3D<T1>& a,const Vector3D<T2>& b)
{
  typedef typename ReturnType<T1,T2>::RetType TR;
  return Vector3D<TR>( a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x );
}

//
//  A operator= B overloads for binary operators
//

template <class T1, class T2> inline
Vector3D<T1>& operator +=(Vector3D<T1>& a,const Vector3D<T2>& b)
{
  a.x+=b.x; a.y+=b.y; a.z+=b.z;
  return a;
}

template <class T1, class T2> inline
Vector3D<T1>& operator -=(Vector3D<T1>& a,const Vector3D<T2>& b)
{
  a.x-=b.x; a.y-=b.y; a.z-=b.z;
  return a;
}

template <class T> inline
Vector3D<T>& operator *=(Vector3D<T>& a,const T F)
{
  a.x*=F;a.y*=F;a.z*=F;
  return a;
}

template <class T> inline
Vector3D<T>& operator /=(Vector3D<T>& a,const T F)
{
  a.x/=F;a.y/=F;a.z/=F;
  return a;
}

//-----------------------------------------------------------------------------
//
//  Relational operators.
//
template <class T> inline
bool operator ==(const Vector3D<T>& a,const Vector3D<T>& b)
{
  return ((a.x==b.x)&&(a.y==b.y)&&(a.z==b.z));
}

template <class T> inline
bool operator !=(const Vector3D<T>& a,const Vector3D<T>& b)
{
  return !(a==b);
}

template <class T> inline
bool operator > (const Vector3D<T>& a,const Vector3D<T>& b)
{
  return (norm(a) > norm(b));
}

template <class T> inline
bool operator < (const Vector3D<T>& a,const Vector3D<T>& b)
{
  return (norm(a) < norm(b));
}

template <class T> inline
bool operator >=(const Vector3D<T>& a,const Vector3D<T>& b)
{
  return (norm(a) >= norm(b));
}

template <class T> inline
bool operator <=(const Vector3D<T>& a,const Vector3D<T>& b)
{
  return (norm(a) <= norm(b));
}

//-----------------------------------------------------------------------------
//
//  Unary operators
//
template <class T> inline
Vector3D<T>  operator -(const Vector3D<T>& a)
{
  return Vector3D<T>(-a.x,-a.y,-a.z);
}

template <class T> inline
Vector3D<T>  operator +(const Vector3D<T>& a)
{
  return a;
}

//-----------------------------------------------------------------------------
//
//  Angle between two vectors (Radians).
//
template <class T> inline
T angle(const Vector3D<T>& a,const Vector3D<T>& b)
{
  return (T)acos( normalize(a) * normalize(b) );
}

//-----------------------------------------------------------------------------
//
//  Angle between two vectors (Degrees).
//
template <class T> inline
T angle_degrees(const Vector3D<T>& a,const Vector3D<T>& b)
{
  return static_cast<T>(acos( normalize(a) * normalize(b) )/pi*180.0);
}

//-----------------------------------------------------------------------------
//
//  Magnitude and normalize.
//
template <class T> inline
T norm(const Vector3D<T>& a)
{
  return (T)sqrt(a*a);
}

template <class T> inline
Vector3D<T>  normalize(const Vector3D<T>& a)  //normalize.
{
  return a/norm(a);
}

//--------------------------------------------------------------------
//
//  Specialized templates for complex data types.
//
template <class T> inline Vector3D<std::complex<T> > conj(const Vector3D<std::complex<T> >& v)
{
  return Vector3D<std::complex<T> >(conj(v.x),conj(v.y),conj(v.z));
}

template <class T> inline Vector3D<T> real(const Vector3D<std::complex<T> >& v)
{
  return Vector3D<T>(real(v.x),real(v.y),real(v.z));
}


template <class T> inline Vector3D<T> imag(const Vector3D<std::complex<T> >& v)
{
  return Vector3D<T>(imag(v.x),imag(v.y),imag(v.z));
}





// vector.h
/*! \class Vector vector.h oml/vector.h
  \brief Numerical container with FORTRAN array symatics.

  Vectors have element indexes ranging from 1...n (by default) and are indexed using the
  (i) syntax, just like FORTRAN arrays.  The base index can also be changed just like in FORTRAN.
  Vectors support efficient copy on write (COW) symmantices, just like Arrays.

  \b Math:

  The special operators for Vector are
  - \c operator* does a dot product (not a direct multiply).
  - for complex value \c V*V will do V-dot-conj(V).
  - use \c DirectMultiply(V,V) to get a Vector of element by element products.
  - \c operator/ is is only allowed for scalars.
  - \c Magnitude(V) and \c !V both return the magnitude of V.
  - \c Normalize(V) returns a unit vector with the same direction as V.

  \nosubgrouping
*/

enum class FillType {Zero,Random,Unit};

template <class T> class Vector
  : public Indexable<T,Vector<T>,Full,Real,VectorShape>
  , public ArrayIndexable<T,Vector<T>,Full,VectorShape>
  , public TStreamableObject<Vector<T> >
{
 public:
  typedef Indexable<T,Vector<T>,Full,Real,VectorShape> IndexableT;
  typedef ArrayIndexable <T,Vector<T>,Full     ,VectorShape> IterableT;
  typedef Ref<T,IndexableT,VectorShape> RefT;
  typedef Ref<T,IterableT ,VectorShape> RefAT;

  /*! \name Constructors/Assignment
  Copy constructor and op=(Vector) are automaically supplied by the compiler.
  */
  //@{
  //!< Vector with size=0;
           Vector() : Vector<T>(VecLimits(1,0)) {};
  //!<  All elements are un-initiated, low index is 1.
  explicit Vector(size_t  size) : Vector<T>(VecLimits(1,size)) {};
  explicit Vector(size_t  size, FillType ft) : Vector<T>(VecLimits(1,size),ft) {};
  explicit Vector(size_t  size, const T&  fillValue) : Vector<T>(VecLimits(1,size),fillValue) {};
  //!<  Specify lower and upper index.
//  explicit Vector(index_t l,index_t h) : Vector<T>(VecLimits(l,h)) {};
  explicit Vector(index_t l,index_t h, FillType ft) : Vector<T>(VecLimits(l,h),ft) {};
  explicit Vector(index_t l,index_t h, const T&  fillValue) : Vector<T>(VecLimits(l,h),fillValue) {};
  Vector(const VecLimits&); //!<  Specify lower and upper index.
  Vector(const VecLimits&, FillType ft); //!<  Specify lower and upper index.
  Vector(const VecLimits&, const T&  fillValue       ); //!<  Specify lower and upper index.
  //! Allows construction from an expression template.
  template <class B,Data D> Vector(const Indexable<T,B,Full,D,VectorShape>&);
  //! Allows assignment from an expression template.
  template <class B,Data D> Vector& operator=(const Indexable<T,B,Full,D,VectorShape>&);

  Vector(const Vector& m);
  Vector& operator=(const Vector&);

#ifdef OML_MOVE_OPS
  Vector(Vector&& m);
  Vector& operator=(Vector&&);
#endif

  VecLimits ReBase(int low);
  VecLimits ReBase(const VecLimits& lim);
  //@}
  void Fill(FillType);
  void Fill(const T& fillValue);
  void FillRandom();

  std::ostream& Write(std::ostream&) const;
  std::istream& Read (std::istream&)      ;
  // Need to disambiguate from expression version.
  friend std::ostream& operator<<(std::ostream& os,const Vector& a)
  {
    return os << static_cast<const TStreamableObject<Vector<T> >& >(a);
  }


  /*! \name Subscripting operators
    FORTRAN style 1-D array subscripting.
    If DEBUG is defined, every index will be checked that it is in range.
    For fast write access with op() make a subscriptor, then the COW check is only done once
    during construction of the Subscriptor.
  */
  //@{
  //! const element acces operator, fast and \e cannot trigger a COW operation.
  T  operator()(index_t) const;
  //! non-const version can trigger a COW operation, and checks for this with every access.
  T& operator()(index_t)      ;
  //@}

  size_t    size     () const; //!<Returns number elements in the Vector.
  VecLimits GetLimits() const; //!<Returns lower and upper limits in a structure.
  index_t   GetLow   () const; //!<Returns lower index.
  index_t   GetHigh  () const; //!<Returns upper index.

  /*! \name Resize functions
    The data can be optionally preserved.
  */
  //@{
  void SetLimits(const VecLimits&, bool preserve=false); //!<Resize from new limits.
  void SetLimits(size_t          , bool preserve=false); //!<Resize from new size.
  void SetLimits(index_t,index_t , bool preserve=false); //!<Resize from new limits.
  //@}

  //! Does V(i)=V(index[i]) for i=low...high.  Used for sorting.
  void   ReIndex(const std::vector<index_t>&);

  /*! \name Sub-vector functions  */
  //@{
  Vector SubVector(const VecLimits&) const; //!< From limits.
  Vector SubVector(size_t          ) const; //!< First n elements.
  Vector SubVector(index_t,index_t ) const; //!< From limits.
  //@}

  /*! \name Iterators.
    Iterators should be STL compatable.
   */
  //@{
  //! Read only iterator.
  typedef typename IterableT::const_iterator  const_iterator;
  //! Read/write iterator.
  typedef typename IterableT::iterator iterator;
  //@}

#if DEBUG
  #define CHECK(i) assert(itsLimits.CheckIndex(i))
#else
  #define CHECK(i)
#endif
  class Subscriptor
  {
   public:
    Subscriptor(Indexable<T,Vector,Full,Real,VectorShape>& a)
      : itsLimits(a.GetLimits())
      , itsPtr(static_cast<Vector*>(&a)->priv_begin()-itsLimits.Low)
      {assert(itsPtr);}

    T& operator()(index_t i) {CHECK(i);return itsPtr[i];}

   private:
    VecLimits itsLimits;
    T*        itsPtr;
  };
#undef CHECK

 private:
  friend class      Indexable<T,Vector,Full,Real,VectorShape>;
  friend class ArrayIndexable<T,Vector,Full     ,VectorShape>;
  friend class Subscriptor;

  const T* priv_begin() const {return &*itsData.begin();} //Required by iterable.
        T* priv_begin()       {return &*itsData.begin();} //Required by iterable.
  void  Check () const; //Check internal consistency between limits and cow.

  VecLimits    itsLimits; //Manages the upper and lower vector limits.
#ifdef OML_USE_STDVEC
  std::vector<T> itsData;
#else
  cow_array<T> itsData;   //Copy-On-Write array for the data.
#endif
};


//-------------------------------------------------------------------------
//
//  Some special operators only for vectors.
//
template <class T, class A, Store M, Data D> inline
std::ostream& operator<<(std::ostream& os,const Indexable<T,A,M,D,VectorShape>& a)
{
  return os << Vector<T>(a);
}


// The define op* as a vector dot (inner) product.
template <class TA, class TB, class A, class B, Store MA, Store MB, Data DA, Data DB> inline
auto operator*(const Indexable<TA,A,MA,DA,VectorShape>& a, const Indexable<TB,B,MB,DB,VectorShape>& b)
{
    typedef typename ReturnType<TA,TB>::RetType TR;
    TR ret(0);
    for (index_t k:a.indices()) ret+=a(k)*b(k);
    return ret;
}

// Calculate the magnitue of a vector.
template <class T, class A, Store M, Data D> inline
T Magnitude(const Indexable<T,A,M,D,VectorShape>& a)
{
  return sqrt(a*a);
}

// Overload op! for magnitue.
template <class T, class A, Store M, Data D> inline
T operator!(const Indexable<T,A,M,D,VectorShape>& a)
{
  return Magnitude(a);
}

// Rescale so that v*v==1.0 .
template <class T, class A, Store M, Data D> inline
void Normalize(Indexable<T,A,M,D,VectorShape>& v)
{
  v/=!v;
}
//-----------------------------------------------------------------------------
//
//  Macro, which expands to an index checking function call,
//  when DEBUG is on
//
#if DEBUG
  #define CHECK(i) itsLimits.CheckIndex(i)
#else
  #define CHECK(i)
#endif

template <class T> inline T Vector<T>::operator()(index_t i) const
{
  CHECK(i);
  return itsData[itsLimits.Offset(i)];
}

template <class T> inline T& Vector<T>::operator()(index_t i)
{
  CHECK(i);
  return itsData[itsLimits.Offset(i)];
}

#undef CHECK

#ifdef DEBUG
#define CHECK Check()
#else
#define CHECK
#endif

template <class T> inline Vector<T>::Vector(const VecLimits& lim,FillType ft)
  : Vector<T>(lim)
  {
    Fill(ft);
  }

template <class T> inline Vector<T>::Vector(const VecLimits& lim,const T& fillValue)
  : Vector<T>(lim)
  {
     Fill(fillValue);
  }

//
//  All other constructors should delegate to this one.
//
template <class T> inline Vector<T>::Vector(const VecLimits& lim)
  : itsLimits(lim          )
  , itsData  (lim.size())
  {
    CHECK;
  }

template <class T> inline Vector<T>::Vector(const Vector<T>& v)
  : itsLimits(v.itsLimits)
  , itsData  (v.itsData)
  {}

template <class T> inline Vector<T>& Vector<T>::operator=(const Vector<T>& v)
{
  itsLimits=v.itsLimits;
  itsData  =v.itsData;
//  std::cout << "Vector<T> move op=" << std::endl;
  return *this;
}
#ifdef OML_MOVE_OPS

template <class T> inline Vector<T>::Vector(Vector<T>&& v)
  : itsLimits(std::move(v.itsLimits))
  , itsData  (std::move(v.itsData))
  {
//    std::cout << "Vector<T> move constructor m.itsData.size()=" << m.itsData.size() << std::endl;
  }

template <class T> inline Vector<T>& Vector<T>::operator=(Vector<T>&& v)
{
  itsLimits=std::move(v.itsLimits);
  itsData  =std::move(v.itsData);
//  std::cout << "Vector<T> move op=" << std::endl;
  return *this;
}
#endif

template <class T> inline VecLimits Vector<T>::ReBase(int low)
{
    VecLimits oldLimits=itsLimits;
    itsLimits.ReBase(low);
    return oldLimits;
}

template <class T> inline VecLimits Vector<T>::ReBase(const VecLimits& lim)
{
    VecLimits oldLimits=itsLimits;
    itsLimits.ReBase(lim.Low);
    return oldLimits;
}

template <class T> inline void Vector<T>::Fill(FillType ft)
{
    switch (ft)
    {
        case (FillType::Zero)   : Fill(T(0));break;
        case (FillType::Random) : ::FillRandom(*this);break;
        case (FillType::Unit)   : Fill(T(1));break;
    }
}

template <> inline Vector3D<double> OMLRand<Vector3D<double> >()
{
  return Vector3D<double>(OMLRand<double>(),OMLRand<double>(),OMLRand<double>());
}

template <>  inline void Vector<Vector3D<double> >::Fill(FillType ft)
{
    switch (ft)
    {
        case (FillType::Zero)   : Fill(Vector3D<double>(double(0),double(0),double(0)));break;
        case (FillType::Random) : ::FillRandom(*this);break;
        case (FillType::Unit)   : Fill(Vector3D<double>(double(1),double(1),double(1)));break;
    }
}

template <class T> inline void Vector<T>::Fill(const T& fillValue)
{
    ::Fill(*this,fillValue);
}

template <class T> inline void Vector<T>::FillRandom()
{
    ::FillRandom(*this);
}


template <class T> inline size_t  Vector<T>::size() const
{
  return GetLimits().size();
}



template <class T> template <class B,Data D> inline
Vector<T>::Vector(const Indexable<T,B,Full,D,VectorShape>& v)
  : itsLimits(v.GetLimits())
  , itsData  (itsLimits.size())
  {
    this->AssignFrom(v); //Choose op[i] or op(i) depending on whether v is abstract.
    CHECK;
  }

#undef CHECK

template <class T> inline  VecLimits Vector<T>::GetLimits() const
{
  return itsLimits;
}

template <class T> inline index_t Vector<T>::GetLow() const
{
  return itsLimits.Low;
}

template <class T> inline index_t Vector<T>::GetHigh() const
{
  return itsLimits.High;
}

template <class T> inline void Vector<T>::SetLimits(size_t size,bool preserve)
{
  SetLimits(VecLimits(size),preserve);
}

template <class T> inline void Vector<T>::SetLimits(index_t l,index_t h,bool preserve)
{
  SetLimits(VecLimits(l,h),preserve);
}

template <class T> inline Vector<T> Vector<T>::SubVector(size_t size) const
{
  return SubVector(VecLimits(size));
}

template <class T> inline Vector<T> Vector<T>::SubVector(index_t l,index_t h) const
{
  return SubVector(VecLimits(l,h));
}

template <class T> template <class B, Data D> inline
Vector<T>& Vector<T>::operator=(const Indexable<T,B,Full,D,VectorShape>& v)
{
  if (size()==0) SetLimits(v.GetLimits());
  this->AssignFrom(v); //Choose op[i] or op(i) depending on whether v is abstract.
  return *this;
}

#ifdef DEBUG
#define CHECK Check()
#else
#define CHECK
#endif



//-----------------------------------------------------------------------------
//
//  Changing size and/or limits.
//
template <class T> void Vector<T>::SetLimits(const VecLimits& theLimits, bool preserve)
{
  theLimits.Check();
  if (itsLimits!=theLimits)
  {
    if (preserve)
    {
      Vector<T> dest(theLimits);

      index_t low =Max(GetLow() ,theLimits.Low );   //Limits of old and new
      index_t high=Min(GetHigh(),theLimits.High);   //data overlap.

      Subscriptor source(*this);         //Destination subscriptor.

      for (index_t i=low;i<=high;i++) dest(i)=source(i); //Transfer any overlaping data.
      *this=dest;
    }
    else
    {
      *this=Vector<T>(theLimits);
    }
  }
  CHECK;
}

//  Concatenate two Vectors into a new Vector.
template <class T> Vector<T> operator&(const Vector<T>& a, const Vector<T>& b)
{
  typedef typename Vector<T>::const_iterator CI;
  VecLimits newlim(a.GetLimits().Low,a.GetLimits().Low+a.size()+b.size()-1);
  Vector<T> ret(newlim);
  typename Vector<T>::iterator i=ret.begin();
  for (CI ab=a.begin();ab!=a.end();ab++,i++) *i=*ab;
  for (CI bb=b.begin();bb!=a.end();bb++,i++) *i=*bb;
  return ret;
}

template <class T> void Vector<T>::ReIndex(const std::vector<index_t>& index)
{
  assert(size()==index.size());

  std::vector<index_t>::const_iterator b=index.begin();
  Vector<T>                      dest(GetLimits());
  iterator                       i=dest.begin();
  for (;b!=index.end();b++,i++) *i=(*this)[*b];
  *this=dest;
}

template <class T> Vector<T> Vector<T>::SubVector(const VecLimits& lim) const
{
  assert(GetLimits().Low  <= lim.Low );
  assert(GetLimits().High >= lim.High);
  Vector<T> ret(lim);
  Subscriptor r(ret);
  for (index_t i=ret.GetLimits().Low;i<=ret.GetLimits().High;i++) r(i)=(*this)(i);
  return ret;
}


#if DEBUG
//-----------------------------------------------------------------------------
//
//  Internal consistency check.
//
template <class T> void Vector<T>::Check() const
{
  assert(itsLimits.Check());
  assert(itsLimits.size()==itsData.size());
}
#endif //DEBUG
//-----------------------------------------------------------------------------
//
//  IO
//
template <class T> std::ostream& Vector<T>::Write(std::ostream& os) const
{
  assert(os);
  if (!this->Pretty())
    os << GetLimits();
   else
  {
    std::streamsize prec=os.precision();
    std::streamsize wid =os.width();
    os << GetLimits();
    os << std::setw(wid) << std::setprecision(prec);
  }
  return ::Write(os,*this);
}

template <class T> std::istream& Vector<T>::Read(std::istream& is)
{
  assert(is);
  VecLimits lim;
  is >> lim;
  if (size()==0)
    SetLimits(lim);
   else
    assert(GetLimits()==lim);

  ::Read(is,*this);
  CHECK;
  assert(is);
  return is;
}



#undef CHECK

// template class Vector<double>;
// template class Vector<std::complex<double>>;

} //export block


//
//	This is an extension of the older implementation of Algorithm M
//	which I previously supplied. The main difference between this
//	version and the old code are:
//
//		+ Andres searched high & low for good constants for
//		  the LCG.
//
//		+ theres more bit chopping going on.
//
//	The following contains his comments.
//
//	agn@UNH.CS.CMU.EDU sez..
//
//	The generator below is based on 2 well known
//	methods: Linear Congruential (LCGs) and Additive
//	Congruential generators (ACGs).
//
//	The LCG produces the longest possible sequence
//	of 32 bit random numbers, each being unique in
//	that sequence (it has only 32 bits of state).
//	It suffers from 2 problems: a) Independence
//	isnt great, that is the (n+1)th number is
//	somewhat related to the preceding one, unlike
//	flipping a coin where knowing the past outcomes
//	dont help to predict the next result.  b)
//	Taking parts of a LCG generated number can be
//	quite non-random: for example, looking at only
//	the least significant byte gives a permuted
//	8-bit counter (that has a period length of only
//	256).  The advantage of an LCA is that it is
//	perfectly uniform when run for the entire period
//	length (and very uniform for smaller sequences
//	too, if the parameters are chosen carefully).
//
//	ACGs have extremly long period lengths and
//	provide good independence.  Unfortunately,
//	uniformity isnt not too great. Furthermore, I
//	didnt find any theoretically analysis of ACGs
//	that addresses uniformity.
//
//	The RNG given below will return numbers
//	generated by an LCA that are permuted under
//	control of a ACG. 2 permutations take place: the
//	4 bytes of one LCG generated number are
//	subjected to one of 16 permutations selected by
//	4 bits of the ACG. The permutation a such that
//	byte of the result may come from each byte of
//	the LCG number. This effectively destroys the
//	structure within a word. Finally, the sequence
//	of such numbers is permuted within a range of
//	256 numbers. This greatly improves independence.
//
//
//  Algorithm M as describes in Knuths "Art of Computer Programming",
//	Vol 2. 1969
//  is used with a linear congruential generator (to get a good uniform
//  distribution) that is permuted with a Fibonacci additive congruential
//  generator to get good independence.
//
//  Bit, byte, and word distributions were extensively tested and pass
//  Chi-squared test near perfect scores (>7E8 numbers tested, Uniformity
//  assumption holds with probability > 0.999)
//
//  Run-up tests for on 7E8 numbers confirm independence with
//  probability > 0.97.
//
//  Plotting random points in 2d reveals no apparent structure.
//
//  Autocorrelation on sequences of 5E5 numbers (A(i) = SUM X(n)*X(n-i),
//	i=1..512)
//  results in no obvious structure (A(i) ~ const).
//
//  Except for speed and memory requirements, this generator outperforms
//  random() for all tests. (random() scored rather low on uniformity tests,
//  while independence test differences were less dramatic).
//
//  AGN would like to..
//  thanks to M.Mauldin, H.Walker, J.Saxe and M.Molloy for inspiration & help.
//
//  And I would (DGC) would like to thank Donald Kunth for AGN for letting me
//  use his extensions in this implementation.
//

//
//	Part of the table on page 28 of Knuth, vol II. This allows us
//	to adjust the size of the table at the expense of shorter sequences.
//

static int randomStateTable[][3] = {
{3,7,16}, {4,9, 32}, {3,10, 32}, {1,11, 32}, {1,15,64}, {3,17,128},
{7,18,128}, {3,20,128}, {2,21, 128}, {1,22, 128}, {5,23, 128}, {3,25, 128},
{2,29, 128}, {3,31, 128}, {13,33, 256}, {2,35, 256}, {11,36, 256},
{14,39,256}, {3,41,256}, {9,49,256}, {3,52,256}, {24,55,256}, {7,57, 256},
{19,58,256}, {38,89,512}, {17,95,512}, {6,97,512}, {11,98,512}, {-1,-1,-1} };

//
// spatial permutation table
//	RANDOM_PERM_SIZE must be a power of two
//

#define RANDOM_PERM_SIZE 64
uint32_t randomPermutations[RANDOM_PERM_SIZE] = {
0xffffffff, 0x00000000,  0x00000000,  0x00000000,  // 3210
0x0000ffff, 0x00ff0000,  0x00000000,  0xff000000,  // 2310
0xff0000ff, 0x0000ff00,  0x00000000,  0x00ff0000,  // 3120
0x00ff00ff, 0x00000000,  0xff00ff00,  0x00000000,  // 1230

0xffff0000, 0x000000ff,  0x00000000,  0x0000ff00,  // 3201
0x00000000, 0x00ff00ff,  0x00000000,  0xff00ff00,  // 2301
0xff000000, 0x00000000,  0x000000ff,  0x00ffff00,  // 3102
0x00000000, 0x00000000,  0x00000000,  0xffffffff,  // 2103

0xff00ff00, 0x00000000,  0x00ff00ff,  0x00000000,  // 3012
0x0000ff00, 0x00000000,  0x00ff0000,  0xff0000ff,  // 2013
0x00000000, 0x00000000,  0xffffffff,  0x00000000,  // 1032
0x00000000, 0x0000ff00,  0xffff0000,  0x000000ff,  // 1023

0x00000000, 0xffffffff,  0x00000000,  0x00000000,  // 0321
0x00ffff00, 0xff000000,  0x00000000,  0x000000ff,  // 0213
0x00000000, 0xff000000,  0x0000ffff,  0x00ff0000,  // 0132
0x00000000, 0xff00ff00,  0x00000000,  0x00ff00ff   // 0123
};

//
//	SEED_TABLE_SIZE must be a power of 2
//
#define SEED_TABLE_SIZE 32
static uint32_t seedTable[SEED_TABLE_SIZE] = {
0xbdcc47e5, 0x54aea45d, 0xec0df859, 0xda84637b,
0xc8c6cb4f, 0x35574b01, 0x28260b7d, 0x0d07fdbf,
0x9faaeeb0, 0x613dd169, 0x5ce2d818, 0x85b9e706,
0xab2469db, 0xda02b0dc, 0x45c60d6e, 0xffe49d10,
0x7224fea3, 0xf9684fc9, 0xfc7ee074, 0x326ce92a,
0x366d13b5, 0x17aaa731, 0xeb83a675, 0x7781cb32,
0x4ec7c92d, 0x7f187521, 0x2cf346b4, 0xad13310f,
0xb89cff2b, 0x12164de1, 0xa865168d, 0x32b56cdf
};

//
//	The LCG used to scramble the ACG
//
//
// LC-parameter selection follows recommendations in
// "Handbook of Mathematical Functions" by Abramowitz & Stegun 10th, edi.
//
// LC_A = 251^2, ~= sqrt(2^32) = 66049
// LC_C = result of a long trial & error series = 3907864577
//

static const uint32_t LC_A = 66049;
static const uint32_t LC_C = 3907864577U;
static inline uint32_t LCG(uint32_t x)
{
    return( x * LC_A + LC_C );
}


ACG::ACG(uint32_t seed, int size)
{
    int l;
    initialSeed = seed;

    //
    //	Determine the size of the state table
    //

    for (l = 0;
	 randomStateTable[l][0] != -1 && randomStateTable[l][1] < size;
	 l++);

    if (randomStateTable[l][1] == -1) {
	l--;
    }

    initialTableEntry = l;

    stateSize = randomStateTable[ initialTableEntry ][ 1 ];
    auxSize = randomStateTable[ initialTableEntry ][ 2 ];

    //
    //	Allocate the state table & the auxillary table in a single malloc
    //

    state = new uint32_t[stateSize + auxSize];
    auxState = &state[stateSize];

    reset();
}

//
//	Initialize the state
//
void
ACG::reset()
{
    uint32_t u;

    if (initialSeed < SEED_TABLE_SIZE) {
	u = seedTable[ initialSeed ];
    } else {
	u = initialSeed ^ seedTable[ initialSeed & (SEED_TABLE_SIZE-1) ];
    }


    j = randomStateTable[ initialTableEntry ][ 0 ] - 1;
    k = randomStateTable[ initialTableEntry ][ 1 ] - 1;

    int i;
    for(i = 0; i < stateSize; i++) {
	state[i] = u = LCG(u);
    }

    for (i = 0; i < auxSize; i++) {
	auxState[i] = u = LCG(u);
    }

    k = u % stateSize;
    int tailBehind = (stateSize - randomStateTable[ initialTableEntry ][ 0 ]);
    j = k - tailBehind;
    if (j < 0) {
	j += stateSize;
    }

    lcgRecurr = u;

    assert(sizeof(double) == 2 * sizeof(int32_t));
}

ACG::~ACG()
{
    if (state) delete [] state;
    state = 0;
    // don't delete auxState, it's really an alias for state.
}

//
//	Returns 32 bits of random information.
//

uint32_t
ACG::asLong()
{
    uint32_t result = state[k] + state[j];
    state[k] = result;
    j = (j <= 0) ? (stateSize-1) : (j-1);
    k = (k <= 0) ? (stateSize-1) : (k-1);

    short int auxIndex = (result >> 24) & (auxSize - 1);
    uint32_t auxACG = auxState[auxIndex];
    auxState[auxIndex] = lcgRecurr = LCG(lcgRecurr);

    //
    // 3c is a magic number. We are doing four masks here, so we
    // do not want to run off the end of the permutation table.
    // This insures that we have always got four entries left.
    //
    uint32_t *perm = & randomPermutations[result & 0x3c];

    result =  *(perm++) & auxACG;
    result |= *(perm++) & ((auxACG << 24)
			   | ((auxACG >> 8)& 0xffffff));
    result |= *(perm++) & ((auxACG << 16)
			   | ((auxACG >> 16) & 0xffff));
    result |= *(perm++) & ((auxACG <<  8)
			   | ((auxACG >> 24) &   0xff));

    return(result);
}


