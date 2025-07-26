
export module oml.unop;

export {
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


}