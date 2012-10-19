// File: Ran250.C

// Modifications for C++ and oml containers Copyright (1994-2004), Jan N. Reimers

#include "oml/ran250.h"
#include "oml/minmax.h"
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
    register int l;
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
    register uint32_t u;

    if (initialSeed < SEED_TABLE_SIZE) {
	u = seedTable[ initialSeed ];
    } else {
	u = initialSeed ^ seedTable[ initialSeed & (SEED_TABLE_SIZE-1) ];
    }


    j = randomStateTable[ initialTableEntry ][ 0 ] - 1;
    k = randomStateTable[ initialTableEntry ][ 1 ] - 1;

    register int i;
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
    register uint32_t auxACG = auxState[auxIndex];
    auxState[auxIndex] = lcgRecurr = LCG(lcgRecurr);

    //
    // 3c is a magic number. We are doing four masks here, so we
    // do not want to run off the end of the permutation table.
    // This insures that we have always got four entries left.
    //
    register uint32_t *perm = & randomPermutations[result & 0x3c];

    result =  *(perm++) & auxACG;
    result |= *(perm++) & ((auxACG << 24)
			   | ((auxACG >> 8)& 0xffffff));
    result |= *(perm++) & ((auxACG << 16)
			   | ((auxACG >> 16) & 0xffff));
    result |= *(perm++) & ((auxACG <<  8)
			   | ((auxACG >> 24) &   0xff));

    return(result);
}
